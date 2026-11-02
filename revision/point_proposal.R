generate_new_design_edcs = function(all_ems, n_requested_points, targets,
                                    cutoff = 3, n = 1, use_edcs = TRUE, mc_cores = 10){
  # REVISION: identical proposal logic to the original version, plus
  #   - counts of proposed / accepted points, returned as attr(result, "acceptance")
  #   - use_edcs = FALSE skips the input- and output-EDCs (Exp. 2)
  #   - no DS_max step when there is no DS_max emulator (Exp. 2)
  #   - mc_cores argument (was fixed at 10)
  # n = 1: max implausibility. n = 2: second largest implausibility
  wave_number = length(all_ems)
  cat("We are at wave:", wave_number, "!\n")

  # get newest emulator's ranges
  # assume that newest emulator is the first
  current_ranges = all_ems[[1]][[1]]$ranges

  # separate between normal emulators and DS_max
  all_ems = unlist(all_ems) # first flatten the list
  ems_normal = all_ems[!names(all_ems) %in% c("DS_max")]
  ems_DS = all_ems[names(all_ems) %in% c("DS_max")]

  if (length(cutoff) == 1) cutoff <- rep(cutoff, length(ems_normal))  # normal emulators

  counts = c(rounds = 0, proposed = 0, pass_edc1 = 0, pass_impl = 0, accepted = 0)

  # now enter point proposal loop
  NROY = data.frame()
  count = 1
  points.factor = 200
  while(nrow(NROY)<n_requested_points){
    if(count > 50){
      stop("Too many iterations, we might deduce that there are no good matches")
    }
    cat("Count:", count, "\n")
    counts[["rounds"]] = counts[["rounds"]] + 1
    # first generate LHD scaled between [-1,1]
    LHD = setNames(
      data.frame(2* (lhs::randomLHS(n_requested_points * points.factor, length(current_ranges)) - 0.5)),
      names(current_ranges)
    )
    counts[["proposed"]] = counts[["proposed"]] + nrow(LHD)
    # Scale each column to the original parameter range
    scaled_lhs <- as.data.frame(
      mapply(
        function(values, par_min, par_max) par_min + (par_max - par_min) * ((values + 1) / 2),
        LHD,
        unlist(lapply(current_ranges, function(x) x[1])),
        unlist(lapply(current_ranges, function(x) x[2]))
      )
    )
    names(scaled_lhs) <- names(current_ranges)

    if (use_edcs) scaled_lhs = assess_EDC1(scaled_lhs, fixed_pars = prior_vals_list)
    counts[["pass_edc1"]] = counts[["pass_edc1"]] + nrow(scaled_lhs)

    # now get the points' implausibility
    cat("now calculating implausibilities for", nrow(scaled_lhs), "points... \n")
    implausibles <- mclapply(seq_along(ems_normal), function(i) {
      current_model <- ems_normal[[i]]
      implaus_fn <- current_model[["implausibility"]]
      output_name <- current_model[["output_name"]]
      target_z <- targets[[output_name]]
      current_cutoff <- cutoff[[i]]
      implaus_fn(scaled_lhs, target_z, current_cutoff) # returns true if point is non implausible
    }, mc.cores = mc_cores)
    # a worker that crashed (usually out of memory) returns an error or NULL instead of a
    # logical vector: recompute those emulators one at a time instead of stopping
    failed <- which(sapply(implausibles, function(x) !is.logical(x) || length(x) != nrow(scaled_lhs)))
    if (length(failed) > 0) {
      message(length(failed), " parallel job(s) failed (memory?); recomputing them serially")
      for (i in failed) {
        current_model <- ems_normal[[i]]
        implausibles[[i]] <- current_model[["implausibility"]](
          scaled_lhs, targets[[current_model[["output_name"]]]], cutoff[[i]])
      }
    }

    d_implausibles <- setNames(
      data.frame(implausibles),
      map_chr(ems_normal, ~.$output_name))
    # an output is not implausible only if it is not implausible under every wave's emulator
    imp_mat <- t(apply(
      d_implausibles, 1,
      function(x) map_lgl(unique(names(x)), ~all(x[names(x) == .]))))
    nonimpl_vector = as.vector(apply(imp_mat, 1, function(x) sum(x) > length(x)-n))
    scaled_lhs = scaled_lhs[nonimpl_vector,]
    cat("we have this many points after implausibility:", nrow(scaled_lhs), "out of", n_requested_points * points.factor, "\n")

    if(nrow(scaled_lhs) > 0 && length(ems_DS) > 0){
      # one-sided implausibility for DS_max
      implausibles_DS <- mclapply(seq_along(ems_DS), function(i) {
        DS_max_implausibility(scaled_lhs, single_ems_DS_flat = ems_DS[[i]], targets, cutoff = 3)
      }, mc.cores = length(ems_DS))
      d_ds_implausibles <- do.call("cbind", implausibles_DS)
      nonimpl_vector_DS = as.vector(apply(d_ds_implausibles, 1, function(x) sum(x) > length(x)-1))
      scaled_lhs = scaled_lhs[nonimpl_vector_DS,]
    }
    counts[["pass_impl"]] = counts[["pass_impl"]] + nrow(scaled_lhs)

    # output-EDCs on the actual DALEC-Crop output
    if(nrow(scaled_lhs) > 1){
      if (use_edcs) scaled_lhs = assess_edcs2(scaled_lhs)
      counts[["accepted"]] = counts[["accepted"]] + nrow(scaled_lhs)
      NROY = rbind(NROY, scaled_lhs)
      cat("Collected ", nrow(NROY), "NROY points. \n")
    }

    count = count + 1
  }

  # once we have enough points, select subset using maximin:
  print("choosing final points using maximin sampling.")
  final_LHD = hmer:::maximin_sample(NROY, n = n_requested_points, nms = names(scaled_lhs))
  attr(final_LHD, "acceptance") = counts
  return(final_LHD)
}

DS_max_implausibility = function(points, single_ems_DS_flat, targets, cutoff = 3){
  # one-sided implausibility (Equation 9): not implausible if z <= E_D[f(x)] + cutoff * SD
  DS_target = targets[["DS_max"]]

  disc_quad <- sum(map_dbl(single_ems_DS_flat$disc, ~.^2)) # model discrepancy
  imp_var <- single_ems_DS_flat$get_cov(points) + DS_target$sigma^2 + disc_quad
  imp <- single_ems_DS_flat$get_exp(points) + cutoff*sqrt(imp_var)

  return(DS_target$val<=imp) # return true if point is not implausible
}
