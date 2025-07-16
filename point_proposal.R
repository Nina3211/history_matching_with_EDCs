generate_new_design_edcs = function(all_ems, n_requested_points, targets, 
                                    cutoff = 3, n = 1){
  # n = 1: max implausibility. n = 2: second largest implausibility
  # find out which wave we are at:
  wave_number = length(all_ems)
  cat("We are at wave:", wave_number, "!\n")
  
  # get newest emulator's ranges
  # assume that newest emulator is the first
  current_ranges = all_ems[[1]][[1]]$ranges

  # separate between normal emulators and DS_max
  all_ems = unlist(all_ems) # first flatten the list
  ems_normal = all_ems[!names(all_ems) %in% c("DS_max")]
  ems_DS = all_ems[names(all_ems) %in% c("DS_max")]
  cat("So many normal emulators:", length(ems_normal), "\n")
  cat("So many DS_max emulators:", length(ems_DS), "\n")
  
  
  if (length(cutoff) == 1) cutoff <- rep(cutoff, length(ems_normal))  # normal emulators
  
  # now enter point proposal loop
  NROY = data.frame()
  count = 1
  points.factor = 200
  while(nrow(NROY)<n_requested_points){
    if(count > 50){
      stop("Too many iterations, we might deduce that there are no good matches")
    }
    cat("Count:", count, "\n")
    # first generate LHD scaled between [-1,1]
    LHD = setNames(
      data.frame(2* (lhs::randomLHS(n_requested_points * points.factor, length(current_ranges)) - 0.5)),
      names(current_ranges)
    )
    # Scale each column to the original parameter range
    scaled_lhs <- as.data.frame(
      mapply(
        function(values, par_min, par_max) par_min + (par_max - par_min) * ((values + 1) / 2),
        LHD,
        unlist(lapply(current_ranges, function(x) x[1])),
        unlist(lapply(current_ranges, function(x) x[2]))
      )
    )

    # Ensure column names are preserved
    names(scaled_lhs) <- names(current_ranges)
    
    scaled_lhs = assess_EDC1(scaled_lhs, fixed_pars = prior_vals_list)
    
    # now get the points' implausibility
    cat("now calculating implausibilities for", nrow(scaled_lhs), "points... \n")
    implausibles <- mclapply(seq_along(ems_normal), function(i) {
      cat(i, "parallel \n")
      current_model <- ems_normal[[i]]
      implaus_fn <- current_model[["implausibility"]]
      output_name <- current_model[["output_name"]]
      target_z <- targets[[output_name]]
      current_cutoff <- cutoff[[i]]

      # Execute the implausibility function and expect a logical outcome
      implaus_fn(scaled_lhs, target_z, current_cutoff) # returns true if point is non implausible
    }, mc.cores = 10)  # this returns a list of vectors with true/ false depending on the implausibility
    
    # putting this into a dataframe
    d_implausibles <- setNames(
      data.frame(implausibles),
      map_chr(ems_normal, ~.$output_name))
    #  For each unique output name, check if all values in that row that match the current output name are TRUE
    # if there is one false, then we can already collapse that output to false, i.e. implausible
    imp_mat <- t(apply(
      d_implausibles, 1,
      function(x) map_lgl(unique(names(x)), ~all(x[names(x) == .]))))
    # calculate true if there are more trues (not implausibles) than max number of trues - n
    # i.e. check if 2nd most implausible point (value doesn't matter) is larger than three
    # i.e. we consider this point to be plausible
    nonimpl_vector = as.vector(apply(imp_mat, 1, function(x) sum(x) > length(x)-n))
    # subset proposed points by implausibility, only keep non implausible points
    scaled_lhs = scaled_lhs[nonimpl_vector,]
    cat("we have this many points after implausibility:", nrow(scaled_lhs), "out of", n_requested_points * points.factor, "\n")
    
    if(nrow(scaled_lhs)> 0){
    # now apply calculate implausibility for DS_max
    implausibles_DS <- mclapply(seq_along(ems_DS), function(i) {
      cat(i, "parallel DS \n")
      current_model <- ems_DS[[i]]
      # Execute the implausibility function and expect a logical outcome
      DS_max_implausibility(scaled_lhs, single_ems_DS_flat = current_model, targets, cutoff = 3) # returns true if point is non implausible
    }, mc.cores = length(ems_DS))  # this returns a list of vectors with true/ false depending on the implausibility
    d_ds_implausibles <- do.call("cbind", implausibles_DS)
    nonimpl_vector_DS = as.vector(apply(d_ds_implausibles, 1, function(x) sum(x) > length(x)-1))
    
    # subset the LHD again according to DS:
    scaled_lhs = scaled_lhs[nonimpl_vector_DS,]
    }
    # now run these proposed points through the actual assess_EDC2 function
    # only do this if this round actually gave us NROY points
    if(nrow(scaled_lhs)> 1){
      scaled_lhs = assess_edcs2(scaled_lhs)
      NROY = rbind(NROY, scaled_lhs)
      cat("Collected ", nrow(NROY), "NROY points. \n")
    }
    


    count = count + 1
    
  }
  
  # once we have enough points, select subset using maximin: 
  print("choosing final points using maximin sampling.")
  final_LHD = hmer:::maximin_sample(NROY, n = n_requested_points, nms = names(scaled_lhs))
  return(final_LHD)
}

DS_max_implausibility = function(points, single_ems_DS_flat, targets, cutoff = 3){
  # implausibility is (E_D[f(x)] - z)/(sqrt(Var_D[f(x)] + Var_md + Var_obs))
  # extract DS_max target
  DS_target = targets[["DS_max"]]
  
  disc_quad <- sum(map_dbl(single_ems_DS_flat$disc, ~.^2)) # model discrepancy
  imp_var <- single_ems_DS_flat$get_cov(points) + DS_target$sigma^2 + disc_quad
  # key difference to usual implausibility is here: instead of taking the abs value, take this
  imp <- single_ems_DS_flat$get_exp(points) + cutoff*sqrt(imp_var)
  
  return(DS_target$val<=imp) # return true if point is not implausible

}
  