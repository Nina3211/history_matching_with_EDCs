# hm_runner.R
# run_history_matching() reproduces the workflow of the paper's scripts
# (synthetic example, ATEC_year1/2_hm_ranges.R) as one function:
#   Wave 0: 720 training + 20 test points (initial_LHD with EDCs)
#   Waves 1-3: emulators -> 360 training + 20 test points -> DALEC-Crop -> new ranges
#   Wave 4: emulators -> 300 final NROY points -> DALEC-Crop (with LAI curves)
# Same emulator settings (order 2, theta = 0.6, beta.var, add = TRUE), same targets,
# same point proposal and maximin selection as before.
#
# Changes compared with the paper's runs (3. also drops runs with non-finite outputs,
# which only occurs without EDCs):
#   1. BUG FIX: discrepancies are matched to the right outputs. hmer assigns them by
#      position in output_names; before, they followed the order of `targets`, so the
#      Ra:GPP and yield:GPP targets got no discrepancy and the SOM targets got theirs.
#   2. BUG FIX (in edc2.R): EDC_SOM1 is now checked on the simulator output.
#   3. Logging only: emulators of all waves, all wave data, per-wave table, emulator
#      validation metrics, DALEC-Crop evaluation counts and time per stage are saved.
#
# Needs: hmer, dplyr, purrr, tibble, parallel, lhs; the CARDAMOM R functions;
# PROJECT / new_PROJECT / drivers; ranges_df_with_priors; and the sourced
# edc1.R, edc2.R, build_initial_LHD.R, point_proposal.R, generate_input_output_DALEC_Crop.R.

library(hmer); library(dplyr); library(purrr); library(tibble); library(parallel)

# ---- counter for DALEC-Crop evaluations (also incremented in edc2.R) -----------------
.hm_log <- new.env()
reset_dalec_counter <- function() assign("n_dalec", 0, envir = .hm_log)
count_dalec <- function(n) assign("n_dalec", get("n_dalec", envir = .hm_log) + n, envir = .hm_log)
reset_dalec_counter()

run_dalec <- function(inputs, ...) {
  count_dalec(nrow(inputs))
  drop_nonfinite(generate_input_output_DALEC_Crop(inputs, ...))
}

# Remove runs whose outputs are infinite or missing. With EDCs these points are already
# removed by the output-EDC check (edc2.R), so this only matters without EDCs (Exp. 2),
# where they would otherwise break the emulator fit. Handles both return formats of
# generate_input_output_DALEC_Crop (data frame, or list with LAI_curves in the final wave).
output_cols <- c("harvest_vector", "EDC_SOM1", "EDC_SOM2", "ra_GPP_ratio",
                 "yield_GPP_ratio", "DS_max")
drop_nonfinite <- function(out) {
  td <- if (is.data.frame(out)) out else out$training_data
  cols <- intersect(c(grep("^date", names(td), value = TRUE), output_cols), names(td))
  keep <- apply(is.finite(as.matrix(td[, cols])), 1, all)
  if (all(keep)) return(out)
  message("Dropping ", sum(!keep), " of ", length(keep), " runs with non-finite outputs")
  if (is.data.frame(out)) return(out[keep, , drop = FALSE])
  lai_runs <- setdiff(names(out$LAI_curves), "date")       # one column per run, same order
  out$training_data <- td[keep, , drop = FALSE]
  out$LAI_curves <- out$LAI_curves[, c(lai_runs[keep], "date")]
  out
}

# ---- saving ------------------------------------------------------------------------
# Writes to a temporary file and renames it only after a successful save, so a full disk
# never leaves a truncated .rds behind (which the drivers would otherwise skip as "done").
safe_save <- function(obj, file) {
  tmp <- paste0(file, ".tmp")
  ok <- tryCatch({ saveRDS(obj, tmp); file.rename(tmp, file) },
                 error = function(e) { message("Saving failed: ", conditionMessage(e)); FALSE })
  if (!isTRUE(ok)) {
    unlink(tmp)
    stop("Could not save ", basename(file), ". Is the disk or your quota full? ",
         "Free some space and rerun: finished runs are skipped.")
  }
  invisible(file)
}

# ---- helpers ------------------------------------------------------------------------
# ranges of the next wave: min/max of the new training points, widened by `buffer`
# (0.025 in the synthetic script, 0.01 in the ATEC scripts)
find_ranges_fct <- function(input_data, original_ranges, buffer) {
  input_data <- input_data[names(original_ranges)]
  ranges <- setNames(purrr::map(names(original_ranges), function(nm) {
    c(max(original_ranges[[nm]][1], min(input_data[, nm]) - buffer * diff(range(input_data[, nm]))),
      min(original_ranges[[nm]][2], max(input_data[, nm]) + buffer * diff(range(input_data[, nm]))))
  }), names(original_ranges))
  hmer:::convertRanges(ranges)
}

box_fraction <- function(ranges, original_ranges) {
  prod(map_dbl(names(original_ranges), ~ diff(ranges[[.x]]) / diff(original_ranges[[.x]])))
}

disc_var <- function(em) sum(map_dbl(em$disc, ~ ifelse(is.null(.x), 0, .x^2)))

validation_metrics <- function(ems, test_data) {
  map_dfr(ems, function(em) {
    x  <- test_data[, names(em$ranges), drop = FALSE]
    e  <- em$get_exp(x); v <- em$get_cov(x); f <- test_data[[em$output_name]]
    se <- (f - e) / sqrt(v)
    tibble(output = em$output_name, n_test = length(f), rmse = sqrt(mean((f - e)^2)),
           mean_abs_std_error = mean(abs(se)), max_abs_std_error = max(abs(se)),
           coverage_3sd = mean(abs(se) <= 3))
  })
}

# mean emulator variance / (observation error variance + discrepancy variance)
variance_ratio <- function(ems, points, targets) {
  map_dbl(ems, ~ mean(.x$get_cov(points[, names(.x$ranges), drop = FALSE])) /
            (targets[[.x$output_name]]$sigma^2 + disc_var(.x)))
}

# ---- main function ----------------------------------------------------------------------
run_history_matching <- function(targets, obs_dates, obs_dates_names, year,
                                 fixed_pars,              # e.g. list(sow_day = 274, harvest_length = 342)
                                 discs_prop = 0.1,
                                 use_edcs = TRUE,         # FALSE = Exp. 2 (hmer proposal, no EDCs)
                                 range_buffer = 0.01,     # 0.025 synthetic, 0.01 ATEC
                                 theta = 0.6,
                                 n_wave0 = 2 * 10 * 36, n_train = 10 * 36,
                                 n_test = 20, n_final = 300, n_waves = 4,
                                 seed = 1, mc_cores = 10) {
  set.seed(seed)
  reset_dalec_counter()
  timings <- list()
  timed <- function(stage, expr) {
    t0 <- Sys.time(); out <- expr
    timings[[stage]] <<- c(timings[[stage]], as.numeric(difftime(Sys.time(), t0, units = "secs")))
    out
  }
  t_start <- Sys.time()

  # helper functions in the utils read this global
  assign("prior_vals_list", fixed_pars, envir = .GlobalEnv)

  if (!use_edcs) targets <- targets[!names(targets) %in% c("EDC_SOM1", "EDC_SOM2", "DS_max")]
  dalec_args <- list(obs_dates = obs_dates, obs_dates_names = obs_dates_names,
                     fixed_pars = fixed_pars, year = year)

  # ---- Wave 0 ----
  wave0 <- timed("point_proposal", list(
    train = initial_LHD(n_wave0, use_edcs, other_fixed_vars = NULL, pars_with_priors = NULL,
                        fixed_vars = names(fixed_pars), ranges_df_with_priors = ranges_df_with_priors),
    test  = initial_LHD(n_test, use_edcs, other_fixed_vars = NULL, pars_with_priors = NULL,
                        fixed_vars = names(fixed_pars), ranges_df_with_priors = ranges_df_with_priors)))
  initial_ranges <- purrr::map(wave0$train$wave0_ranges, ~ as.numeric(unlist(.x)))

  data_list <- list(wave0 = timed("simulator_runs", list(
    training_data = do.call(run_dalec, c(list(wave0$train$initial_LHD), dalec_args)),
    testing_data  = do.call(run_dalec, c(list(wave0$test$initial_LHD),  dalec_args)))))

  output_names <- intersect(names(targets), colnames(data_list[[1]]$training_data))
  # BUG FIX: one discrepancy per output, in the order of output_names
  model_discs <- purrr::map(output_names, ~ list(internal = discs_prop * targets[[.x]]$val, external = NULL))

  ranges_list <- list(initial_ranges)
  ems_list <- list(); wave_rows <- list(); validation <- list()

  # point proposal: with EDCs, the paper's EDC proposal (point_proposal.R); without EDCs
  # (Exp. 2), hmer's standard generate_new_design(), as in the paper's Exp. 2 script.
  # hmer does not report acceptance counts, so those columns are NA for Exp. 2.
  propose <- function(n_points) {
    if (use_edcs) {
      generate_new_design_edcs(all_ems = ems_list, n_requested_points = n_points,
                               targets = targets, cutoff = 3, n = 1,
                               use_edcs = TRUE, mc_cores = mc_cores)
    } else {
      out <- generate_new_design(ems_list, n_points = n_points, z = targets,
                                 cutoff = 3, n = 1, verbose = FALSE)
      attr(out, "acceptance") <- c(rounds = NA, proposed = NA, pass_edc1 = NA,
                                   pass_impl = NA, accepted = NA)
      out
    }
  }

  # ---- Waves 1, ..., n_waves ----
  for (wave in seq_len(n_waves)) {
    final_wave <- wave == n_waves
    n_dalec_before <- get("n_dalec", envir = .hm_log)

    new_ems <- timed("emulator_fitting", emulator_from_data(
      input_data = data_list[[wave]]$training_data, output_names = output_names,
      specified_priors = list(hyper_p = rep(theta, length(output_names))),
      ranges = ranges_list[[wave]], beta.var = TRUE, check.ranges = FALSE,
      targets = targets[output_names], add = TRUE, order = 2,
      discrepancies = model_discs, verbose = FALSE))
    ems_list <- c(list(new_ems), ems_list)

    validation[[wave]] <- timed("diagnostics",
      validation_metrics(new_ems, data_list[[wave]]$testing_data) %>% mutate(wave = wave, .before = 1))

    new_train <- timed("point_proposal", propose(if (final_wave) n_final else n_train))
    acc <- attr(new_train, "acceptance")
    if (!final_wave) {
      new_test <- timed("point_proposal", propose(n_test))
      acc <- acc + attr(new_test, "acceptance")
    }

    vr <- variance_ratio(new_ems, new_train, targets)
    bf <- box_fraction(ranges_list[[wave]], initial_ranges)

    data_list[[wave + 1]] <- timed("simulator_runs", if (final_wave) {
      list(training_data = do.call(run_dalec, c(list(new_train), dalec_args, list(is_final_wave = TRUE))))
    } else {
      list(training_data = do.call(run_dalec, c(list(new_train), dalec_args)),
           testing_data  = do.call(run_dalec, c(list(new_test),  dalec_args)))
    })

    wave_rows[[wave]] <- tibble(
      wave = wave,
      n_train = nrow(data_list[[wave]]$training_data),
      n_test  = nrow(data_list[[wave]]$testing_data),
      proposal_rounds = acc[["rounds"]],
      points_proposed = acc[["proposed"]],
      pass_input_edcs = acc[["pass_edc1"]],
      pass_implausibility = acc[["pass_impl"]],
      pass_output_edcs = acc[["accepted"]],
      box_fraction = bf,
      proportion_retained = bf * acc[["accepted"]] / acc[["proposed"]],
      dalec_evaluations = get("n_dalec", envir = .hm_log) - n_dalec_before,
      max_variance_ratio = max(vr),
      worst_target = names(vr)[which.max(vr)])

    if (!final_wave) ranges_list[[wave + 1]] <- find_ranges_fct(new_train, ranges_list[[wave]], range_buffer)
  }

  timings$total <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  list(emulators = ems_list,        # newest wave first
       data = data_list,            # data[[1]] = Wave 0, ..., data[[n_waves + 1]] = final 300 points
       ranges = ranges_list,
       wave_table = bind_rows(wave_rows),
       validation = bind_rows(validation),
       timings = map_dbl(timings, sum),
       dalec_evaluations_total = get("n_dalec", envir = .hm_log),
       settings = list(targets = targets, discs_prop = discs_prop, use_edcs = use_edcs,
                       range_buffer = range_buffer, theta = theta, seed = seed,
                       fixed_pars = fixed_pars, mc_cores = mc_cores,
                       date = Sys.time(), session = sessionInfo()))
}
