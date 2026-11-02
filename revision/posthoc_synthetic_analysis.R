# posthoc_synthetic_analysis.R
# Analyses of an existing synthetic history matching run that do NOT change the
# reported results: they only evaluate the saved emulators and data.
#   1. truth check: is the true parameter vector ruled out at any wave? (R2.4)
#   2. per-wave table: design sizes, bounding-box fraction, estimated proportion
#      of the parameter space retained, emulator variance ratio (R2.4, stopping rule)
#   3. emulator validation metrics per target and wave (R2.4)
#   4. correlation length sensitivity, theta = 0.4 / 0.6 / 0.8 (R2.4)
#
# Input: the list saved at the end of the synthetic script,
#   this_plots_data = list(emulators = ems_list, data = data_list, priors = ...)
# ems_list has the newest wave first; data_list[[w]] holds the runs used to
# build the Wave w emulators.
#
# Needs the same set-up as the synthetic script (hmer, CARDAMOM functions,
# vars_DALEC.A3.C3.H2.M1, edc1.R for assess_EDC1, prior_vals_list).

library(hmer); library(dplyr); library(purrr); library(tibble)

# Location of the scratch folder: "~/scratch" on the cluster, e.g. "Z:/scratch" on Windows.
if (!exists("scratch")) scratch <- "~/scratch"

# ---- set-up (same as the synthetic script) ----------------------------------------
# All code is loaded from this folder (relative to the CARDAMOM working directory).
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
source(paste0(code_dir, "parameter_ranges.R"))
source(paste0(code_dir, "edc1.R"))
source("./R_functions/load_all_cardamom_functions.r")
load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year1/infofile.RData")
prior_vals_list <- list(sow_day = 274, harvest_length = 342)

# ---- choose the run -----------------------------------------------------------------
# A file saved by run_synthetic_all.R, e.g. exp1a_md10_seed1.rds, edcs_md20_seed1.rds
if (!exists("run_file")) run_file <- paste0(scratch, "/synthetic_revision/exp1a_md10_seed1.rds")

res        <- readRDS(run_file)
discs_prop <- res$settings$discs_prop
synth     <- readRDS(paste0(scratch, "/synth_obs_feb.RDS"))
true_pars <- synth$true_pars
targets   <- res$settings$targets          # targets used in this run (no EDC targets for Exp. 2)

ems_list  <- res$emulators
data_list <- res$data
n_waves   <- length(ems_list)
ems_wave  <- function(w) ems_list[[n_waves - w + 1]]          # emulators of wave w
original_ranges <- ems_wave(1)[[1]]$ranges

# ---- helpers -------------------------------------------------------------------
disc_var <- function(em) sum(map_dbl(em$disc, ~ ifelse(is.null(.x), 0, .x^2)))

# TRUE if a point is not implausible under the emulators of waves 1..w
# (same rule as generate_new_design_edcs: I <= cutoff, one-sided for DS_max)
not_implausible <- function(points, w, cutoff = 3) {
  keep <- rep(TRUE, nrow(points))
  for (j in seq_len(w)) {
    for (em in ems_wave(j)) {
      x <- points[, names(em$ranges), drop = FALSE]
      if (em$output_name == "DS_max") {
        sd_tot <- sqrt(em$get_cov(x) + targets$DS_max$sigma^2 + disc_var(em))
        keep <- keep & (targets$DS_max$val <= em$get_exp(x) + cutoff * sd_tot)
      } else {
        keep <- keep & em$implausibility(x, targets[[em$output_name]], cutoff)
      }
    }
  }
  keep
}

uniform_in_box <- function(n, ranges) {
  as.data.frame(purrr::map(ranges, ~ runif(n, .x[1], .x[2])))
}

box_fraction <- function(ranges) {
  prod(map_dbl(names(original_ranges), ~ diff(ranges[[.x]]) / diff(original_ranges[[.x]])))
}

# ---- 1. truth check ----------------------------------------------------------------
x_true <- as.data.frame(as.list(true_pars))
truth_check <- map_dfr(seq_len(n_waves), function(w) {
  rng <- ems_wave(w)[[1]]$ranges
  in_box <- all(map_lgl(names(rng), ~ x_true[[.x]] >= rng[[.x]][1] & x_true[[.x]] <= rng[[.x]][2]))
  imp <- map_dbl(ems_wave(w)[names(ems_wave(w)) != "DS_max"],
                 ~ .x$implausibility(x_true[, names(.x$ranges), drop = FALSE],
                                     targets[[.x$output_name]]))
  i_star <- NA_real_                       # no DS_max emulator in Exp. 2
  if ("DS_max" %in% names(ems_wave(w))) {
    em <- ems_wave(w)[["DS_max"]]
    xd <- x_true[, names(em$ranges), drop = FALSE]
    i_star <- (targets$DS_max$val - em$get_exp(xd)) /
      sqrt(em$get_cov(xd) + targets$DS_max$sigma^2 + disc_var(em))
  }
  tibble(wave = w, truth_in_ranges = in_box, I_max = max(imp),
         worst_target = names(imp)[which.max(imp)], I_star_DS_max = i_star,
         truth_retained = in_box & max(imp) <= 3 & (is.na(i_star) | i_star <= 3))
})
print(truth_check)

# ---- 2. per-wave table ---------------------------------------------------------------
# The reruns log this table during the run (res$wave_table), including the share of the
# parameter space retained. The Monte Carlo estimate below is slow (it evaluates all
# emulators at 2e5 points per wave) and only needed for runs without a logged table;
# set estimate_volume <- TRUE before sourcing to use it anyway.
if (!exists("estimate_volume")) estimate_volume <- is.null(res$wave_table)
if (!estimate_volume) {
  wave_table <- res$wave_table
} else {
  # proportion retained = (volume of the wave's bounding box) x (fraction of uniform
  # points in that box passing the input-EDCs and the emulators of waves 1..w).
  # The output-EDC check on DALEC-Crop runs is not included, so this slightly
  # overestimates the NROY volume.
  set.seed(1)
  n_mc <- 2e5                      # increase if a wave accepts very few points
  wave_table <- map_dfr(seq_len(n_waves), function(w) {
    rng <- ems_wave(w)[[1]]$ranges
    pts <- uniform_in_box(n_mc, rng)
    pts_edc <- if (isTRUE(res$settings$use_edcs)) assess_EDC1(pts, fixed_pars = prior_vals_list) else pts
    pass <- if (nrow(pts_edc) > 0) sum(unlist(purrr::map(split(seq_len(nrow(pts_edc)),
                                                    ceiling(seq_len(nrow(pts_edc)) / 2e4)),
                                                ~ not_implausible(pts_edc[.x, , drop = FALSE], w))))
            else 0
    ratios <- map_dbl(ems_wave(w), ~ mean(.x$get_cov(data_list[[w]]$testing_data[, names(.x$ranges)])) /
                        (targets[[.x$output_name]]$sigma^2 + disc_var(.x)))
    tibble(wave = w,
           n_train = nrow(data_list[[w]]$training_data),
           n_test = nrow(data_list[[w]]$testing_data),
           box_fraction = box_fraction(rng),
           share_passing_in_box = pass / n_mc,
           proportion_retained = box_fraction(rng) * pass / n_mc,
           max_variance_ratio = max(ratios),
           worst_target = names(ratios)[which.max(ratios)])
  })
  print(wave_table)
}

# ---- 3. emulator validation metrics ---------------------------------------------------
validation <- map_dfr(seq_len(n_waves), function(w) {
  test <- data_list[[w]]$testing_data
  map_dfr(ems_wave(w), function(em) {
    x <- test[, names(em$ranges), drop = FALSE]
    e <- em$get_exp(x); v <- em$get_cov(x); f <- test[[em$output_name]]
    se <- (f - e) / sqrt(v)
    tibble(wave = w, output = em$output_name, rmse = sqrt(mean((f - e)^2)),
           mean_abs_std_error = mean(abs(se)), max_abs_std_error = max(abs(se)),
           coverage_3sd = mean(abs(se) <= 3))
  })
})
print(validation, n = 60)

# ---- 4. theta sensitivity (refits Wave 1 emulators only) ---------------------------------
output_names <- names(ems_wave(1))
model_discs <- purrr::map(output_names, ~ list(internal = discs_prop * targets[[.x]]$val,
                                         external = NULL))   # same (fixed) construction as hm_runner.R
theta_check <- map_dfr(c(0.4, 0.6, 0.8), function(th) {
  ems <- emulator_from_data(
    input_data = data_list[[1]]$training_data, output_names = output_names,
    specified_priors = list(hyper_p = rep(th, length(output_names))),
    ranges = original_ranges, beta.var = TRUE, check.ranges = FALSE,
    targets = targets[output_names], add = TRUE, order = 2, discrepancies = model_discs)
  test <- data_list[[1]]$testing_data
  map_dfr(ems, function(em) {
    x <- test[, names(em$ranges), drop = FALSE]
    se <- (test[[em$output_name]] - em$get_exp(x)) / sqrt(em$get_cov(x))
    tibble(theta = th, output = em$output_name,
           mean_abs_std_error = mean(abs(se)), coverage_3sd = mean(abs(se) <= 3))
  })
})
print(theta_check %>% group_by(theta) %>%
        summarise(mean_abs_std_error = mean(mean_abs_std_error),
                  min_coverage_3sd = min(coverage_3sd), .groups = "drop"))

# ---- save ------------------------------------------------------------------------------
saveRDS(list(run_file = run_file, discs_prop = discs_prop, truth_check = truth_check,
             wave_table = wave_table, validation = validation, theta_check = theta_check),
        sub("\\.rds$", "_posthoc.rds", run_file))
