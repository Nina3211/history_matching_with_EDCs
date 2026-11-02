# timing_benchmark.R
# Cost per evaluation: one DALEC-Crop run vs one emulator evaluation (R2.8).
# Uses the saved synthetic Exp. 1a run; nothing is calibrated. Single core throughout.
# Run on the cluster from the CARDAMOM folder:
#   source("synth_experiment/revision/timing_benchmark.R")

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
if (!exists("n_bench"))  n_bench  <- 1000          # number of parameter sets to time

# set-up as in run_synthetic_all.R (all runs exist, so nothing is rerun)
source(paste0(code_dir, "run_synthetic_all.R"))
library(dplyr)

res <- readRDS(file.path(scratch, "synthetic_revision", "exp1a_md10_seed1.rds"))
calibrated <- setdiff(vars_DALEC.A3.C3.H2.M1, names(fixed_pars))
final <- res$data[[length(res$data)]]$training_data$training_data
pts <- final[sample(nrow(final), n_bench, replace = TRUE), calibrated]
rownames(pts) <- NULL

time_it <- function(expr) { t0 <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t0 }

# 1. DALEC-Crop through the R wrapper used by history matching (Fortran + R overhead)
t_dalec <- time_it(generate_input_output_DALEC_Crop(pts, obs_dates = obs_dates,
                     obs_dates_names = obs_dates_names, fixed_pars = fixed_pars, year = year))

# 2. one wave's emulators (11 outputs): expectation and variance
ems_newest <- res$emulators[[1]]
t_em_wave <- time_it(for (em in ems_newest) {
  x <- pts[, names(em$ranges), drop = FALSE]; em$get_exp(x); em$get_cov(x) })

# 3. what the point proposal actually does at Wave 4: implausibility under the emulators of
#    all four waves (44 emulators), as in point_proposal.R
all_ems <- unlist(res$emulators)
all_ems <- all_ems[!names(all_ems) %in% "DS_max"]
t_impl_all <- time_it(for (em in all_ems)
  em$implausibility(pts[, names(em$ranges), drop = FALSE], targets[[em$output_name]], 3))

per_ms <- function(t) 1000 * t / n_bench
cat(sprintf("\nPer parameter set (single core, %d sets):\n", n_bench))
cat(sprintf("  DALEC-Crop run (R wrapper):                 %8.2f ms\n", per_ms(t_dalec)))
cat(sprintf("  one wave of emulators (11 outputs):         %8.2f ms\n", per_ms(t_em_wave)))
cat(sprintf("  implausibility, all 4 waves (44 emulators): %8.2f ms\n", per_ms(t_impl_all)))
cat(sprintf("  CARDAMOM, implied per MCMC iteration:       %8.2f ms  (17 h / 1e8 iterations)\n",
            17 * 3600 * 1000 / 1e8))

# how many evaluations each method needs per calibration
wt <- res$wave_table
emulator_point_evals <- sum(wt$pass_input_edcs)        # points whose implausibility was computed
cat(sprintf("\nPer calibration:\n  history matching: %s DALEC-Crop runs, %s points evaluated by the emulators\n",
            format(res$dalec_evaluations_total, big.mark = ","),
            format(emulator_point_evals, big.mark = ",")))
cat("  CARDAMOM:         1e8 DALEC-Crop runs (recommended run length)\n")
