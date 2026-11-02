# run_synthetic_all.R
# Reruns every synthetic history matching experiment with the bug fixes.
#   Exp. 1a  EDCs, 10% discrepancy                  (paper)
#   Exp. 1b  EDCs, 0% discrepancy                   (paper)
#   Exp. 2   no EDCs, 10% discrepancy               (paper; hmer's generate_new_design(),
#            as in the original Exp. 2 script; neither bug affected it, so the rerun
#            should agree with the original apart from random variation)
#   EDCs with 5%, 20%, 30% discrepancy              (sensitivity, R1.6 / R2.6)
#   Exp. 1a repeated with seeds 2-6                 (timing variability, R2.8)
# 11 calibrations, about 20-30 minutes each. Finished runs are skipped, so the
# script can be restarted after an interruption.
# Afterwards, run posthoc_synthetic_analysis.R on each saved file.

library(hmer); library(dplyr); library(purrr); library(tidyr)

out_dir <- "~/scratch/synthetic_revision/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- set-up (as in the synthetic script) ----------------------------------------------
# All code is loaded from this folder (relative to the CARDAMOM working directory).
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
source(paste0(code_dir, "parameter_ranges.R"))
source(paste0(code_dir, "edc1.R"))
source(paste0(code_dir, "edc2.R"))                       # fixed version
source(paste0(code_dir, "build_initial_LHD.R"))
source(paste0(code_dir, "point_proposal.R"))             # revision version
source(paste0(code_dir, "generate_input_output_DALEC_Crop.R"))
source(paste0(code_dir, "hm_runner.R"))
source("./R_functions/load_all_cardamom_functions.r")

year <- 2017
obs_dates_names <- c("date20170420", "date20170516", "date20170530",
                     "date20170612", "date20170626", "date20170710")
obs_dates <- as.Date(c("2017-04-20", "2017-05-16", "2017-05-30",
                       "2017-06-12", "2017-06-26", "2017-07-10"))
fixed_pars <- list(sow_day = 274, harvest_length = 342)
load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year1/infofile.RData")
load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_year2_jan25/RESULTS_PROCESSED/ES_plot1_y2_2016_2018_parameters.RData")
new_PROJECT <- PROJECT

synth     <- readRDS("~/scratch/synth_obs_feb.RDS")
true_pars <- synth$true_pars
targets   <- synth$LAI_targets
targets$ra_GPP_ratio    <- list(val = 0.54, sigma = 0.12 * (1 / sqrt(3)))
targets$yield_GPP_ratio <- list(val = 0.38, sigma = 0.087 * (1 / sqrt(3)))
targets$EDC_SOM1        <- list(val = 0, sigma = log(10) / 3)
targets$EDC_SOM2        <- list(val = 0, sigma = 0.2 / 3)
targets$DS_max          <- list(val = 1.9, sigma = 1)

# ---- the runs -----------------------------------------------------------------------------
runs <- tribble(
  ~name,                 ~discs_prop, ~use_edcs, ~seed,
  "exp1a_md10_seed1",    0.10,        TRUE,      1,
  "exp1b_md00_seed1",    0.00,        TRUE,      1,
  "exp2_noedcs_md10",    0.10,        FALSE,     1,
  "edcs_md05_seed1",     0.05,        TRUE,      1,
  "edcs_md20_seed1",     0.20,        TRUE,      1,
  "edcs_md30_seed1",     0.30,        TRUE,      1,
  "exp1a_md10_seed2",    0.10,        TRUE,      2,
  "exp1a_md10_seed3",    0.10,        TRUE,      3,
  "exp1a_md10_seed4",    0.10,        TRUE,      4,
  "exp1a_md10_seed5",    0.10,        TRUE,      5,
  "exp1a_md10_seed6",    0.10,        TRUE,      6)

for (i in seq_len(nrow(runs))) {
  r <- runs[i, ]
  file <- paste0(out_dir, r$name, ".rds")
  if (file.exists(file)) { message("Skipping ", r$name, " (already done)"); next }
  message("\n==== ", r$name, " (", i, " of ", nrow(runs), ") ====")

  res <- tryCatch(
    run_history_matching(targets = targets, obs_dates = obs_dates,
                         obs_dates_names = obs_dates_names, year = year,
                         fixed_pars = fixed_pars, discs_prop = r$discs_prop,
                         use_edcs = r$use_edcs, range_buffer = 0.025,   # as in the synthetic script
                         seed = r$seed),
    error = function(e) { message("FAILED: ", conditionMessage(e)); NULL })
  if (is.null(res)) next

  res$true_pars <- true_pars
  res$true_yield <- synth$true_yield      # NULL if not stored in synth_obs_feb.RDS
  safe_save(res, file)    # emulators are kept: the truth check and theta sensitivity need them
}
