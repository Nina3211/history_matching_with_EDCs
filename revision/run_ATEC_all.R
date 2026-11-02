# run_ATEC_all.R
# Reruns history matching for every ATEC plot with the bug fixes.
# Replaces ATEC_year1_hm_ranges.R and ATEC_year2_hm_ranges.R; settings as before
# (full parameter ranges, 10% discrepancy, 1% range buffer).
#
# 100 calibrations at 20-30 minutes each: run the two seasons in two separate R
# sessions, e.g.   year <- 2017; source("run_ATEC_all.R")
#                  year <- 2018; source("run_ATEC_all.R")
# Finished plots are skipped, so a session can be restarted after an interruption.

library(hmer); library(dplyr); library(purrr); library(tidyr)

if (!exists("year"))     year     <- 2017   # 2017 = 2016-2017 season, 2018 = 2017-2018 season
if (!exists("plot_ids")) plot_ids <- 1:50

out_dir <- paste0("~/scratch/ATEC_revision_", year, "/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# gC m-2 per t ha-1 for the observed yields (as in the old scripts: 0.48 * 100).
# Make Appendix E.1 and Figure 6 use the same factor.
yield_factor <- 0.48 * 100

# ---- set-up (as in the ATEC scripts) -----------------------------------------------------
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

if (year == 2017) {
  obs_dates_names <- c("date20170420", "date20170516", "date20170530",
                       "date20170612", "date20170626", "date20170710")
  fixed_pars <- list(sow_day = 274, harvest_length = 342)
  load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year1/infofile.RData")
} else {
  obs_dates_names <- c("date20180507", "date20180613", "date20180619", "date20180704")
  fixed_pars <- list(sow_day = 273, harvest_length = 329)
  load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year2/infofile.RData")
}
obs_dates <- as.Date(sub("date", "", obs_dates_names), format = "%Y%m%d")

load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_year2_jan25/RESULTS_PROCESSED/ES_plot1_y2_2016_2018_parameters.RData")
new_PROJECT <- PROJECT

atec <- readr::read_csv("synth_experiment/synth_jan/trial_plot_data_2017_2018.csv",
                        show_col_types = FALSE) %>%
  mutate(LAI_unc = 0.12 * LAI, yield_conv = yield_factor * YIELD)

make_targets <- function(plot_id) {
  d <- atec %>% filter(PLOT_ID == plot_id, YEAR == year) %>% arrange(DOY)
  targets <- map2(d$LAI, d$LAI_unc, ~ list(val = .x, sigma = .y))
  names(targets) <- obs_dates_names
  c(targets, list(
    ra_GPP_ratio    = list(val = 0.54, sigma = 0.12 * (1 / sqrt(3))),
    yield_GPP_ratio = list(val = 0.38, sigma = 0.087 * (1 / sqrt(3))),
    EDC_SOM1        = list(val = 0, sigma = log(10) / 3),
    EDC_SOM2        = list(val = 0, sigma = 0.2 / 3),
    DS_max          = list(val = 1.9, sigma = 1)))
}

# ---- run every plot -------------------------------------------------------------------------
prefix <- paste0("HM_year", year - 2016, "_revision_plot")   # HM_year1_revision_plot5.rds etc.
for (plot_id in plot_ids) {
  file <- paste0(out_dir, prefix, plot_id, ".rds")
  if (file.exists(file)) { message("Skipping plot ", plot_id, " (already done)"); next }
  message("\n==== Season ", year - 1, "-", year, ", plot ", plot_id, " ====")

  res <- tryCatch(
    run_history_matching(targets = make_targets(plot_id), obs_dates = obs_dates,
                         obs_dates_names = obs_dates_names, year = year,
                         fixed_pars = fixed_pars, discs_prop = 0.1, use_edcs = TRUE,
                         range_buffer = 0.01, seed = plot_id),
    error = function(e) { message("FAILED plot ", plot_id, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) next

  rows <- atec %>% filter(PLOT_ID == plot_id, YEAR == year)
  res$plot_info <- list(plot_id = plot_id, year = year,
                        true_yield = unique(rows$yield_conv),     # same name as in the old files
                        N_level = unique(rows$FERTILISER_LEVEL),
                        yield_factor = yield_factor)
  res$emulators <- NULL   # not needed for the ATEC analysis; keeps files small
  safe_save(res, file)
}
