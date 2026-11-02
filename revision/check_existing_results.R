# check_existing_results.R
# Summarises the ATEC history matching results (by default the reruns); nothing is run.
# Part 1: EDC check on the final NROY points (should be all zero after the fix)
# Part 2: metrics for the reviewers (R2.7), by season and nitrogen level
#
# Each saved file is read ONCE; only the final NROY points and plot_info are kept
# in memory (the emulators are dropped straight away).
# Saved files contain $data (final wave only, or all waves) and
# $plot_info = list(plot_id, true_yield, N_level); get_td() handles both layouts.

library(dplyr); library(purrr); library(tidyr)   # purrr >= 1.0 for .progress and list_rbind()

# Location of the scratch folder: "~/scratch" on the cluster, e.g. "Z:/scratch" on Windows.
if (!exists("scratch")) scratch <- "~/scratch"
# All code is loaded from this folder (relative to the CARDAMOM working directory).
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
source(paste0(code_dir, "parameter_ranges.R"))   # ranges_df_with_priors

# The rerun results from run_ATEC_all.R. To check the old runs instead, use
# folders ~/scratch/ATEC_year1 / ATEC_year2 and prefixes HM_year1/2_no_ranges_june_plot.
# observed yield in gC m-2 per t ha-1, as in Figure 6 (0.425 * 100); recomputed here from the
# trial data rather than taken from plot_info$true_yield, which the HM scripts stored with 48
if (!exists("yield_factor")) yield_factor <- 42.5

if (!exists("which_runs")) which_runs <- "revision"     # "revision" (reruns) or "old" (paper's runs)
if (which_runs == "revision") {
  dir_2017 <- paste0(scratch, "/ATEC_revision_2017");  prefix_2017 <- "HM_year1_revision_plot"
  dir_2018 <- paste0(scratch, "/ATEC_revision_2018");  prefix_2018 <- "HM_year2_revision_plot"
} else {
  dir_2017 <- paste0(scratch, "/ATEC_year1");  prefix_2017 <- "HM_year1_no_ranges_june_plot"
  dir_2018 <- paste0(scratch, "/ATEC_year2");  prefix_2018 <- "HM_year2_no_ranges_june_plot"
}

files <- bind_rows(
  tibble(year = 2017, file = list.files(dir_2017, full.names = TRUE,
                                        pattern = paste0("^", prefix_2017, "\\d+\\.rds$"))),
  tibble(year = 2018, file = list.files(dir_2018, full.names = TRUE,
                                        pattern = paste0("^", prefix_2018, "\\d+\\.rds$"))))
# seasons to include: 2017 = 2016-2017, 2018 = 2017-2018
# "years" is also a function in the chron package (loaded with CARDAMOM): only use a numeric setting
if (!exists("years") || !is.numeric(years)) years <- c(2017, 2018)
files <- files %>% filter(year %in% years)
print(count(files, year))   # should be 50 per season

obs_dates_names <- list(
  "2017" = c("date20170420", "date20170516", "date20170530", "date20170612", "date20170626", "date20170710"),
  "2018" = c("date20180507", "date20180613", "date20180619", "date20180704"))

atec <- readr::read_csv("synth_experiment/synth_jan/trial_plot_data_2017_2018.csv", show_col_types = FALSE)
orig <- ranges_df_with_priors %>% filter(par_max > par_min)

# the data frame of the final NROY points, for either way the results were saved:
#   data = list(training_data = list(training_data, LAI_curves))   (final wave only)
#   data = data_list (wave0, wave1, ..., final wave)               (all waves)
get_td <- function(d) {
  x <- d$data
  if (!"training_data" %in% names(x)) x <- x[[length(x)]]   # all waves: take the final one
  y <- x$training_data
  if (is.data.frame(y)) return(y)
  if (is.list(y) && is.data.frame(y$training_data)) return(y$training_data)
  stop("Could not find the final NROY points; check str(d, max.level = 3)")
}

# ---- read every file once and compute all per-plot quantities ------------------------------
summarise_file <- function(year, file) {
  d <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(d)) {                     # truncated file, e.g. from a full disk
    warning("Skipping unreadable file ", basename(file), ": delete it and rerun that plot")
    unreadable <<- c(unreadable, file)
    return(NULL)
  }
  td   <- get_td(d)
  info <- d$plot_info
  # LAI curves of the final NROY points (date + one column per point), if saved
  fin <- if ("training_data" %in% names(d$data)) d$data else d$data[[length(d$data)]]
  curves <- if (is.list(fin$training_data)) fin$training_data$LAI_curves else NULL
  lai_band_width <- if (is.null(curves)) NA else
    mean(4 * apply(as.matrix(curves[, setdiff(names(curves), "date")]), 1, sd))
  # logged by the reruns (NA for files from the old scripts)
  run_log <- list(
    proportion_retained = if (!is.null(d$wave_table)) tail(d$wave_table$proportion_retained, 1) else NA,
    dalec_evaluations   = if (!is.null(d$dalec_evaluations_total)) d$dalec_evaluations_total else NA,
    minutes             = if (!is.null(d$timings)) d$timings[["total"]] / 60 else NA)
  rm(d)
  if (is.null(info$plot_id)) {          # not a saved ATEC plot (e.g. a test or synthetic run)
    warning("Skipping ", basename(file), ": no plot_info$plot_id")
    return(NULL)
  }

  y     <- td$harvest_vector
  dates <- obs_dates_names[[as.character(year)]]
  lai_obs  <- atec %>% filter(PLOT_ID == info$plot_id, YEAR == year) %>% arrange(DOY) %>% pull(LAI)
  lai_sim  <- as.matrix(td[, dates])
  lai_mean <- colMeans(lai_sim)
  lai_sd   <- apply(lai_sim, 2, sd)
  p <- orig %>% filter(par_name %in% names(td))

  list(
    edc = tibble(year = year, file = basename(file), n_points = nrow(td),
                 n_fail_SOM1 = sum(abs(td$EDC_SOM1) > log(10)),
                 n_fail_SOM2 = sum(abs(td$EDC_SOM2) > 0.2),     # should be 0
                 n_fail_DS   = sum(td$DS_max < 1.9)),           # should be 0
    plot = tibble(year = year, plot_id = info$plot_id, N_level = info$N_level,
                  obs_yield = yield_factor * unique(atec$YIELD[atec$PLOT_ID == info$plot_id & atec$YEAR == year]),
                  yield_mean = mean(y), yield_lo = mean(y) - 2 * sd(y), yield_hi = mean(y) + 2 * sd(y),
                  yield_min = min(y), yield_max = max(y),
                  lai_rmse = sqrt(mean((lai_mean - lai_obs)^2)),
                  lai_in_band = mean(abs(lai_obs - lai_mean) <= 2 * lai_sd),
                  lai_band_width = lai_band_width,      # mean daily width of mean +/- 2 SD
                  proportion_retained = run_log$proportion_retained,
                  dalec_evaluations = run_log$dalec_evaluations,
                  minutes = run_log$minutes),
    ranges = tibble(year = year, plot_id = info$plot_id, par_name = p$par_name,
                    frac = map2_dbl(p$par_name, p$par_max - p$par_min, ~ diff(range(td[[.x]])) / .y)))
}

unreadable <- character(0)
message("Reading ", nrow(files), " saved plots (each file once)")
per_file <- pmap(files, summarise_file, .progress = "Reading plots") %>% compact()
if (length(unreadable) > 0) {
  message("\nUnreadable files (delete them; run_ATEC_all.R will then redo these plots):")
  message(paste(" ", unreadable, collapse = "\n"))
}

som1_check <- purrr::map(per_file, "edc")    %>% list_rbind()
plots      <- purrr::map(per_file, "plot")   %>% list_rbind()
ranges_all <- purrr::map(per_file, "ranges") %>% list_rbind()

# ---- Part 1: EDC_SOM1 ----------------------------------------------------------------------
print(som1_check %>% group_by(year) %>%
        summarise(plots = n(), plots_with_SOM1_failures = sum(n_fail_SOM1 > 0),
                  points_failing = sum(n_fail_SOM1), points_total = sum(n_points)))
# With the fixed edc2.R this should be 0 for every plot.

# ---- Part 2: metrics -----------------------------------------------------------------------
by_level <- plots %>% group_by(year, N_level) %>% summarise(
  n_plots     = n(),
  bias        = mean(yield_mean - obs_yield),
  rmse        = sqrt(mean((yield_mean - obs_yield)^2)),
  in_band     = mean(obs_yield >= yield_lo & obs_yield <= yield_hi),
  in_range    = mean(obs_yield >= yield_min & obs_yield <= yield_max),
  band_width  = mean(yield_hi - yield_lo),
  lai_rmse    = mean(lai_rmse),
  lai_in_band = mean(lai_in_band),
  lai_band_width = mean(lai_band_width),
  .groups = "drop")
print(by_level)

by_year <- plots %>% group_by(year) %>%
  summarise(rmse = sqrt(mean((yield_mean - obs_yield)^2)),
            in_band = mean(obs_yield >= yield_lo & obs_yield <= yield_hi),
            spearman = cor(yield_mean, obs_yield, method = "spearman"),
            lai_band_width = mean(lai_band_width),
            n_plots = n(),
            median_proportion_retained = median(proportion_retained),
            median_dalec_evaluations = median(dalec_evaluations),
            median_minutes = median(minutes), max_minutes = max(minutes),
            .groups = "drop")
print(by_year)

# marginal range of the final NROY sample as a fraction of the original range
range_frac <- ranges_all %>% group_by(year, par_name) %>%
  summarise(median_frac = median(frac), .groups = "drop") %>% arrange(year, median_frac)
print(range_frac, n = 80)

saveRDS(list(som1_check = som1_check, plots = plots, by_level = by_level,
             by_year = by_year, range_frac = range_frac, ranges_all = ranges_all),
        paste0(scratch, "/ATEC_metrics_", which_runs, ".rds"))

# CARDAMOM: build the same `plots` table from its posterior samples and forward runs,
# add method = "CARDAMOM" / "HM", and bind_rows() before summarising.
