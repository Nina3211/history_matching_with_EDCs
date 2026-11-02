# summarise_CARDAMOM.R
# The CARDAMOM half of the field-data comparison (R2.5, R2.7, R3.3):
#   1. R-hat (potential scale reduction factor) per plot and parameter, from the 3 x 100
#      stored samples (Appendix A.3)
#   2. forward runs of the 300 posterior samples per plot (as in the Figure 5/6 code) and
#      the same per-plot metrics as check_existing_results.R computes for history matching
#   3. one combined table by method, season and nitrogen level
#
# Run check_existing_results.R (which_runs = "revision") first; this script reads its output.
# Needs the CARDAMOM working directory and DALEC-Crop, like the history matching runs.
# Run it on the cluster: it runs DALEC-Crop 300 times per plot.

library(dplyr); library(purrr); library(tidyr); library(tibble)

if (!exists("scratch"))      scratch      <- "~/scratch"
if (!exists("code_dir"))     code_dir     <- "synth_experiment/revision/"
if (!exists("yield_factor")) yield_factor <- 42.5      # as in Figure 6 and check_existing_results.R

# CARDAMOM output folders and file tags per season.
# 2016-2017 is the folder used in the Figure 5/6 code. CHECK the 2017-2018 folder name:
cardamom_runs <- list(
  "2017" = list(dir = "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_plots_year1_jan17/RESULTS_PROCESSED/",
                tag = "y1"),
  "2018" = list(dir = "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_year2_jan25/RESULTS_PROCESSED/",
                tag = "y2"))

source(paste0(code_dir, "parameter_ranges.R"))
source(paste0(code_dir, "generate_input_output_DALEC_Crop.R"))
source("./R_functions/load_all_cardamom_functions.r")
# Season set-up exactly as in run_ATEC_all.R: the season's infofile, then the drivers file,
# then new_PROJECT. Done separately for each season, as PROJECT differs between seasons.
drivers_file <- "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_year2_jan25/RESULTS_PROCESSED/ES_plot1_y2_2016_2018_parameters.RData"
setup_season_project <- function(infofile) {
  load(infofile, envir = .GlobalEnv)
  load(drivers_file, envir = .GlobalEnv)
  assign("new_PROJECT", get("PROJECT", envir = .GlobalEnv), envir = .GlobalEnv)
}

atec <- readr::read_csv("synth_experiment/synth_jan/trial_plot_data_2017_2018.csv", show_col_types = FALSE)
orig <- ranges_df_with_priors %>% filter(par_max > par_min)
fixed_names <- c("sow_day", "harvest_length")

season_setup <- list(
  "2017" = list(dates = c("date20170420", "date20170516", "date20170530",
                          "date20170612", "date20170626", "date20170710"),
                fixed_pars = list(sow_day = 274, harvest_length = 342),
                infofile = "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year1/infofile.RData"),
  "2018" = list(dates = c("date20180507", "date20180613", "date20180619", "date20180704"),
                fixed_pars = list(sow_day = 273, harvest_length = 329),
                infofile = "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year2/infofile.RData"))

# R-hat as in Appendix A.3 (Gelman and Rubin): chains = list of (samples x parameters) matrices
rhat <- function(chains) {
  n <- nrow(chains[[1]])
  out <- sapply(colnames(chains[[1]]), function(p) {
    x <- sapply(chains, function(ch) ch[, p])
    B <- n * var(colMeans(x)); W <- mean(apply(x, 2, var))
    sqrt(((n - 1) / n * W + B / n) / W)
  })
  out[is.finite(out)]          # drops parameters that are fixed in all chains
}

summarise_cardamom_plot <- function(file, year) {
  setup <- season_setup[[as.character(year)]]
  plot_id <- as.integer(sub(".*ES_plot(\\d+)_.*", "\\1", basename(file)))
  e <- new.env(); load(file, envir = e)
  P <- e$parameters                                   # [parameter, sample, chain]
  n_par <- length(vars_DALEC.A3.C3.H2.M1)
  chains <- lapply(seq_len(dim(P)[3]), function(k) {
    m <- t(P[1:n_par, , k]); colnames(m) <- vars_DALEC.A3.C3.H2.M1; m })
  calibrated <- setdiff(vars_DALEC.A3.C3.H2.M1, fixed_names)
  rh <- rhat(lapply(chains, function(m) m[, calibrated]))

  pars <- as.data.frame(do.call(rbind, chains))
  out <- generate_input_output_DALEC_Crop(
    pars %>% select(-all_of(fixed_names)),
    obs_dates = as.Date(sub("date", "", setup$dates), format = "%Y%m%d"),
    obs_dates_names = setup$dates, fixed_pars = setup$fixed_pars,
    year = year, is_final_wave = TRUE)
  td <- out$training_data
  y  <- td$harvest_vector
  curves <- out$LAI_curves
  lai_band_width <- mean(4 * apply(as.matrix(curves[, setdiff(names(curves), "date")]), 1, sd))

  rows <- atec %>% filter(PLOT_ID == plot_id, YEAR == year)
  lai_obs  <- rows %>% arrange(DOY) %>% pull(LAI)
  lai_sim  <- as.matrix(td[, setup$dates])
  lai_mean <- colMeans(lai_sim); lai_sd <- apply(lai_sim, 2, sd)
  p <- orig %>% filter(par_name %in% calibrated)

  list(
    plot = tibble(method = "CARDAMOM", year = year, plot_id = plot_id,
                  N_level = unique(rows$FERTILISER_LEVEL),
                  obs_yield = yield_factor * unique(rows$YIELD),
                  yield_mean = mean(y), yield_lo = mean(y) - 2 * sd(y), yield_hi = mean(y) + 2 * sd(y),
                  yield_min = min(y), yield_max = max(y),
                  lai_rmse = sqrt(mean((lai_mean - lai_obs)^2)),
                  lai_in_band = mean(abs(lai_obs - lai_mean) <= 2 * lai_sd),
                  lai_band_width = lai_band_width,      # mean daily width of mean +/- 2 SD
                  rhat_max = max(rh), rhat_share_below_1.2 = mean(rh < 1.2)),
    rhat = tibble(year = year, plot_id = plot_id, par_name = names(rh), rhat = rh),
    # width of the central 95% credible interval as a fraction of the original range
    ranges = tibble(method = "CARDAMOM", year = year, plot_id = plot_id, par_name = p$par_name,
                    frac = map2_dbl(p$par_name, p$par_max - p$par_min,
                                    ~ diff(quantile(pars[[.x]], c(0.025, 0.975))) / .y)))
}

# ---- run the seasons ---------------------------------------------------------------------------
# "years" is also a function in the chron package (loaded with CARDAMOM): only use a numeric setting
if (!exists("years") || !is.numeric(years)) years <- c(2017, 2018)   # 2017 = 2016-2017, 2018 = 2017-2018
per_plot <- list()
for (yr in intersect(names(cardamom_runs), as.character(years))) {
  setup_season_project(season_setup[[yr]]$infofile)
  files <- list.files(cardamom_runs[[yr]]$dir, full.names = TRUE,
                      pattern = paste0("^ES_plot\\d+_", cardamom_runs[[yr]]$tag, ".*parameters\\.RData$"))
  message(yr, ": ", length(files), " CARDAMOM plots")   # should be 50
  per_plot[[yr]] <- purrr::map(files, summarise_cardamom_plot, year = as.integer(yr),
                        .progress = paste("CARDAMOM", yr))
}
per_plot <- unlist(per_plot, recursive = FALSE)
cardamom_plots  <- purrr::map(per_plot, "plot")   %>% list_rbind()
cardamom_rhat   <- purrr::map(per_plot, "rhat")   %>% list_rbind()
cardamom_ranges <- purrr::map(per_plot, "ranges") %>% list_rbind()

# ---- R-hat summary (Appendix A.3) -----------------------------------------------------------------
rhat_summary <- cardamom_rhat %>% group_by(year) %>%
  summarise(plot_parameter_pairs = n(), share_below_1.2 = mean(rhat < 1.2),
            median_rhat = median(rhat), max_rhat = max(rhat), .groups = "drop")
print(rhat_summary)
print(cardamom_rhat %>% filter(rhat >= 1.2) %>% count(par_name, sort = TRUE))   # which parameters

# ---- combine with history matching --------------------------------------------------------------
hm <- readRDS(paste0(scratch, "/ATEC_metrics_revision.rds"))
hm_plots <- hm$plots %>% filter(year %in% years) %>% mutate(method = "HM")
# Each method keeps every plot it calibrated (e.g. 2016-2017 plot 41: CARDAMOM only). in_both marks plots calibrated by both methods;
# paired comparisons (bootstrap differences) use only those.
common <- inner_join(distinct(hm_plots, year, plot_id), distinct(cardamom_plots, year, plot_id),
                     by = c("year", "plot_id")) %>% mutate(in_both = TRUE)
only_one <- full_join(distinct(hm_plots, year, plot_id) %>% mutate(HM = TRUE),
                      distinct(cardamom_plots, year, plot_id) %>% mutate(CARDAMOM = TRUE),
                      by = c("year", "plot_id")) %>% filter(is.na(HM) | is.na(CARDAMOM))
if (nrow(only_one) > 0) {
  message("Plots calibrated by one method only:")
  print(only_one)
}
plots_both <- bind_rows(hm_plots, cardamom_plots) %>%
  left_join(common, by = c("year", "plot_id")) %>%
  mutate(in_both = !is.na(in_both))

by_level_both <- plots_both %>% group_by(method, year, N_level) %>% summarise(
  n_plots     = n(),
  bias        = mean(yield_mean - obs_yield),
  rmse        = sqrt(mean((yield_mean - obs_yield)^2)),
  in_band     = mean(obs_yield >= yield_lo & obs_yield <= yield_hi),
  in_range    = mean(obs_yield >= yield_min & obs_yield <= yield_max),
  band_width  = mean(yield_hi - yield_lo),
  lai_rmse    = mean(lai_rmse),
  lai_in_band = mean(lai_in_band),
  lai_band_width = mean(lai_band_width),
  .groups = "drop") %>% arrange(year, N_level, method)
print(by_level_both, n = 40, width = Inf)

by_year_both <- plots_both %>% group_by(method, year) %>% summarise(
  n_plots = n(),
  rmse = sqrt(mean((yield_mean - obs_yield)^2)),
  in_band = mean(obs_yield >= yield_lo & obs_yield <= yield_hi),
  band_width = mean(yield_hi - yield_lo),
  spearman = cor(yield_mean, obs_yield, method = "spearman"),
  lai_rmse = mean(lai_rmse), lai_in_band = mean(lai_in_band),
  lai_band_width = mean(lai_band_width),
  .groups = "drop")
print(by_year_both, width = Inf)

# parameter constraint: HM marginal NROY range vs CARDAMOM 95% credible interval
ranges_both <- bind_rows(hm$ranges_all %>% filter(year %in% years) %>% mutate(method = "HM"),
                         cardamom_ranges) %>%
  group_by(method, year, par_name) %>% summarise(median_frac = median(frac), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = median_frac) %>% arrange(year, HM)
print(ranges_both, n = 80)

saveRDS(list(cardamom_plots = cardamom_plots, cardamom_rhat = cardamom_rhat,
             rhat_summary = rhat_summary, plots_both = plots_both,
             by_level_both = by_level_both, by_year_both = by_year_both,
             ranges_both = ranges_both),
        paste0(scratch, "/ATEC_comparison_HM_CARDAMOM.rds"))
readr::write_csv(by_level_both, paste0(scratch, "/table_ATEC_by_level.csv"))
readr::write_csv(by_year_both,  paste0(scratch, "/table_ATEC_by_year.csv"))
