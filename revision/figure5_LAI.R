# figure5_LAI.R
# Figure 5: LAI curves for one plot per season, history matching (rerun) vs CARDAMOM.
#   HM:        the 300 final NROY points' LAI curves, saved in the rerun files
#   CARDAMOM:  the 3 x 100 posterior samples, run forward through DALEC-Crop
#              (as in the original Figure 5 code)
# Colours, line types and theme come from figure_scheme.R (colour-blind safe).
# Bands are the pointwise mean +/- 2 SD (CARDAMOM: approximate 95% posterior predictive
# band; HM: NROY envelope). Observations are shown with +/- 2 SD of the 12% observation error.
# Run from the CARDAMOM working directory on the cluster (needs DALEC-Crop).

library(dplyr); library(tidyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"

source(paste0(code_dir, "figure_scheme.R"))
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

# same settings as summarise_CARDAMOM.R (CHECK the 2017-2018 CARDAMOM folder there too)
season <- list(
  "2017" = list(dates = c("date20170420", "date20170516", "date20170530",
                          "date20170612", "date20170626", "date20170710"),
                fixed_pars = list(sow_day = 274, harvest_length = 342),
                infofile = "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year1/infofile.RData",
                cardamom_dir = "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_plots_year1_jan17/RESULTS_PROCESSED/",
                tag = "y1", xlim = as.Date(c("2016-09-29", "2017-09-08"))),
  "2018" = list(dates = c("date20180507", "date20180613", "date20180619", "date20180704"),
                fixed_pars = list(sow_day = 273, harvest_length = 329),
                infofile = "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year2/infofile.RData",
                cardamom_dir = "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_year2_jan25/RESULTS_PROCESSED/",
                tag = "y2", xlim = as.Date(c("2017-09-29", "2018-08-26"))))

# pointwise mean and band from a LAI_curves data frame (date + one column per run)
lai_band <- function(curves, method) {
  curves %>% pivot_longer(-date, values_to = "LAI") %>%
    group_by(date) %>%
    summarise(mean = mean(LAI), lo = mean(LAI) - 2 * sd(LAI), hi = mean(LAI) + 2 * sd(LAI),
              .groups = "drop") %>%
    mutate(method = method)
}

# The bands of each plot are cached in <scratch>/figure5_cache/, so DALEC-Crop only runs the
# first time a plot is drawn; later runs (e.g. to change fonts) only redraw. Delete the
# folder to recompute.
cache_dir <- file.path(scratch, "figure5_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

lai_panel <- function(plot_id, year, title = paste("Plot", plot_id), base_size = 18) {
  s <- season[[as.character(year)]]
  cache_file <- file.path(cache_dir, paste0("bands_", year, "_plot", plot_id, ".rds"))
  if (file.exists(cache_file)) {
    d <- readRDS(cache_file)
  } else {
    d <- lai_panel_data(plot_id, year)
    saveRDS(d, cache_file)
  }
  lai_panel_plot(d$bands, d$obs, s, title, base_size)
}

lai_panel_data <- function(plot_id, year) {
  s <- season[[as.character(year)]]
  setup_season_project(s$infofile)
  obs_dates <- as.Date(sub("date", "", s$dates), format = "%Y%m%d")

  # history matching: saved curves of the final NROY points
  hm_file <- file.path(scratch, paste0("ATEC_revision_", year),
                       paste0("HM_year", year - 2016, "_revision_plot", plot_id, ".rds"))
  hm <- readRDS(hm_file)
  hm_curves <- hm$data[[length(hm$data)]]$training_data$LAI_curves

  # CARDAMOM: forward runs of the posterior samples
  cd_file <- list.files(s$cardamom_dir, full.names = TRUE,
                        pattern = paste0("^ES_plot", plot_id, "_", s$tag, ".*parameters\\.RData$"))
  e <- new.env(); load(cd_file, envir = e)
  n_par <- length(vars_DALEC.A3.C3.H2.M1)
  pars <- as.data.frame(t(do.call(cbind, lapply(seq_len(dim(e$parameters)[3]),
                                                 function(k) e$parameters[1:n_par, , k]))))
  colnames(pars) <- vars_DALEC.A3.C3.H2.M1
  cd_curves <- generate_input_output_DALEC_Crop(
    pars %>% select(-c("sow_day", "harvest_length")),
    obs_dates = obs_dates, obs_dates_names = s$dates,
    fixed_pars = s$fixed_pars, year = year, is_final_wave = TRUE)$LAI_curves

  bands <- bind_rows(lai_band(hm_curves, "HM"), lai_band(cd_curves, "CARDAMOM")) %>%
    mutate(method = factor(method, levels = c("CARDAMOM", "HM")))
  obs <- atec %>% filter(PLOT_ID == plot_id, YEAR == year) %>% arrange(DOY) %>%
    transmute(date = obs_dates, LAI, sd = 0.12 * LAI)
  list(bands = bands, obs = obs)
}

lai_panel_plot <- function(bands, obs, s, title, base_size = 18) {
  ggplot(bands, aes(x = date)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = method, colour = method, linetype = method),
                alpha = 0.2, linewidth = 0.5) +
    geom_line(aes(y = mean, colour = method, linetype = method), linewidth = 1.2) +
    geom_errorbar(data = obs, aes(ymin = LAI - 2 * sd, ymax = LAI + 2 * sd), width = 6, colour = "black") +
    geom_point(data = obs, aes(y = LAI), size = 2.5, colour = "black") +
    scale_colour_methods("Calibration method") + scale_fill_methods("Calibration method") +
    scale_linetype_methods("Calibration method") +
    scale_x_date(date_labels = "%b\n%Y") +   # two-line date labels, so they do not overlap
    coord_cartesian(xlim = s$xlim) +
    labs(x = "Date", y = expression(LAI~(m^2~m^-2)), title = title) +
    theme_paper(base_size = base_size)
}
