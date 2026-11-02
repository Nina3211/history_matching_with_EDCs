# make_figures_synthetic.R
# Synthetic-experiment figures from the reruns (run_synthetic_all.R), in the shared
# colour-blind-safe scheme (figure_scheme.R):
#   Figure 2:  A) LAI bands and B) yield at harvest, CARDAMOM and HM Exp. 1a / 1b / 2
#   Figure 4:  six selected parameters: HM NROY samples, CARDAMOM posterior, prior, truth
#   Figure 10: all calibrated parameters, split over two pages (Reviewer 2, comment 11)
#   Figure 3:  active parameters per emulator in the final wave (Exp. 1a)
#   Figure 9:  active parameters per emulator for every wave (Exp. 1a)
#   Table:     LAI and yield summaries for all methods incl. CARDAMOM (table_synthetic_comparison.csv)
# Saves PDFs to <scratch>/synthetic_revision/.
#
# Run from the CARDAMOM working directory on the cluster: it runs DALEC-Crop for the
# CARDAMOM posterior samples (synth_feb run) and for the true parameter vector.

library(dplyr); library(tidyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"

# CARDAMOM calibration of the synthetic experiment (as in the original plotting code)
cardamom_synth_dir  <- "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/synth_feb/RESULTS_PROCESSED/"
cardamom_synth_file <- list.files(cardamom_synth_dir, pattern = "parameters\\.RData$", full.names = TRUE)[1]
if (is.na(cardamom_synth_file)) stop("No CARDAMOM parameters file found in ", cardamom_synth_dir)

run_dir <- file.path(scratch, "synthetic_revision")
hm_runs <- c("HM 1a: EDCs and md" = "exp1a_md10_seed1",
             "HM 1b: EDCs"        = "exp1b_md00_seed1",
             "HM 2: md"           = "exp2_noedcs_md10")
method_levels <- c("CARDAMOM", names(hm_runs))

source(paste0(code_dir, "figure_scheme.R"))
source(paste0(code_dir, "parameter_ranges.R"))
source(paste0(code_dir, "generate_input_output_DALEC_Crop.R"))
source("./R_functions/load_all_cardamom_functions.r")
load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year1/infofile.RData")
load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_year2_jan25/RESULTS_PROCESSED/ES_plot1_y2_2016_2018_parameters.RData")
new_PROJECT <- PROJECT

obs_dates_names <- c("date20170420", "date20170516", "date20170530",
                     "date20170612", "date20170626", "date20170710")
obs_dates  <- as.Date(sub("date", "", obs_dates_names), format = "%Y%m%d")
fixed_pars <- list(sow_day = 274, harvest_length = 342)
fixed_names <- names(fixed_pars)
calibrated <- setdiff(vars_DALEC.A3.C3.H2.M1, fixed_names)

synth <- readRDS(file.path(scratch, "synth_obs_feb.RDS"))
true_pars <- synth$true_pars

run_forward <- function(pars) {
  generate_input_output_DALEC_Crop(as.data.frame(pars)[, calibrated, drop = FALSE],
                                   obs_dates = obs_dates, obs_dates_names = obs_dates_names,
                                   fixed_pars = fixed_pars, year = 2017, is_final_wave = TRUE)
}

# ---- collect parameters, LAI curves and yields for every method ----------------------------------
message("Reading history matching runs")
hm <- lapply(hm_runs, function(r) readRDS(file.path(run_dir, paste0(r, ".rds"))))
hm_final <- lapply(hm, function(res) res$data[[length(res$data)]]$training_data)

message("Running CARDAMOM posterior samples and the truth through DALEC-Crop")
e <- new.env(); load(cardamom_synth_file, envir = e)
n_par <- length(vars_DALEC.A3.C3.H2.M1)
cd_pars <- as.data.frame(t(do.call(cbind, lapply(seq_len(dim(e$parameters)[3]),
                                                  function(k) e$parameters[1:n_par, , k]))))
colnames(cd_pars) <- vars_DALEC.A3.C3.H2.M1
# forward runs are cached, so redrawing the figures (e.g. for fonts) does not rerun DALEC-Crop;
# delete forward_runs_cache.rds to recompute
forward_cache <- file.path(run_dir, "forward_runs_cache.rds")
if (file.exists(forward_cache)) {
  cached_runs <- readRDS(forward_cache)
  cd_out <- cached_runs$cd_out
} else {
  cd_out <- run_forward(cd_pars)
}

# R-hat across the three CARDAMOM chains (Appendix A.3), as in summarise_CARDAMOM.R
rhat <- function(chains) {
  n <- nrow(chains[[1]])
  out <- sapply(colnames(chains[[1]]), function(p) {
    x <- sapply(chains, function(ch) ch[, p])
    B <- n * var(colMeans(x)); W <- mean(apply(x, 2, var))
    sqrt(((n - 1) / n * W + B / n) / W)
  })
  out[is.finite(out)]
}
cd_chains <- lapply(seq_len(dim(e$parameters)[3]), function(k) {
  m <- t(e$parameters[1:n_par, , k]); colnames(m) <- vars_DALEC.A3.C3.H2.M1; m[, calibrated] })
synthetic_rhat <- rhat(cd_chains)
message("Synthetic CARDAMOM R-hat: ", sum(synthetic_rhat < 1.2), " of ", length(synthetic_rhat),
        " parameters below 1.2 (median ", round(median(synthetic_rhat), 3),
        ", max ", round(max(synthetic_rhat), 3), ")")
if (any(synthetic_rhat >= 1.2)) print(sort(synthetic_rhat[synthetic_rhat >= 1.2], decreasing = TRUE))
readr::write_csv(tibble::tibble(par_name = names(synthetic_rhat), rhat = synthetic_rhat),
                 file.path(run_dir, "table_synthetic_rhat.csv"))

true_df  <- as.data.frame(as.list(true_pars))
names(true_df) <- if (is.null(names(true_pars))) vars_DALEC.A3.C3.H2.M1 else names(true_pars)
if (exists("cached_runs")) {
  true_out <- cached_runs$true_out
} else {
  true_out <- run_forward(true_df[c(1, 1), , drop = FALSE])   # twice: DALEC-Crop wrapper fails for a single row
  saveRDS(list(cd_out = cd_out, true_out = true_out), forward_cache)
}
true_yield <- true_out$training_data$harvest_vector[1]

lai_curves <- c(list(CARDAMOM = cd_out$LAI_curves), lapply(hm_final, `[[`, "LAI_curves"))
yields     <- c(list(CARDAMOM = cd_out$training_data$harvest_vector),
                lapply(hm_final, function(x) x$training_data$harvest_vector))
par_samples <- c(list(CARDAMOM = cd_pars), lapply(hm_final, `[[`, "training_data"))

# ---- comparison table: all methods, synthetic experiment ---------------------------------------------
# (the synthetic counterpart of the field-data table; uses the forward runs above, nothing new)
obs_vals <- sapply(synth$LAI_targets[obs_dates_names], `[[`, "val")
lai_at_obs <- c(list(CARDAMOM = cd_out$training_data[, obs_dates_names]),
                lapply(hm_final, function(x) x$training_data[, obs_dates_names]))
synthetic_comparison <- bind_rows(lapply(names(yields), function(m) {
  y <- yields[[m]]; curves <- lai_curves[[m]]
  lai_sim <- as.matrix(lai_at_obs[[m]])
  tibble(method = m,
         n_points = length(y),
         lai_band_width = mean(4 * apply(as.matrix(curves[, setdiff(names(curves), "date")]), 1, sd)),
         lai_rmse = sqrt(mean((colMeans(lai_sim) - obs_vals)^2)),
         yield_mean = mean(y),
         yield_bias = mean(y) - true_yield,
         yield_band_width = 4 * sd(y),
         true_yield_in_band = abs(true_yield - mean(y)) <= 2 * sd(y),
         true_yield_in_range = true_yield >= min(y) & true_yield <= max(y))
})) %>% mutate(lai_width_ratio_to_CARDAMOM = lai_band_width / lai_band_width[method == "CARDAMOM"],
               yield_width_ratio_to_CARDAMOM = yield_band_width / yield_band_width[method == "CARDAMOM"])
message("True yield: ", round(true_yield, 1), " gC m-2")
print(synthetic_comparison, width = Inf)
readr::write_csv(synthetic_comparison, file.path(run_dir, "table_synthetic_comparison.csv"))

# ---- shared settings: original figure style, colour-blind-safe colours ----------------------------
# Colours follow figure_scheme.R (CARDAMOM amber, HM 1a blue, HM 1b purple, HM 2 magenta);
# layout, line styles and themes follow the original manuscript figures.
method_labels <- c("CARDAMOM" = "CARDAMOM", "HM 1a: EDCs and md" = "HM 1a: EDCs and md",
                   "HM 1b: EDCs" = "HM 1b: EDCs, no md", "HM 2: md" = "HM 2: no EDCs, md")
original_theme <- theme_bw(base_size = 16) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key = element_rect(fill = NA, colour = NA),   # plain line keys, no boxes
        legend.key.width = unit(1.6, "cm"), legend.text = element_text(size = 15),
        plot.tag = element_text(face = "bold", size = 18))

# ---- Figure 2: A) LAI bands, B) yield at harvest -----------------------------------------------------
# The data for Figure 2 are saved to figure2_data.rds, and figure2_synthetic.R draws the figure
# from that file; to change only the layout later, run figure2_synthetic.R on its own.
saveRDS(list(
  bands = bind_rows(lapply(names(lai_curves), function(m) {
    lai_curves[[m]] %>% pivot_longer(-date, values_to = "LAI") %>% group_by(date) %>%
      summarise(mean = mean(LAI), lo = mean(LAI) - 2 * sd(LAI), hi = mean(LAI) + 2 * sd(LAI),
                .groups = "drop") %>% mutate(method = m)
  })),
  obs = tibble(date = obs_dates,
               LAI = sapply(synth$LAI_targets[obs_dates_names], `[[`, "val"),
               sd  = sapply(synth$LAI_targets[obs_dates_names], `[[`, "sigma")),
  yield_df = tibble(method = names(yields), mean = sapply(yields, mean), sd = sapply(yields, sd)),
  true_yield = true_yield),
  file.path(run_dir, "figure2_data.rds"))
source(paste0(code_dir, "figure2_synthetic.R"))

# ---- Figures 4 and 10: parameters ---------------------------------------------------------------------
# The data are saved to figure4_10_data.rds and figure4_10_parameters.R draws the figures from
# that file; to change only the layout later, run figure4_10_parameters.R on its own.
saveRDS(list(
  calibrated = calibrated,
  par_long = bind_rows(lapply(names(par_samples), function(m)
    par_samples[[m]][, calibrated] %>% pivot_longer(everything(), names_to = "par_name") %>%
      mutate(method = m))),
  truth_long = tibble(par_name = calibrated, value = unlist(true_df[1, calibrated])),
  prior_curves = ranges_df_with_priors %>%
    filter(par_name %in% calibrated, is.finite(prior_mean), is.finite(prior_unc),
           prior_unc > 0, prior_mean != -9999) %>%
    rowwise() %>%
    reframe(par_name = par_name,
            value = seq(par_min, par_max, length.out = 200),
            density = dnorm(seq(par_min, par_max, length.out = 200), prior_mean, prior_unc))),
  file.path(run_dir, "figure4_10_data.rds"))
source(paste0(code_dir, "figure4_10_parameters.R"))

# ---- Figures 3 and 9: active parameters (Exp. 1a), original style -------------------------------------
ems <- hm[["HM 1a: EDCs and md"]]$emulators          # newest wave first
n_waves <- length(ems)
output_order <- c(obs_dates_names, "EDC_SOM1", "EDC_SOM2", "ra_GPP_ratio", "yield_GPP_ratio", "DS_max")
active <- bind_rows(lapply(seq_len(n_waves), function(j) {
  wave <- n_waves - j + 1
  bind_rows(lapply(ems[[j]], function(em)
    tibble(wave = wave, output = em$output_name, par_name = names(em$ranges),
           active = em$active_vars)))
})) %>% mutate(par_name = factor(par_name, levels = calibrated),
               output = factor(output, levels = output_order))

# Figure 3: final wave, grey (active) and black (inactive), as in the original
fig3 <- ggplot(filter(active, wave == n_waves), aes(x = par_name, y = output,
                                                    fill = ifelse(active, "Active", "Inactive"))) +
  geom_tile(colour = "grey30") +
  scale_fill_manual(name = NULL, values = c(Active = "grey80", Inactive = "grey15")) +
  labs(x = "Parameter", y = "Output") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13),
        axis.text.y = element_text(size = 14), legend.text = element_text(size = 16),
        legend.position = "bottom", panel.grid = element_blank())
# a narrower canvas, so the text is not shrunk when the figure is scaled to the text width
ggsave(file.path(run_dir, "Figure3_revision.pdf"), fig3, width = 10, height = 6.5)

# Figure 9: all waves, green (active) and black (inactive), as in the original;
# waves stacked in one column on a full page, so the labels stay large
fig9 <- ggplot(active, aes(x = par_name, y = output, fill = active)) +
  geom_tile(colour = "grey30") +
  facet_wrap(~ wave, ncol = 1, labeller = labeller(wave = function(w) paste("Wave", w))) +
  scale_fill_manual(name = NULL, values = c(`FALSE` = "black", `TRUE` = "#33CC33"),
                    labels = c(`FALSE` = "Inactive", `TRUE` = "Active")) +
  labs(x = "Parameter", y = "Output") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y = element_text(size = 11), legend.text = element_text(size = 15),
        legend.position = "bottom", panel.grid = element_blank(),
        strip.text = element_text(hjust = 0, size = 15, face = "bold"))
ggsave(file.path(run_dir, "Figure9_revision.pdf"), fig9, width = 9, height = 15)

message("Saved Figures 2, 3, 4, 9 and 10 to ", run_dir)
