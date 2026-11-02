# figure4_10_parameters.R
# Figure 4 (six selected parameters) and Figures 10a/10b (all calibrated parameters):
# smoothed NROY samples (history matching), CARDAMOM posterior, CARDAMOM prior, true value.
# Redraws from saved data only, so it does not run DALEC-Crop:
#   1. <scratch>/synthetic_revision/figure4_10_data.rds, if it exists (written by
#      make_figures_synthetic.R or by an earlier run of this script);
#   2. otherwise it is built from the three history matching result files, the CARDAMOM
#      parameter file of the synthetic run, synth_obs_feb.RDS and parameter_ranges.R, and
#      saved as figure4_10_data.rds for next time.
#
# Usage (from the CARDAMOM working directory):
#   source("synth_experiment/revision/figure4_10_parameters.R")

library(dplyr); library(tidyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
source(paste0(code_dir, "figure_scheme.R"))

run_dir <- file.path(scratch, "synthetic_revision")
fig410_file <- file.path(run_dir, "figure4_10_data.rds")

hm_runs <- c("HM 1a: EDCs and md" = "exp1a_md10_seed1",
             "HM 1b: EDCs"        = "exp1b_md00_seed1",
             "HM 2: md"           = "exp2_noedcs_md10")
method_levels <- c("CARDAMOM", names(hm_runs))
method_labels <- c("CARDAMOM" = "CARDAMOM", "HM 1a: EDCs and md" = "HM 1a: EDCs and md",
                   "HM 1b: EDCs" = "HM 1b: EDCs, no md", "HM 2: md" = "HM 2: no EDCs, md")

# ---- data --------------------------------------------------------------------------------------
if (file.exists(fig410_file)) {
  fig410 <- readRDS(fig410_file)
} else {
  message("Building figure4_10_data.rds from the saved runs (no DALEC-Crop runs)")
  source(paste0(code_dir, "parameter_ranges.R"))       # parameter names, ranges and priors
  fixed_names <- c("sow_day", "harvest_length")
  calibrated <- setdiff(vars_DALEC.A3.C3.H2.M1, fixed_names)

  hm_pars <- lapply(hm_runs, function(r) {
    res <- readRDS(file.path(run_dir, paste0(r, ".rds")))
    res$data[[length(res$data)]]$training_data$training_data
  })
  cardamom_synth_dir <- "CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/synth_feb/RESULTS_PROCESSED/"
  cardamom_synth_file <- list.files(cardamom_synth_dir, pattern = "parameters\\.RData$", full.names = TRUE)[1]
  e <- new.env(); load(cardamom_synth_file, envir = e)
  n_par <- length(vars_DALEC.A3.C3.H2.M1)
  cd_pars <- as.data.frame(t(do.call(cbind, lapply(seq_len(dim(e$parameters)[3]),
                                                    function(k) e$parameters[1:n_par, , k]))))
  colnames(cd_pars) <- vars_DALEC.A3.C3.H2.M1
  par_samples <- c(list(CARDAMOM = cd_pars), hm_pars)

  true_pars <- readRDS(file.path(scratch, "synth_obs_feb.RDS"))$true_pars
  true_df <- as.data.frame(as.list(true_pars))
  names(true_df) <- if (is.null(names(true_pars))) vars_DALEC.A3.C3.H2.M1 else names(true_pars)

  fig410 <- list(
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
              density = dnorm(seq(par_min, par_max, length.out = 200), prior_mean, prior_unc)))
  saveRDS(fig410, fig410_file)
}
calibrated   <- fig410$calibrated
par_long     <- fig410$par_long %>% mutate(method = factor(method, levels = method_levels))
truth_long   <- fig410$truth_long
prior_curves <- fig410$prior_curves

# ---- figures -----------------------------------------------------------------------------------
truth_colour <- "grey50"
# legend keys: a vertical line for the true parameter value, a horizontal line for all
# other entries (every layer draws its key for every entry, so all layers use this function)
draw_key_truth <- function(data, params, size) {
  if (identical(as.character(data$colour), truth_colour)) draw_key_vline(data, params, size)
  else draw_key_path(data, params, size)
}

fig10_page <- function(pars, ncol = 4, order = calibrated) {
  series <- c(method_levels, "Prior", "True parameter value")
  cols   <- c(method_colours[method_levels], "Prior" = "black", "True parameter value" = truth_colour)
  ltys   <- c(method_linetypes[method_levels], "Prior" = "dotted", "True parameter value" = "solid")
  labs_all <- c(method_labels, "Prior" = "Prior", "True parameter value" = "True parameter value")
  ggplot(filter(par_long, par_name %in% pars), aes(x = value)) +
    geom_vline(data = filter(truth_long, par_name %in% pars),
               aes(xintercept = value, colour = "True parameter value", linetype = "True parameter value"),
               linewidth = 0.7, key_glyph = draw_key_truth) +
    geom_line(data = filter(prior_curves, par_name %in% pars),
              aes(y = density, colour = "Prior", linetype = "Prior"), linewidth = 0.7,
              key_glyph = draw_key_truth) +
    geom_density(aes(colour = method, linetype = method), linewidth = 0.8, adjust = 1.2,
                 key_glyph = draw_key_truth) +
    facet_wrap(~ factor(par_name, levels = order), scales = "free", ncol = ncol, strip.position = "bottom") +
    scale_colour_manual(name = NULL, values = cols, breaks = series, labels = labs_all) +
    scale_linetype_manual(name = NULL, values = ltys, breaks = series, labels = labs_all) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
    labs(x = NULL, y = "Density") +
    theme_grey(base_size = 16) +
    theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"),
          legend.key = element_rect(fill = NA, colour = NA),   # plain line keys, no boxes
          legend.text = element_text(size = 15), axis.text = element_text(size = 11),
          strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(size = 15))
}

# Figure 4: six selected parameters
fig4_pars <- c("LMA", "frac_GPP_resp_auto", "RDRSHMAX", "NUE", "N_interc", "decomp_rate")
fig4 <- fig10_page(fig4_pars, ncol = 2, order = fig4_pars)
ggsave(file.path(run_dir, "Figure4_revision.pdf"), fig4, width = 10, height = 10)

# Figures 10a and 10b: all calibrated parameters, half on each page
half <- ceiling(length(calibrated) / 2)
ggsave(file.path(run_dir, "Figure10a_revision.pdf"), fig10_page(calibrated[1:half]), width = 11, height = 14)
ggsave(file.path(run_dir, "Figure10b_revision.pdf"),
       fig10_page(calibrated[(half + 1):length(calibrated)]), width = 11, height = 14)
message("Saved Figure4_revision.pdf, Figure10a_revision.pdf and Figure10b_revision.pdf to ", run_dir)
