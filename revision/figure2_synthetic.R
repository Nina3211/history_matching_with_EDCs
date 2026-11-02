# figure2_synthetic.R
# Figure 2: A) LAI bands and B) yield at harvest for CARDAMOM and HM Exp. 1a / 1b / 2.
# Redraws from saved data only, so it does not run DALEC-Crop or load CARDAMOM:
#   1. <scratch>/synthetic_revision/figure2_data.rds, if it exists (written by
#      make_figures_synthetic.R or by an earlier run of this script);
#   2. otherwise it is built from forward_runs_cache.rds (CARDAMOM and truth forward runs,
#      written by make_figures_synthetic.R), the three history matching result files and
#      synth_obs_feb.RDS, and saved as figure2_data.rds for next time.
#
# Usage (from the CARDAMOM working directory):
#   source("synth_experiment/revision/figure2_synthetic.R")
# Options, set before sourcing:
#   fig2_width_B  relative width of panel B (default 1.1; panel A is 1.15)
#   fig2_labels_B "full" (two-line method names, default) or "short" (HM 1a, HM 1b, HM 2)

library(dplyr); library(tidyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
source(paste0(code_dir, "figure_scheme.R"))

run_dir <- file.path(scratch, "synthetic_revision")
fig2_file <- file.path(run_dir, "figure2_data.rds")

hm_runs <- c("HM 1a: EDCs and md" = "exp1a_md10_seed1",
             "HM 1b: EDCs"        = "exp1b_md00_seed1",
             "HM 2: md"           = "exp2_noedcs_md10")
method_levels <- c("CARDAMOM", names(hm_runs))
method_labels <- c("CARDAMOM" = "CARDAMOM", "HM 1a: EDCs and md" = "HM 1a: EDCs and md",
                   "HM 1b: EDCs" = "HM 1b: EDCs, no md", "HM 2: md" = "HM 2: no EDCs, md")

# ---- data --------------------------------------------------------------------------------------
if (file.exists(fig2_file)) {
  fig2_data <- readRDS(fig2_file)
} else {
  forward_cache <- file.path(run_dir, "forward_runs_cache.rds")
  if (!file.exists(forward_cache))
    stop("Neither figure2_data.rds nor forward_runs_cache.rds found in ", run_dir,
         ": run make_figures_synthetic.R once first")
  message("Building figure2_data.rds from the saved runs (no DALEC-Crop runs)")
  runs <- readRDS(forward_cache)
  hm_final <- lapply(hm_runs, function(r) {
    res <- readRDS(file.path(run_dir, paste0(r, ".rds")))
    res$data[[length(res$data)]]$training_data
  })
  lai_curves <- c(list(CARDAMOM = runs$cd_out$LAI_curves), lapply(hm_final, `[[`, "LAI_curves"))
  yields     <- c(list(CARDAMOM = runs$cd_out$training_data$harvest_vector),
                  lapply(hm_final, function(x) x$training_data$harvest_vector))
  obs_dates_names <- c("date20170420", "date20170516", "date20170530",
                       "date20170612", "date20170626", "date20170710")
  synth <- readRDS(file.path(scratch, "synth_obs_feb.RDS"))
  fig2_data <- list(
    bands = bind_rows(lapply(names(lai_curves), function(m) {
      lai_curves[[m]] %>% pivot_longer(-date, values_to = "LAI") %>% group_by(date) %>%
        summarise(mean = mean(LAI), lo = mean(LAI) - 2 * sd(LAI), hi = mean(LAI) + 2 * sd(LAI),
                  .groups = "drop") %>% mutate(method = m)
    })),
    obs = tibble(date = as.Date(sub("date", "", obs_dates_names), format = "%Y%m%d"),
                 LAI = sapply(synth$LAI_targets[obs_dates_names], `[[`, "val"),
                 sd  = sapply(synth$LAI_targets[obs_dates_names], `[[`, "sigma")),
    yield_df = tibble(method = names(yields), mean = sapply(yields, mean), sd = sapply(yields, sd)),
    true_yield = runs$true_out$training_data$harvest_vector[1])
  saveRDS(fig2_data, fig2_file)
}
bands <- fig2_data$bands %>% mutate(method = factor(method, levels = method_levels))
obs <- fig2_data$obs
yield_df <- fig2_data$yield_df %>% mutate(method = factor(method, levels = method_levels))
true_yield <- fig2_data$true_yield

# ---- figure ------------------------------------------------------------------------------------
if (!exists("fig2_width_B"))  fig2_width_B  <- 1.1
if (!exists("fig2_labels_B")) fig2_labels_B <- "full"
labels_B <- if (fig2_labels_B == "short") {
  c("CARDAMOM" = "CARDAMOM", "HM 1a: EDCs and md" = "HM 1a", "HM 1b: EDCs" = "HM 1b", "HM 2: md" = "HM 2")
} else {
  c("CARDAMOM" = "CARDAMOM", "HM 1a: EDCs and md" = "HM 1a:\nEDCs, md",
    "HM 1b: EDCs" = "HM 1b:\nEDCs, no md", "HM 2: md" = "HM 2:\nno EDCs, md")
}

original_theme <- theme_bw(base_size = 16) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key = element_rect(fill = NA, colour = NA),   # plain line keys, no boxes
        legend.key.width = unit(1.6, "cm"), legend.text = element_text(size = 15),
        plot.tag = element_text(face = "bold", size = 18))

fig2a <- ggplot(bands, aes(x = date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = method), fill = "grey40", alpha = 0.07,
              show.legend = FALSE) +
  geom_line(aes(y = lo, colour = method, linetype = method), linewidth = 0.5, show.legend = FALSE) +
  geom_line(aes(y = hi, colour = method, linetype = method), linewidth = 0.5, show.legend = FALSE) +
  geom_line(aes(y = mean, colour = method, linetype = method), linewidth = 1.4, key_glyph = "path") +
  geom_errorbar(data = obs, aes(ymin = LAI - 2 * sd, ymax = LAI + 2 * sd), width = 8, colour = "black") +
  geom_point(data = obs, aes(y = LAI), size = 1.5, colour = "black") +
  scale_colour_methods(name = NULL, labels = method_labels) +
  scale_linetype_manual(name = NULL, labels = method_labels, values = method_linetypes[method_levels]) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
  coord_cartesian(xlim = as.Date(c("2016-10-01", "2017-09-08"))) +
  labs(x = "Date", y = expression(LAI~(m^2~m^-2)), tag = "A") +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  original_theme

fig2b <- ggplot(yield_df, aes(x = method, y = mean, colour = method)) +
  geom_hline(yintercept = true_yield, colour = "black", linewidth = 1.2) +
  geom_errorbar(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd), width = 0.3, linewidth = 2) +
  geom_point(size = 3.5) +
  scale_colour_methods(name = NULL, labels = method_labels) +
  # methods spread across the full width of the panel, with room between the labels
  scale_x_discrete(labels = labels_B, expand = expansion(add = 0.5)) +
  labs(x = NULL, y = expression("Yield at harvest (gC m"^{-2}*")"), tag = "B") +
  original_theme +
  theme(legend.position = "none", axis.text.x = element_text(size = 13, lineheight = 0.9))

# panel B gets about as much width as panel A, so its labels have room
fig2 <- ggpubr::ggarrange(fig2a, fig2b, widths = c(1.15, fig2_width_B), common.legend = TRUE, legend = "bottom")
ggsave(file.path(run_dir, "Figure2_revision.pdf"), fig2, width = 12, height = 5.5)
message("Saved Figure2_revision.pdf to ", run_dir)
