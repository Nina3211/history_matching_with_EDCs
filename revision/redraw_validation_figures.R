# redraw_validation_figures.R
# Redraws Figures 12 and 13 (emulator validation) with larger text, from the CSV files the
# earlier scripts saved. Nothing is refitted or rerun:
#   Figure 12 (synthetic Exp. 1a): <scratch>/synthetic_revision/validation_points_synthetic.csv
#                                   (written by validation_and_stopping.R)
#   Figure 13 (field data):        <scratch>/table_validation_ATEC.csv
#                                   (written by figure_validation_ATEC.R)
# If a CSV is missing, it is recomputed from the saved results (emulators or validation logs),
# which takes a minute or two but still runs no calibration or DALEC-Crop.
#
# Usage (from the CARDAMOM working directory):
#   source("synth_experiment/revision/redraw_validation_figures.R")
# Text sizes can be changed without editing the code, e.g.
#   val_sizes <- list(base = 20, strips = 15, axis = 12); source(...)

library(dplyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
source(paste0(code_dir, "figure_scheme.R"))

if (!exists("val_sizes")) val_sizes <- list(base = 18, strips = 14, axis = 11)
sz <- val_sizes
syn_dir <- file.path(scratch, "synthetic_revision")

# ---- Figure 12: synthetic emulator validation, one page per wave --------------------------
val_csv <- file.path(syn_dir, "validation_points_synthetic.csv")
if (file.exists(val_csv)) {
  val_points <- readr::read_csv(val_csv, show_col_types = FALSE)
} else {
  message("validation_points_synthetic.csv not found: recomputing from the saved emulators (no reruns)")
  library(hmer)
  main <- readRDS(file.path(syn_dir, "exp1a_md10_seed1.rds"))
  ems <- main$emulators; n_waves <- length(ems)
  val_points <- bind_rows(lapply(seq_len(n_waves), function(w) {
    test <- main$data[[w]]$testing_data
    bind_rows(lapply(ems[[n_waves - w + 1]], function(em) {
      x <- test[, names(em$ranges), drop = FALSE]
      tibble(wave = w, output = em$output_name, simulator = test[[em$output_name]],
             expectation = as.numeric(em$get_exp(x)), sd = sqrt(as.numeric(em$get_cov(x))))
    }))
  })) %>% mutate(within_3sd = ifelse(abs(simulator - expectation) <= 3 * sd, "Within 3 SD", "Outside 3 SD"))
  readr::write_csv(val_points, val_csv)
}
output_order <- c(sort(unique(grep("^date", val_points$output, value = TRUE))),
                  "EDC_SOM1", "EDC_SOM2", "ra_GPP_ratio", "yield_GPP_ratio", "DS_max")
val_points <- val_points %>%
  mutate(output = factor(output, levels = intersect(output_order, unique(output))),
         within_3sd = factor(within_3sd, levels = c("Within 3 SD", "Outside 3 SD")))

val_plot <- function(w) {
  ggplot(filter(val_points, wave == w), aes(x = simulator, y = expectation)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50") +
    geom_errorbar(aes(ymin = expectation - 3 * sd, ymax = expectation + 3 * sd, colour = within_3sd),
                  width = 0, linewidth = 0.6) +
    geom_point(aes(colour = within_3sd, shape = within_3sd), size = 2.5) +
    facet_wrap(~ output, scales = "free", ncol = 4) +
    scale_colour_manual(name = NULL, values = c("Within 3 SD" = "#648FFF", "Outside 3 SD" = "#DC267F"),
                        drop = FALSE) +
    scale_shape_manual(name = NULL, values = c("Within 3 SD" = 16, "Outside 3 SD" = 17), drop = FALSE) +
    scale_x_continuous(n.breaks = 4) + scale_y_continuous(n.breaks = 4) +
    labs(x = "DALEC-Crop output at test point", y = "Emulator expectation ± 3 SD",
         title = paste("Wave", w, "emulators (synthetic Exp. 1a)")) +
    theme_paper(base_size = sz$base) +
    theme(strip.text = element_text(size = sz$strips), axis.text = element_text(size = sz$axis),
          legend.text = element_text(size = sz$base))
}
waves <- sort(unique(val_points$wave))
pdf(file.path(syn_dir, "FigureS_validation_synthetic.pdf"), width = 11, height = 9)
for (w in waves) print(val_plot(w))
dev.off()
# one file per wave as well, for the sub-figures (a), (b), ... in LaTeX
for (w in waves) ggsave(file.path(syn_dir, paste0("FigureS_validation_synthetic_wave", w, ".pdf")),
                        val_plot(w), width = 11, height = 9)
message("Saved Figure 12 (FigureS_validation_synthetic*.pdf) to ", syn_dir)

# ---- Figure 13: field-data emulator validation --------------------------------------------
atec_csv <- file.path(scratch, "table_validation_ATEC.csv")
if (!file.exists(atec_csv)) {
  message("table_validation_ATEC.csv not found: recomputing it with figure_validation_ATEC.R (reads the saved logs only)")
  source(paste0(code_dir, "figure_validation_ATEC.R"))
}
output_labels <- c("LAI (all dates)", "Ra:GPP ratio", "Yield:GPP ratio", "EDC SOM1", "EDC SOM2", "DS max")
val_sum <- readr::read_csv(atec_csv, show_col_types = FALSE) %>%
  mutate(output = factor(output, levels = output_labels))

season_cols   <- c("2016-2017" = "#648FFF", "2017-2018" = "#FFB000")
season_shapes <- c("2016-2017" = 16, "2017-2018" = 17)
dodge <- position_dodge(0.5)
atec_theme <- theme_paper(base_size = sz$base) +
  theme(strip.text = element_text(size = sz$strips), axis.text = element_text(size = sz$axis + 1),
        legend.text = element_text(size = sz$base), legend.title = element_text(size = sz$base))

p_cov <- ggplot(val_sum, aes(x = factor(wave), y = coverage, colour = season, shape = season)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = cov_lo, ymax = cov_hi), width = 0, linewidth = 0.9, position = dodge) +
  geom_point(size = 3.5, position = dodge) +
  facet_wrap(~ output, nrow = 1, labeller = label_wrap_gen(width = 10)) +
  scale_colour_manual(name = "Season", values = season_cols) +
  scale_shape_manual(name = "Season", values = season_shapes) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Test points within\nE ± 3 SD") +
  atec_theme

p_mase <- ggplot(val_sum, aes(x = factor(wave), y = mase, colour = season, shape = season)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40") +
  geom_point(size = 3.5, position = dodge) +
  facet_wrap(~ output, nrow = 1) +
  scale_colour_manual(name = "Season", values = season_cols) +
  scale_shape_manual(name = "Season", values = season_shapes) +
  labs(x = "Wave", y = "Mean absolute\nstandardised error") +
  atec_theme + theme(strip.text = element_blank())

fig13 <- ggpubr::ggarrange(p_cov, p_mase, ncol = 1, heights = c(1, 1), align = "v",
                           common.legend = TRUE, legend = "bottom")
ggsave(file.path(scratch, "FigureS_validation_ATEC.pdf"), fig13, width = 12, height = 8)
message("Saved Figure 13 (FigureS_validation_ATEC.pdf) to ", scratch)
