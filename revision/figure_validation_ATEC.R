# figure_validation_ATEC.R
# Emulator validation over all field calibrations (Appendix D.4). Reads only the validation
# metrics saved with each ATEC run; nothing is refitted or rerun.
# Per season, output and wave, the test points of all plots are pooled (20 per plot and
# emulator; the six or four LAI dates pooled as one output), giving:
#   - coverage: share of test points within E +/- 3 SD, with a 95% Jeffreys interval
#   - mean absolute standardised error (about 0.8 for well-calibrated emulators)
# Saves FigureS_validation_ATEC.pdf and table_validation_ATEC.csv in the scratch folder.

library(dplyr); library(tidyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
source(paste0(code_dir, "figure_scheme.R"))

files <- c(list.files(file.path(scratch, "ATEC_revision_2017"), pattern = "plot\\d+\\.rds$", full.names = TRUE),
           list.files(file.path(scratch, "ATEC_revision_2018"), pattern = "plot\\d+\\.rds$", full.names = TRUE))
message("Reading ", length(files), " ATEC files")
val <- bind_rows(lapply(files, function(f) {
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  r$validation %>% mutate(year = r$plot_info$year, plot_id = r$plot_info$plot_id)
}))

output_labels <- c("LAI (all dates)" = "LAI (all dates)", ra_GPP_ratio = "Ra:GPP ratio",
                   yield_GPP_ratio = "Yield:GPP ratio", EDC_SOM1 = "EDC SOM1",
                   EDC_SOM2 = "EDC SOM2", DS_max = "DS max")

val_sum <- val %>%
  mutate(output = ifelse(grepl("^date", output), "LAI (all dates)", output),
         season = ifelse(year == 2017, "2016-2017", "2017-2018"),
         n_in = round(coverage_3sd * n_test)) %>%
  group_by(season, output, wave) %>%
  summarise(n_plots = n_distinct(plot_id), n_points = sum(n_test), n_in = sum(n_in),
            coverage = n_in / n_points,
            mase = weighted.mean(mean_abs_std_error, n_test), .groups = "drop") %>%
  mutate(cov_lo = qbeta(0.025, n_in + 0.5, n_points - n_in + 0.5),
         cov_hi = qbeta(0.975, n_in + 0.5, n_points - n_in + 0.5),
         output = factor(output_labels[output], levels = output_labels))

print(val_sum, n = 60, width = Inf)
readr::write_csv(val_sum, file.path(scratch, "table_validation_ATEC.csv"))

season_cols   <- c("2016-2017" = "#648FFF", "2017-2018" = "#FFB000")
season_shapes <- c("2016-2017" = 16, "2017-2018" = 17)
dodge <- position_dodge(0.4)

p_cov <- ggplot(val_sum, aes(x = factor(wave), y = coverage, colour = season, shape = season)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = cov_lo, ymax = cov_hi), width = 0, linewidth = 0.8, position = dodge) +
  geom_point(size = 3.5, position = dodge) +
  facet_wrap(~ output, nrow = 1, labeller = label_wrap_gen(width = 10)) +
  scale_colour_manual(name = "Season", values = season_cols) +
  scale_shape_manual(name = "Season", values = season_shapes) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Test points within\nE ± 3 SD") +
  theme_paper(base_size = 18)

p_mase <- ggplot(val_sum, aes(x = factor(wave), y = mase, colour = season, shape = season)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40") +
  geom_point(size = 3.5, position = dodge) +
  facet_wrap(~ output, nrow = 1, labeller = label_wrap_gen(width = 10)) +
  scale_colour_manual(name = "Season", values = season_cols) +
  scale_shape_manual(name = "Season", values = season_shapes) +
  labs(x = "Wave", y = "Mean absolute\nstandardised error") +
  theme_paper(base_size = 18) + theme(strip.text = element_blank())

fig <- ggpubr::ggarrange(p_cov, p_mase, ncol = 1, heights = c(1, 1), align = "v",
                         common.legend = TRUE, legend = "bottom")
ggsave(file.path(scratch, "FigureS_validation_ATEC.pdf"), fig, width = 12, height = 8)
message("Saved FigureS_validation_ATEC.pdf and table_validation_ATEC.csv to ", scratch)
