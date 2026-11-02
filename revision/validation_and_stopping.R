# validation_and_stopping.R
# Emulator validation and the stopping rule (R2.4), from the saved runs only: no
# simulator or history matching runs are repeated.
#   1. Synthetic Exp. 1a: emulator prediction +/- 3 SD vs simulator output at the test
#      points, every target and wave (uses the saved emulators and test data)
#   2. ATEC: validation metrics of every plot (saved during the runs), summarised over plots
#   3. Stopping rule: share of the parameter space retained and emulator variance ratio
#      per wave, for the synthetic runs and the ATEC plots, plus a table for the text
# Saves PDFs and CSVs to the scratch folder.

library(dplyr); library(tidyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
source(paste0(code_dir, "figure_scheme.R"))
library(hmer)   # emulator objects need hmer loaded to be used

syn_dir <- file.path(scratch, "synthetic_revision")
out_dir <- syn_dir

# ---- 1. synthetic validation plots ------------------------------------------------------------
main <- readRDS(file.path(syn_dir, "exp1a_md10_seed1.rds"))
ems  <- main$emulators                         # newest wave first
n_waves <- length(ems)

val_points <- bind_rows(lapply(seq_len(n_waves), function(w) {
  test <- main$data[[w]]$testing_data         # Wave w emulators were validated on these
  bind_rows(lapply(ems[[n_waves - w + 1]], function(em) {
    x <- test[, names(em$ranges), drop = FALSE]
    tibble(wave = w, output = em$output_name,
           simulator = test[[em$output_name]],
           expectation = as.numeric(em$get_exp(x)), sd = sqrt(as.numeric(em$get_cov(x))))
  }))
})) %>% mutate(within_3sd = abs(simulator - expectation) <= 3 * sd,
               within_3sd = factor(ifelse(within_3sd, "Within 3 SD", "Outside 3 SD"),
                                   levels = c("Within 3 SD", "Outside 3 SD")))

val_plot <- function(w) {
  ggplot(filter(val_points, wave == w), aes(x = simulator, y = expectation)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50") +
    geom_errorbar(aes(ymin = expectation - 3 * sd, ymax = expectation + 3 * sd, colour = within_3sd),
                  width = 0) +
    geom_point(aes(colour = within_3sd, shape = within_3sd), size = 2) +
    facet_wrap(~ output, scales = "free", ncol = 4) +
    scale_colour_manual(name = NULL, values = c("Within 3 SD" = "#648FFF", "Outside 3 SD" = "#DC267F")) +
    scale_shape_manual(name = NULL, values = c("Within 3 SD" = 16, "Outside 3 SD" = 17)) +
    labs(x = "DALEC-Crop output at test point", y = "Emulator expectation +/- 3 SD",
         title = paste("Wave", w, "emulators (synthetic Exp. 1a)")) +
    theme_paper(base_size = 11)
}
pdf(file.path(out_dir, "FigureS_validation_synthetic.pdf"), width = 13, height = 9)
for (w in seq_len(n_waves)) print(val_plot(w))
dev.off()

# ---- 2. ATEC validation metrics over plots ----------------------------------------------------
atec_files <- c(list.files(file.path(scratch, "ATEC_revision_2017"), pattern = "plot\\d+\\.rds$", full.names = TRUE),
                list.files(file.path(scratch, "ATEC_revision_2018"), pattern = "plot\\d+\\.rds$", full.names = TRUE))
message("Reading ", length(atec_files), " ATEC files")
atec_logs <- lapply(atec_files, function(f) {
  res <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  list(validation = res$validation %>% mutate(year = res$plot_info$year, plot_id = res$plot_info$plot_id),
       waves = res$wave_table %>% mutate(year = res$plot_info$year, plot_id = res$plot_info$plot_id))
})
atec_logs <- Filter(Negate(is.null), atec_logs)
atec_val   <- bind_rows(lapply(atec_logs, `[[`, "validation"))
atec_waves <- bind_rows(lapply(atec_logs, `[[`, "waves"))

atec_val_plot <- atec_val %>%
  mutate(output = ifelse(grepl("^date", output), "LAI (all dates)", output),
         season = ifelse(year == 2017, "2016-2017", "2017-2018")) %>%
  ggplot(aes(x = factor(wave), y = coverage_3sd, fill = season)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~ output, ncol = 3) +
  scale_fill_manual(name = "Season", values = c("2016-2017" = "#648FFF", "2017-2018" = "#FFB000")) +
  labs(x = "Wave", y = "Share of test points within E +/- 3 SD",
       caption = "One value per plot, emulator and wave (20 test points each); dashed line: 0.95") +
  theme_paper(base_size = 12)
ggsave(file.path(out_dir, "FigureS_validation_ATEC.pdf"), atec_val_plot, width = 12, height = 9)

# ---- 3. stopping rule ---------------------------------------------------------------------------
syn_files <- list.files(syn_dir, pattern = "^(exp1a|exp1b|edcs).*seed\\d\\.rds$", full.names = TRUE)
syn_waves <- bind_rows(lapply(syn_files, function(f) {
  r <- readRDS(f)
  r$wave_table %>% mutate(run = sub("\\.rds$", "", basename(f)))
}))

stop_long <- bind_rows(
  syn_waves %>% transmute(source = "Synthetic runs", wave, proportion_retained, max_variance_ratio),
  atec_waves %>% transmute(source = ifelse(year == 2017, "ATEC 2016-2017", "ATEC 2017-2018"),
                           wave, proportion_retained, max_variance_ratio)) %>%
  pivot_longer(c(proportion_retained, max_variance_ratio), names_to = "measure") %>%
  mutate(measure = recode(measure,
                          proportion_retained = "Share of parameter space retained",
                          max_variance_ratio = "Largest emulator variance / (obs. + discrepancy variance)"))

stop_plot <- ggplot(stop_long, aes(x = factor(wave), y = value, fill = source)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  geom_hline(data = tibble(measure = "Largest emulator variance / (obs. + discrepancy variance)", y = 1),
             aes(yintercept = y), linetype = "dashed", colour = "grey40") +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  scale_y_log10() +
  scale_fill_manual(name = NULL, values = c("Synthetic runs" = "#785EF0",
                                            "ATEC 2016-2017" = "#648FFF", "ATEC 2017-2018" = "#FFB000")) +
  labs(x = "Wave", y = NULL) +
  theme_paper(base_size = 12)
ggsave(file.path(out_dir, "FigureS_stopping_rule.pdf"), stop_plot, width = 13, height = 6)

# table for the text: how much each wave still reduces the space, and the variance ratio
stop_table <- bind_rows(
  syn_waves %>% filter(grepl("^exp1a", run)) %>% mutate(source = "Synthetic Exp. 1a (6 seeds)"),
  atec_waves %>% mutate(source = ifelse(year == 2017, "ATEC 2016-2017", "ATEC 2017-2018"))) %>%
  group_by(source, wave) %>%
  summarise(median_proportion_retained = median(proportion_retained),
            median_max_variance_ratio = median(max_variance_ratio),
            share_runs_ratio_below_1 = mean(max_variance_ratio < 1),
            .groups = "drop") %>%
  group_by(source) %>% arrange(wave, .by_group = TRUE) %>%
  mutate(reduction_vs_previous_wave = 1 - median_proportion_retained / lag(median_proportion_retained)) %>%
  ungroup()
print(stop_table, n = 20, width = Inf)
readr::write_csv(stop_table, file.path(out_dir, "table_stopping_rule.csv"))
readr::write_csv(val_points, file.path(out_dir, "validation_points_synthetic.csv"))
