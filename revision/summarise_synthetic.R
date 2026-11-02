# summarise_synthetic.R
# Combines the synthetic runs (run_synthetic_all.R) and their post-hoc analyses
# (posthoc_synthetic_analysis.R, *_posthoc.rds) into the tables for the paper:
#   A  all runs: NROY size, truth retained, LAI and yield bands, runs and time
#      -> Sections 4.2 and 4.2.5 (discrepancy sensitivity), R1.6 / R2.6
#   B  run time per stage over the six Exp. 1a seeds        -> Section 4.2.4, R2.8
#   C  per-wave table for Exp. 1a                            -> appendix, R2.4
#   D  emulator validation for Exp. 1a                       -> appendix, R2.4
#   E  correlation length sensitivity                        -> appendix, R2.4
# Saves everything to synthetic_tables.rds and one CSV per table.

library(dplyr); library(purrr); library(tidyr); library(tibble)

# Location of the scratch folder: "~/scratch" on the cluster, e.g. "Z:/scratch" on Windows.
if (!exists("scratch")) scratch <- "~/scratch"

run_dir <- paste0(scratch, "/synthetic_revision/")

# true yield at harvest of the synthetic truth (gC m-2): taken from the saved runs if
# synth_obs_feb.RDS stored it, otherwise set it here (it is the black line in Figure 2B)
true_yield_manual <- 347.9858

run_files <- list.files(run_dir, pattern = "^(exp|edcs).*\\.rds$", full.names = TRUE)
run_files <- run_files[!grepl("_posthoc\\.rds$", run_files)]
message("Found ", length(run_files), " runs")

# ---- A: one row per run ---------------------------------------------------------------------
summarise_run <- function(file) {
  res  <- readRDS(file)
  name <- sub("\\.rds$", "", basename(file))
  fin  <- res$data[[length(res$data)]]$training_data
  y    <- fin$training_data$harvest_vector
  lai  <- as.matrix(select(fin$LAI_curves, -date))
  lai_sd <- apply(lai, 1, sd)
  true_yield <- if (!is.null(res$true_yield)) res$true_yield else true_yield_manual

  ph_file <- sub("\\.rds$", "_posthoc.rds", file)
  truth <- if (file.exists(ph_file)) all(readRDS(ph_file)$truth_check$truth_retained) else NA

  tibble(run = name,
         experiment = case_when(!res$settings$use_edcs ~ "Exp. 2 (no EDCs)",
                                res$settings$discs_prop == 0 ~ "Exp. 1b (EDCs, 0%)",
                                TRUE ~ paste0("EDCs, ", 100 * res$settings$discs_prop, "%")),
         discrepancy = res$settings$discs_prop,
         seed = res$settings$seed,
         n_final = length(y),
         proportion_retained = tail(res$wave_table$proportion_retained, 1),
         truth_retained = truth,
         lai_band_width = mean(4 * lai_sd),          # mean daily width of mean +/- 2 SD
         yield_mean = mean(y),
         yield_band_width = 4 * sd(y),
         true_yield_in_band = abs(true_yield - mean(y)) <= 2 * sd(y),
         dalec_evaluations = res$dalec_evaluations_total,
         minutes = res$timings[["total"]] / 60)
}
table_A <- purrr::map(run_files, summarise_run, .progress = "Runs") %>% list_rbind() %>%
  arrange(desc(discrepancy == 0.1 & seed == 1), discrepancy, seed)
print(table_A, width = Inf)

# ---- B: run time per stage, six Exp. 1a seeds ------------------------------------------------
seed_files <- run_files[grepl("^exp1a_md10_seed", basename(run_files))]
table_B <- purrr::map(seed_files, ~ as_tibble(as.list(readRDS(.x)$timings / 60))) %>% list_rbind() %>%
  pivot_longer(everything(), names_to = "stage", values_to = "minutes") %>%
  group_by(stage) %>%
  summarise(mean_minutes = mean(minutes), sd_minutes = sd(minutes), n_runs = n(), .groups = "drop")
print(table_B)

# ---- C, D, E: Exp. 1a (seed 1) ------------------------------------------------------------------
main <- readRDS(file.path(run_dir, "exp1a_md10_seed1.rds"))
table_C <- main$wave_table
table_D <- main$validation
ph_main <- file.path(run_dir, "exp1a_md10_seed1_posthoc.rds")
table_E <- if (file.exists(ph_main)) {
  readRDS(ph_main)$theta_check %>% group_by(theta) %>%
    summarise(mean_abs_std_error = mean(mean_abs_std_error),
              min_coverage_3sd = min(coverage_3sd),
              mean_coverage_3sd = mean(coverage_3sd), .groups = "drop")
} else NULL
print(table_C, width = Inf)
print(table_E)

# share of proposed points that needed a DALEC-Crop run (Section 3.5.2)
screen_share <- with(table_C, sum(pass_implausibility) / sum(points_proposed))
message("Share of proposed points run through DALEC-Crop (waves 1-4): ",
        round(100 * screen_share, 1), "%")

# ---- save ------------------------------------------------------------------------------------------
tables <- list(A_runs = table_A, B_timing = table_B, C_waves = table_C,
               D_validation = table_D, E_theta = table_E, screen_share = screen_share)
saveRDS(tables, file.path(run_dir, "synthetic_tables.rds"))
iwalk(tables[c("A_runs", "B_timing", "C_waves", "D_validation", "E_theta")],
      ~ if (!is.null(.x)) readr::write_csv(.x, file.path(run_dir, paste0("table_", .y, ".csv"))))
