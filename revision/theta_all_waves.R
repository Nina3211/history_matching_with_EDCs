# theta_all_waves.R
# Correlation length sensitivity (R2.4) at every wave: refits each wave's emulators on that
# wave's saved training data with theta = 0.4 / 0.6 / 0.8 and validates them on that wave's
# saved test points. Only emulators are refitted; no DALEC-Crop or history matching runs.

library(hmer); library(dplyr); library(purrr); library(tibble)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("run_file")) run_file <- file.path(scratch, "synthetic_revision", "exp1a_md10_seed1.rds")
if (!exists("thetas"))   thetas   <- c(0.4, 0.6, 0.8)
# FALSE (default): validate on each wave's 20 independent test points.
# TRUE: also use the next wave's simulator runs. These were selected by the emulators being
# validated (biased in their favour, and cannot reveal wrongly ruled-out points), so they
# are only an informal accuracy check inside the retained region, not for the paper.
if (!exists("extended_test")) extended_test <- FALSE

res <- readRDS(run_file)
targets    <- res$settings$targets
discs_prop <- res$settings$discs_prop
n_waves    <- length(res$emulators)
output_names <- names(res$emulators[[n_waves]])        # Wave 1 emulators (last in the list)
model_discs  <- purrr::map(output_names, ~ list(internal = discs_prop * targets[[.x]]$val,
                                                external = NULL))

theta_waves <- bind_rows(lapply(seq_len(n_waves), function(w) {
  train <- res$data[[w]]$training_data
  test  <- res$data[[w]]$testing_data
  if (extended_test) {
    # also validate on the next wave's simulator runs, which these emulators have not seen
    nxt <- res$data[[w + 1]]
    nxt_train <- if (is.data.frame(nxt$training_data)) nxt$training_data else nxt$training_data$training_data
    cols <- intersect(names(test), names(nxt_train))
    test <- bind_rows(test[, cols], nxt_train[, cols],
                      if (!is.null(nxt$testing_data)) nxt$testing_data[, cols])
  }
  bind_rows(lapply(thetas, function(th) {
    message("Wave ", w, ", theta = ", th)
    ems <- emulator_from_data(
      input_data = train, output_names = output_names,
      specified_priors = list(hyper_p = rep(th, length(output_names))),
      ranges = res$ranges[[w]], beta.var = TRUE, check.ranges = FALSE,
      targets = targets[output_names], add = TRUE, order = 2,
      discrepancies = model_discs, verbose = FALSE)
    bind_rows(lapply(ems, function(em) {
      x  <- test[, names(em$ranges), drop = FALSE]
      se <- (test[[em$output_name]] - em$get_exp(x)) / sqrt(em$get_cov(x))
      tibble(wave = w, theta = th, output = em$output_name, n_test = length(se),
             mean_abs_std_error = mean(abs(se)), coverage_3sd = mean(abs(se) <= 3))
    }))
  }))
}))

# why theta matters little: distance from each test point to its nearest training point,
# in the scaled active parameters of each emulator, and the implied largest correlation
nn_distances <- bind_rows(lapply(seq_len(n_waves), function(w) {
  train <- res$data[[w]]$training_data; test <- res$data[[w]]$testing_data
  bind_rows(lapply(res$emulators[[n_waves - w + 1]], function(em) {
    act <- names(em$ranges)[em$active_vars]
    scale01 <- function(df) as.matrix(as.data.frame(lapply(act, function(p)
      2 * (df[[p]] - em$ranges[[p]][1]) / diff(em$ranges[[p]]) - 1)))
    tr <- scale01(train); te <- scale01(test)
    d <- apply(te, 1, function(z) min(sqrt(colSums((t(tr) - z)^2))))
    tibble(wave = w, output = em$output_name, n_active = length(act),
           median_nn_distance = median(d),
           max_corr_theta_0.8 = exp(-(min(d) / 0.8)^2))
  }))
}))
print(nn_distances %>% group_by(wave) %>%
        summarise(median_active = median(n_active),
                  median_nn_distance = median(median_nn_distance),
                  largest_correlation_theta_0.8 = max(max_corr_theta_0.8), .groups = "drop"))

theta_summary <- theta_waves %>% group_by(wave, theta) %>%
  summarise(n_test = max(n_test),
            mean_abs_std_error = mean(mean_abs_std_error),
            min_coverage_3sd = min(coverage_3sd),
            mean_coverage_3sd = mean(coverage_3sd), .groups = "drop")
print(theta_summary, n = 30)

out <- sub("\\.rds$", "_theta_all_waves.csv", run_file)
readr::write_csv(theta_waves, out)
message("Saved per-output results to ", out)
