# figure6_yield.R
# Figure 6: yield at harvest for every plot, by nitrogen level (HM vs CARDAMOM vs observed).
# Self-contained: reads only ATEC_comparison_HM_CARDAMOM.rds (from summarise_CARDAMOM.R),
# so it can be sourced or pasted into R on its own; it does not run DALEC-Crop.
# Also sourced by make_figures_ATEC.R.
#
# Usage (from the CARDAMOM working directory):
#   source("synth_experiment/revision/figure6_yield.R")
# Text sizes can be changed without editing the code, e.g.
#   fig6_sizes <- list(base = 20, plot_numbers = 16, strips = 18, legend = 20)

library(dplyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
if (!exists("years") || !is.numeric(years)) years <- c(2017, 2018)
source(paste0(code_dir, "figure_scheme.R"))

# text sizes (points, before the figure is scaled to the page)
if (!exists("fig6_sizes")) fig6_sizes <- list(base = 18, plot_numbers = 15, strips = 17, legend = 19)

comp <- readRDS(file.path(scratch, "ATEC_comparison_HM_CARDAMOM.rds"))
season_titles <- c("2017" = "2016-2017", "2018" = "2017-2018")

# panel titles on two lines, so they stay large without overlapping
fertiliser_labels <- c(
  "1" = "atop(Level~1, (0~kg~N~ha^{-1}))",   "2" = "atop(Level~2, (50~kg~N~ha^{-1}))",
  "3" = "atop(Level~3, (100~kg~N~ha^{-1}))", "4" = "atop(Level~4, (150~kg~N~ha^{-1}))",
  "5" = "atop(Level~5, (200~kg~N~ha^{-1}))")

method_names <- c(HM = "HM", CARDAMOM = "CARDAMOM", ATEC = "Observed (ATEC)")

plot_yield <- function(yr, title, caption = NULL, sz = fig6_sizes) {
  d <- comp$plots_both %>% filter(year == yr)
  if (nrow(d) == 0) stop("No results for ", yr, " in ATEC_comparison_HM_CARDAMOM.rds: ",
                         "run check_existing_results.R and summarise_CARDAMOM.R for that year first")
  obs <- d %>% distinct(plot_id, N_level, obs_yield) %>%
    transmute(plot_id, N_level, method = "ATEC", yield_mean = obs_yield)
  pts <- bind_rows(d %>% select(plot_id, N_level, method, yield_mean, yield_lo, yield_hi), obs) %>%
    mutate(x = as.character(plot_id))

  # "Mean" column per nitrogen level: mean prediction and average band per method;
  # for the observed yields, mean +/- 2 SD across the plots of the level
  means <- pts %>% group_by(N_level, method) %>%
    summarise(yield_lo = if (first(method) == "ATEC") mean(yield_mean) - 2 * sd(yield_mean) else mean(yield_lo),
              yield_hi = if (first(method) == "ATEC") mean(yield_mean) + 2 * sd(yield_mean) else mean(yield_hi),
              yield_mean = mean(yield_mean), .groups = "drop") %>%
    mutate(x = "Mean")
  pts <- bind_rows(pts, means) %>%
    mutate(method = factor(method, levels = c("HM", "CARDAMOM", "ATEC")),
           x = factor(x, levels = c(as.character(sort(unique(d$plot_id))), "Mean")))

  # plots with a CARDAMOM result but no history matching result
  no_hm <- d %>% group_by(plot_id, N_level) %>%
    summarise(has_hm = any(method == "HM"), .groups = "drop") %>%
    filter(!has_hm) %>% mutate(x = factor(as.character(plot_id), levels = levels(pts$x)))

  # shading behind the Mean column
  shade <- pts %>% filter(x != "Mean") %>% group_by(N_level) %>%
    summarise(pos = n_distinct(x) + 1, .groups = "drop")

  ggplot(pts, aes(x = x, y = yield_mean, colour = method)) +
    scale_x_discrete() +
    geom_rect(data = shade, inherit.aes = FALSE,
              aes(xmin = pos - 0.5, xmax = pos + 0.5, ymin = -Inf, ymax = Inf), fill = "grey90") +
    geom_errorbar(data = filter(pts, method != "ATEC" | x == "Mean"),
                  aes(ymin = yield_lo, ymax = yield_hi, linewidth = method),
                  width = 0.5, position = position_dodge(0.6)) +
    geom_point(aes(shape = method, fill = method), colour = "black",
               size = 3, stroke = 0.5, position = position_dodge(0.6)) +
    geom_text(data = no_hm, aes(x = x, y = 680, label = "no NROY\npoints"), inherit.aes = FALSE,
              size = 4, lineheight = 0.85) +
    facet_wrap(~ N_level, scales = "free_x", ncol = 5,
               labeller = labeller(N_level = fertiliser_labels, .default = label_parsed)) +
    scale_colour_methods(labels = method_names) +
    scale_fill_methods(labels = method_names) +
    scale_shape_methods(labels = method_names) +
    scale_linewidth_manual(values = c(HM = 0.9, CARDAMOM = 1.8, ATEC = 0.9), guide = "none") +
    guides(shape = guide_legend(override.aes = list(size = 5))) +   # larger legend symbols
    labs(x = "Plot number", y = expression("Yield at harvest (gC m"^{-2}*")"), title = title,
         caption = caption) +
    coord_cartesian(ylim = c(0, 700)) +
    theme_paper(base_size = sz$base) +
    theme(axis.text.x  = element_text(size = sz$plot_numbers, angle = 90, vjust = 0.5, hjust = 1),
          strip.text   = element_text(size = sz$strips),
          legend.text  = element_text(size = sz$legend),
          legend.title = element_text(size = sz$legend),
          legend.key.size = unit(1.2, "cm"),
          plot.caption = element_text(size = sz$base - 4))
}

# ---- draw and save ------------------------------------------------------------------------
yrs <- as.character(sort(years))
note_mean <- paste("Shaded column: mean over the plots of each nitrogen level",
                   "(methods: mean prediction and average band; observed: mean +/- 2 SD across plots).")
note_nroy <- "\"no NROY points\": history matching found no parameter sets consistent with the observations (Section 5.3)."
fig6_notes <- if (length(yrs) == 1) {
  setNames(list(paste(note_mean, if (yrs == "2017") note_nroy)), yrs)
} else {
  setNames(lapply(yrs, function(y) if (y == "2017") note_nroy else if (y == tail(yrs, 1)) note_mean else NULL), yrs)
}
fig6_panels <- lapply(yrs, function(y) plot_yield(as.integer(y), season_titles[[y]], caption = fig6_notes[[y]]))
fig6 <- if (length(fig6_panels) == 1) fig6_panels[[1]] else
  ggpubr::ggarrange(plotlist = fig6_panels, nrow = length(fig6_panels), labels = LETTERS[seq_along(yrs)],
                    font.label = list(size = fig6_sizes$base + 2, face = "bold"),
                    common.legend = TRUE, legend = "bottom")
fig6_file <- file.path(scratch, paste0("Figure6_revision_", paste(yrs, collapse = "_"), ".pdf"))
ggsave(fig6_file, fig6, width = 14, height = 6.5 * length(yrs))
message("Saved ", fig6_file)
