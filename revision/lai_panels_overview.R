# lai_panels_overview.R
# Draws Figure 5-style LAI panels (history matching vs CARDAMOM) for any set of plots, to
# compare candidates for Figure 5. By default: a random subset of two plots per nitrogen
# level and season (plots calibrated by both methods). Each panel title shows the plot's LAI RMSE and yield
# error per method, so typical and unusual plots are easy to spot.
#
# Run from the CARDAMOM working directory on the cluster (CARDAMOM's samples are run
# through DALEC-Crop: a few seconds to a minute per plot).
#
# Usage:
#   source("synth_experiment/revision/lai_panels_overview.R")              # 2 random plots per level
#   plots_to_draw <- list("2017" = c(22, 34, 39), "2018" = c(3, 37))      # or choose plots
#   source("synth_experiment/revision/lai_panels_overview.R")

library(dplyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"

source(paste0(code_dir, "figure5_LAI.R"))   # lai_panel(), season settings, colour scheme

comp <- readRDS(file.path(scratch, "ATEC_comparison_HM_CARDAMOM.rds"))
d <- comp$plots_both %>% filter(in_both)

# default selection: a random subset of plots from every nitrogen level, per season
# (plots calibrated by both methods only). Set per_level or selection_seed before sourcing
# to change it, or set plots_to_draw directly.
if (!exists("per_level"))      per_level      <- 2
if (!exists("selection_seed")) selection_seed <- 1
if (!exists("plots_to_draw")) {
  set.seed(selection_seed)
  sel <- d %>% distinct(year, N_level, plot_id) %>%
    group_by(year, N_level) %>% slice_sample(n = per_level) %>% ungroup() %>%
    arrange(year, N_level, plot_id)
  plots_to_draw <- lapply(split(sel$plot_id, sel$year), identity)
}
level_of <- d %>% distinct(year, plot_id, N_level)

# per-plot errors for the panel titles
errors <- d %>%
  mutate(yield_err = yield_mean - obs_yield) %>%
  select(year, plot_id, method, lai_rmse, yield_err)

panel_title <- function(plot_id, year) {
  e <- errors %>% filter(year == !!year, plot_id == !!plot_id)
  fmt <- function(m) {
    r <- e[e$method == m, ]
    sprintf("%s: LAI RMSE %.2f, yield error %+.0f", m, r$lai_rmse, r$yield_err)
  }
  lev <- level_of$N_level[level_of$year == year & level_of$plot_id == plot_id][1]
  rate <- c("1" = 0, "2" = 50, "3" = 100, "4" = 150, "5" = 200)[as.character(lev)]
  # three short lines, so the title can be large without running off the panel
  paste0(ifelse(year == 2017, "2016-2017", "2017-2018"), ", Level ", lev, " (", rate, " kg N/ha), plot ", plot_id, "\n",
         fmt("HM"), "\n", fmt("CARDAMOM"))
}

panels <- list()
for (yr in names(plots_to_draw)) {
  for (pid in plots_to_draw[[yr]]) {
    message("Drawing plot ", pid, " (", yr, ")")
    # bands come from <scratch>/figure5_cache/ when available, so redrawing does not rerun DALEC-Crop
    p <- tryCatch(lai_panel(pid, as.integer(yr), title = panel_title(pid, as.integer(yr)), base_size = 18),
                  error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
    if (!is.null(p)) panels[[paste(yr, pid)]] <- p + theme(plot.title = element_text(size = 14, lineheight = 1.0))
  }
}

# four panels per page, one PDF. The pages are wide and short (about 1.6:1), so two of them
# stacked fit on one page of the paper (8 plots per page). Change page_width / page_height
# before sourcing to adjust.
if (!exists("page_width"))  page_width  <- 12
if (!exists("page_height")) page_height <- 7.5
out_file <- file.path(scratch, "LAI_panels_overview.pdf")
pdf(out_file, width = page_width, height = page_height)
for (i in seq(1, length(panels), by = 4)) {
  page <- panels[i:min(i + 3, length(panels))]
  print(ggpubr::ggarrange(plotlist = page, ncol = 2, nrow = 2,
                          common.legend = TRUE, legend = "right"))
}
dev.off()
message("Saved ", length(panels), " panels to ", out_file)

# the same errors as a table, sorted within each season
print(errors %>% filter(plot_id %in% unlist(plots_to_draw)) %>%
        tidyr::pivot_wider(names_from = method, values_from = c(lai_rmse, yield_err)) %>%
        arrange(year, lai_rmse_HM), width = Inf)
