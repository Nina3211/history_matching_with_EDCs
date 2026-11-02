# make_figures_ATEC.R
# Creates the field-data figures from the reruns:
#   Figure 5: LAI curves for one plot per season (HM vs CARDAMOM)
#   Figure 6: yield at harvest for every plot, by nitrogen level (HM vs CARDAMOM vs observed)
# Colours, shapes and line types come from figure_scheme.R (colour-blind safe).
# Saves PDFs to the scratch folder. Only seasons listed in `years` are drawn, so this can
# be run for 2018 now and for both seasons once 2017 is complete.
#
# Before running:
#   - check_existing_results.R and summarise_CARDAMOM.R have been run for the same years
#     (Figure 6 uses their combined table ATEC_comparison_HM_CARDAMOM.rds)
#   - run from the CARDAMOM working directory on the cluster (Figure 5 runs DALEC-Crop
#     for CARDAMOM's posterior samples)
#
# Usage:  years <- 2018;  source("synth_experiment/revision/make_figures_ATEC.R")

library(dplyr); library(tidyr); library(ggplot2)

if (!exists("scratch"))  scratch  <- "~/scratch"
if (!exists("code_dir")) code_dir <- "synth_experiment/revision/"
# "years" is also a function in the chron package (loaded with CARDAMOM): only use a numeric setting
if (!exists("years") || !is.numeric(years)) years    <- c(2017, 2018)

# plots shown in Figure 5 (as in the paper: plot 22 in 2016-2017, plot 3 in 2017-2018)
figure5_plots <- c("2017" = 22, "2018" = 3)
season_titles <- c("2017" = "2016-2017", "2018" = "2017-2018")

# lai_panel(), the season settings and the shared colour scheme (figure_scheme.R)
source(paste0(code_dir, "figure5_LAI.R"))

# ---- Figure 6: yield at harvest (figure6_yield.R draws and saves it) --------------------
message("Figure 6 (yield)")
source(paste0(code_dir, "figure6_yield.R"))
panel_labels <- LETTERS[seq_along(yrs)]


message("Figure 5 (LAI): the first run per plot runs DALEC-Crop; later runs use figure5_cache/")
fig5_panels <- lapply(yrs, function(y)
  lai_panel(figure5_plots[[y]], as.integer(y),
            title = paste0("Plot ", figure5_plots[[y]], " (", season_titles[[y]], ")")))
fig5 <- if (length(fig5_panels) == 1) fig5_panels[[1]] else
  ggpubr::ggarrange(plotlist = fig5_panels, ncol = length(fig5_panels), labels = panel_labels,
                    common.legend = TRUE, legend = "bottom")
ggsave(file.path(scratch, paste0("Figure5_revision_", paste(yrs, collapse = "_"), ".pdf")),
       fig5, width = 6.5 * length(yrs), height = 6)

message("Saved Figure5_revision_", paste(yrs, collapse = "_"), ".pdf and Figure6_revision_",
        paste(yrs, collapse = "_"), ".pdf in ", scratch)
