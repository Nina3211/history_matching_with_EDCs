# figure_scheme.R
# One colour / shape / line-type scheme for every figure in the paper.
# Colours: IBM colour-blind-safe palette (distinguishable with the common types of colour
# vision deficiency). Shapes and line types are added as a second cue, so the figures
# remain readable in greyscale.
#
# Use in any ggplot:
#   source(paste0(code_dir, "figure_scheme.R"))
#   ggplot(...) + scale_colour_methods() + scale_fill_methods() +
#     scale_shape_methods() + scale_linetype_methods() + theme_paper()
# The method names in the data must match the names below.

library(ggplot2)

method_colours <- c(
  "CARDAMOM"           = "#FFB000",   # amber
  "HM"                 = "#648FFF",   # blue
  "HM 1a: EDCs and md" = "#648FFF",   # blue (the main history matching run)
  "HM 1b: EDCs"        = "#785EF0",   # purple
  "HM 2: md"           = "#DC267F",   # magenta
  "ATEC"               = "black",     # observations
  "Observations"       = "black",
  "Truth"              = "black",
  "Prior"              = "grey45")

method_shapes <- c(
  "CARDAMOM" = 24, "HM" = 21, "HM 1a: EDCs and md" = 21, "HM 1b: EDCs" = 22,
  "HM 2: md" = 25, "ATEC" = 23, "Observations" = 16, "Truth" = 23, "Prior" = 4)

method_linetypes <- c(
  "CARDAMOM" = "dashed", "HM" = "solid", "HM 1a: EDCs and md" = "solid",
  "HM 1b: EDCs" = "longdash", "HM 2: md" = "dotdash",
  "ATEC" = "solid", "Observations" = "solid", "Truth" = "solid", "Prior" = "dotted")

scale_colour_methods   <- function(name = "Method", ...) scale_colour_manual(name = name, values = method_colours, ...)
scale_fill_methods     <- function(name = "Method", ...) scale_fill_manual(name = name, values = method_colours, ...)
scale_shape_methods    <- function(name = "Method", ...) scale_shape_manual(name = name, values = method_shapes, ...)
scale_linetype_methods <- function(name = "Method", ...) scale_linetype_manual(name = name, values = method_linetypes, ...)

# larger text throughout (Reviewer 2: axis, legend and panel labels were hard to read)
theme_paper <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.key.width = unit(1.5, "cm"),
          strip.background = element_rect(fill = "grey92"),
          plot.title = element_text(face = "bold"))
}
