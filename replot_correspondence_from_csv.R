#!/usr/bin/env Rscript
# =============================================================================
# Re-draw the cluster-correspondence figures from the CSVs written by
# cluster_correspondence_labeltransfer.R, without repeating the label transfer.
#
# Produces three orderings of the same matrix:
#   (a) fig_correspondence_<CLASS>                  rows AND columns seriated
#   (b) fig_correspondence_byClusterNumber_<CLASS>  rows AND columns native
#   (c) fig_correspondence_rowNumberColSorted_<CLASS>
#          rows in native cross-condition cluster order (0, 1, 2, ...),
#          columns reordered so the best matches fall on a diagonal.
#       This is the one to use when the y axis has to stay readable against
#       the manuscript's cluster numbering but the structure should still
#       be visible.
#
# plus sourcedata_fig_correspondence_<CLASS>.csv, the long-form table of every
# value that appears in the panels.
#
# Usage:
#   OUTDIR=correspondence_out CLASS=Ex  Rscript replot_correspondence_from_csv.R
#   OUTDIR=correspondence_out CLASS=Inh Rscript replot_correspondence_from_csv.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})

OUTDIR <- Sys.getenv("OUTDIR", "correspondence_out")
CLASS  <- Sys.getenv("CLASS",  "Ex")
CLASS_LONG <- c(Ex = "excitatory (vGluT2+)", Inh = "inhibitory (vGat+)")[CLASS]
if (is.na(CLASS_LONG)) CLASS_LONG <- CLASS

TH_FRAC_PRIMARY <- 0.50
TH_SCORE_MIN    <- 0.40

f <- function(stem, ext = "csv") file.path(OUTDIR, sprintf("%s_%s.%s", stem, CLASS, ext))

M_fwd <- as.matrix(read.csv(f("correspondence_matrix_forward"),
                            row.names = 1, check.names = FALSE))
M_rev <- as.matrix(read.csv(f("correspondence_matrix_reverse"),
                            row.names = 1, check.names = FALSE))
tab   <- read.csv(f("correspondence_table"), check.names = FALSE,
                  colClasses = c(query_cluster = "character"))

stopifnot(nrow(M_fwd) == nrow(tab))

best_fwd <- colnames(M_fwd)[apply(M_fwd, 1, which.max)]; names(best_fwd) <- rownames(M_fwd)
best_rev <- colnames(M_rev)[apply(M_rev, 1, which.max)]; names(best_rev) <- rownames(M_rev)

# --- tidy tile table ---------------------------------------------------------
df0 <- as.data.frame(as.table(M_fwd))
colnames(df0) <- c("query_cluster", "ref_label", "frac")
df0$query_cluster <- as.character(df0$query_cluster)
df0$ref_label     <- as.character(df0$ref_label)

key     <- paste(df0$query_cluster, df0$ref_label)
fwd_hit <- paste(rownames(M_fwd), best_fwd)
rev_hit <- paste(best_rev, rownames(M_rev))
df0$mark <- dplyr::case_when(
  key %in% fwd_hit & key %in% rev_hit ~ "reciprocal best",
  key %in% fwd_hit ~ "one-directional best",
  key %in% rev_hit ~ "one-directional best",
  TRUE                                ~ NA_character_)

# legend order: reciprocal first, then one-directional (forward-only or reverse-only)
MARK_LEVELS <- c("reciprocal best", "one-directional best")
df0$mark <- factor(df0$mark, levels = MARK_LEVELS)

bar_df <- data.frame(
  query_cluster = tab$query_cluster,
  n_cells       = tab$n_cells,
  med_score     = tab$median_pred_score,
  frac_best     = tab$frac_best,
  stringsAsFactors = FALSE)

# --- source data -------------------------------------------------------------
sourcedata <- rbind(
  data.frame(panel = "A (heatmap)", query_cluster = df0$query_cluster,
             ref_label = df0$ref_label, value_type = "fraction_of_query_cluster",
             value = df0$frac, mark = df0$mark, stringsAsFactors = FALSE),
  data.frame(panel = "B (bar)", query_cluster = bar_df$query_cluster,
             ref_label = NA_character_, value_type = "median_prediction_score",
             value = bar_df$med_score, mark = NA_character_),
  data.frame(panel = "C (bar)", query_cluster = bar_df$query_cluster,
             ref_label = NA_character_, value_type = "fraction_in_best_matching_label",
             value = bar_df$frac_best, mark = NA_character_),
  data.frame(panel = "n", query_cluster = bar_df$query_cluster,
             ref_label = NA_character_, value_type = "n_cells_in_query_cluster",
             value = as.numeric(bar_df$n_cells), mark = NA_character_))
write.csv(sourcedata, f("sourcedata_fig_correspondence"), row.names = FALSE)

# --- plotting ----------------------------------------------------------------
# Figures are drawn at their final printed size so that the point sizes set
# here are the point sizes that appear on the page; nothing is rescaled when
# the panel is placed into the assembled figure. TILE_MM is the centre-to-centre
# tile spacing, held constant across all four correspondence panels so they form
# a visually consistent set.
TICK_EVERY <- as.integer(Sys.getenv("TICK_EVERY", "10"))
TILE_MM    <- as.numeric(Sys.getenv("TILE_MM", "1.283"))
BASE_PT    <- as.numeric(Sys.getenv("BASE_PT", "7"))   # Nature Portfolio: 5-7 pt body text
YLAB_PT    <- as.numeric(Sys.getenv("YLAB_PT", "4"))   # y-axis cluster numbers
# Margins reserved for axis labels, titles and legend (mm). These depend on the
# text, not on the figure size, so they were measured once from a rendered
# panel by calibrate_panel_size.py and fixed here; with them the tile pitch
# comes out at TILE_MM for every figure.
.cal <- list(Ex = c(37.4, 11.4), Inh = c(37.4, 11.3))[[CLASS]]
MARGIN_W_MM <- as.numeric(Sys.getenv("MARGIN_W_MM", .cal[1]))
MARGIN_H_MM <- as.numeric(Sys.getenv("MARGIN_H_MM", .cal[2]))
BAR_FRAC    <- as.numeric(Sys.getenv("BAR_FRAC", "0.22"))  # bar height / heatmap

# All text, tick marks and panel borders are pure black (theme_bw defaults to
# dark greys); tick marks are 0.3 mm long at the final printed size.
base_theme <- theme_bw(base_size = BASE_PT) +
  theme(panel.grid = element_blank(),
        text = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.15),
        axis.ticks.length = unit(0.3, "mm"),
        axis.line = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(3.2, "mm"),
        legend.key.width  = unit(3.2, "mm"),
        legend.spacing.y  = unit(1, "mm"),
        legend.text  = element_text(size = BASE_PT, colour = "black"),
        legend.title = element_text(size = BASE_PT, colour = "black"),
        legend.margin = margin(0, 0, 0, 1, "mm"),
        plot.margin = margin(1, 1, 1, 1, "mm"))

# ---------------------------------------------------------------------------
# Panel layout, matching the arrangement used in the manuscript figures:
#   x = query cluster in native numbering, ticked every TICK_EVERY;
#   y = cross-species subtype, ordered so the best matches run down the
#       diagonal, every label shown;
#   a confidence bar aligned above the heatmap, sharing the x axis;
#   no tile borders and no threshold lines -- the cut-offs are given in the
#   Methods and drawing them here would imply they are standard values.
# ---------------------------------------------------------------------------

make_fig <- function(x_ord, y_ord, conf_var, conf_lab, y_title) {
  d <- df0
  d$query_cluster <- factor(d$query_cluster, levels = x_ord)
  d$ref_label     <- factor(d$ref_label, levels = rev(y_ord))
  b <- bar_df; b$query_cluster <- factor(b$query_cluster, levels = x_ord)

  xbreaks <- x_ord[suppressWarnings(as.integer(x_ord)) %% TICK_EVERY == 0]
  if (!length(xbreaks)) xbreaks <- x_ord[seq(1, length(x_ord), by = TICK_EVERY)]

  pA <- ggplot(d, aes(query_cluster, ref_label, fill = frac)) +
    geom_tile() +
    # best matches are shown as a red frame around the tile: full red for
    # reciprocal best, pale red for one-directional (forward or reverse) best
    geom_tile(data = subset(d, !is.na(mark)),
              aes(colour = mark), fill = NA, linewidth = 0.3,
              width = 0.88, height = 0.88) +
    scale_fill_gradientn(
      colours = c("white", "#d9e7f5", "#6baed6", "#2171b5", "#08306b"),
      limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1.0"),
      name = "Fraction of\nquery cluster") +
    scale_colour_manual(values = c("reciprocal best" = "#d62728",
                                   "one-directional best" = "#f4a6a6"),
                        name = NULL,
                        drop = FALSE) +
    # a tick at every cluster; numeric labels only every TICK_EVERY so they
    # do not collide at the 1.28 mm tile pitch
    scale_x_discrete(breaks = x_ord,
                     labels = ifelse(x_ord %in% xbreaks, x_ord, ""),
                     expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    # fix the stacking order of the two legend blocks so every panel matches
    guides(colour = guide_legend(order = 1,
                                 override.aes = list(fill = "#6baed6",
                                                     linewidth = 0.4)),
           fill = guide_colourbar(order = 2)) +
    labs(x = NULL, y = y_title) +
    base_theme +
    theme(axis.text.x = element_text(size = BASE_PT),
          axis.text.y = element_text(size = YLAB_PT, margin = margin(r = 0.4, unit = "mm")),
          axis.title.y = element_text(size = BASE_PT),
          panel.border = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.25))

  pTop <- ggplot(b, aes(query_cluster, .data[[conf_var]])) +
    geom_col(fill = "grey70", width = 0.85) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1.0"), expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = conf_lab) +
    base_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = BASE_PT),
          axis.title.y = element_text(size = BASE_PT),
          panel.border = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.25))

  pTop / pA + plot_layout(heights = c(BAR_FRAC, 1), guides = "collect")
}

# ---------------------------------------------------------------------------
# Column order for a FIXED row order.
#
# Walks the rows in the given (native, numeric) order and places each row's
# best-matching column next in the sequence if it has not been placed yet.
# The best match of cluster 0 therefore sits leftmost, that of cluster 1 next,
# and so on, producing a staircase diagonal while the y axis keeps the
# manuscript's cluster numbering. Columns that are never a best match are
# appended at the end, ordered by the row they match best.
# ---------------------------------------------------------------------------
cols_for_fixed_rows <- function(M, row_order) {
  placed <- character(0)
  for (r in row_order) {
    bc <- colnames(M)[which.max(M[r, ])]
    if (!(bc %in% placed)) placed <- c(placed, bc)
  }
  rest <- setdiff(colnames(M), placed)
  if (length(rest)) {
    br <- rownames(M)[apply(M[, rest, drop = FALSE], 2, which.max)]
    rest <- rest[order(match(br, row_order))]
  }
  c(placed, rest)
}

save_fig <- function(p, stem) {
  # drawn at final printed size: tile pitch is TILE_MM, so the point sizes set
  # above are the point sizes that appear on the page
  w <- TILE_MM * nrow(M_fwd) + MARGIN_W_MM
  h <- TILE_MM * ncol(M_fwd) * (1 + BAR_FRAC) + MARGIN_H_MM
  cat(sprintf("%s: page %.1f x %.1f mm; heatmap requested %.1f x %.1f mm\n",
            stem, w, h, TILE_MM * nrow(M_fwd), TILE_MM * ncol(M_fwd)))
  sv <- function(ext, ...) ggsave(file.path(OUTDIR, sprintf("%s_%s.%s", stem, CLASS, ext)),
                                  p, width = w, height = h, units = "mm",
                                  bg = "white", limitsize = FALSE, ...)
  sv("png", dpi = 600)
  sv("pdf")
  # EPS for journal submission. The PostScript device has no alpha channel,
  # so the figure is built from fully opaque fills and strokes only; cairo_ps
  # is used when available because it embeds fonts more reliably.
  eps_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_ps else "eps"
  sv("eps", device = eps_dev, fallback_resolution = 600)
}

row_native <- rownames(M_fwd)        # already in natural cluster order
col_native <- colnames(M_fwd)

# Column seriation against the NATIVE row order: each column is placed at the
# weighted mean row index of its own mass, so the best matches run top-left to
# bottom-right while the y axis keeps the manuscript's cluster numbering.
# Columns carrying no mass at all are parked at the right-hand edge.
save_fig(make_fig(row_native, cols_for_fixed_rows(M_fwd, row_native),
                  conf_var = "med_score", conf_lab = "median score",
                  y_title  = "Cross-species cluster No."),
         "fig_correspondence_rowNumberColSorted")

message("done -> ", normalizePath(OUTDIR))
