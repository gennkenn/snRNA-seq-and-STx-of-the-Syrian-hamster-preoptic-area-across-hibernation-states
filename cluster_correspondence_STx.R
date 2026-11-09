#!/usr/bin/env Rscript
# =============================================================================
# Correspondence between the Xenium spatial clusters (Figs. 7-9) and the
# cross-species atlas subtypes (Figs. 1-3), Syrian hamster preoptic area.
#
# Unlike the snRNA-seq pairings, no new label transfer is needed here: the
# Xenium objects already carry an RCTD assignment (`predicted.celltype`,
# e1-e44 / i1-i48) alongside their own clustering (`seurat_clusters`), so the
# contingency table of those two columns is the sufficient statistic for the
# whole analysis. This script therefore works from the small counts CSVs
# produced by extract_stx_contingency.R, not from the multi-GB objects.
#
# Two normalisations of the same table are reported:
#   forward  (row-normalised)    what fraction of an STx cluster carries a
#                                given atlas label  -> identity of the cluster
#   reverse  (column-normalised) what fraction of an atlas label's cells sits
#                                in a given STx cluster -> where the subtype went
# A pair that is the best match in BOTH directions is a reciprocal best hit.
#
# Cells RCTD could not assign (label NA) are excluded from both normalisations
# but reported per cluster as `frac_unassigned` — with a 100-gene panel this is
# the honest confidence measure, and it replaces the median prediction score
# used for the snRNA-seq pairings (RCTD's spot_class was not retained here).
#
# Inputs  : counts_<CLASS>.csv  (seurat_clusters x predicted.celltype counts)
# Outputs : correspondence_table_STx_<CLASS>.csv
#           correspondence_matrix_forward_STx_<CLASS>.csv
#           correspondence_matrix_reverse_STx_<CLASS>.csv
#           sourcedata_fig_correspondence_STx_<CLASS>.csv
#           fig_correspondence_STx[_byClusterNumber|_rowNumberColSorted]_<CLASS>
#               .png / .pdf / .eps
#
# Usage: INDIR=. OUTDIR=stx_out CLASS=EX Rscript cluster_correspondence_STx.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})

INDIR  <- Sys.getenv("INDIR",  ".")
OUTDIR <- Sys.getenv("OUTDIR", "stx_out")
CLASS  <- Sys.getenv("CLASS",  "EX")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

CLASS_LONG <- c(EX = "excitatory", IN = "inhibitory")[CLASS]
if (is.na(CLASS_LONG)) CLASS_LONG <- CLASS

# Calling thresholds, identical to the snRNA-seq pairings so that the three
# correspondence tables can be concatenated and read with one set of rules.
TH_FRAC_PRIMARY <- 0.50   # >= 50% of assigned cells on one label
TH_FRAC_SECOND  <- 0.20   # a label holding >= 20% counts as an extra partner
TH_FRAC_MIN     <- 0.30   # below this the cluster is unresolved
TH_ASSIGNED_MIN <- 0.30   # below this the cluster is flagged low-confidence

natural_levels <- function(x) {
  u <- unique(as.character(x)); u <- u[!is.na(u)]
  m <- regmatches(u, regexpr("[0-9]+$", u))
  if (length(m) == length(u) && all(nzchar(m))) {
    u[order(sub("[0-9]+$", "", u), as.integer(m))]
  } else sort(u)
}

# ---------------------------------------------------------------------------
# 1. Load the contingency table and split off the unassigned column
# ---------------------------------------------------------------------------
C <- as.matrix(read.csv(file.path(INDIR, sprintf("counts_%s.csv", CLASS)),
                        row.names = 1, check.names = FALSE))
storage.mode(C) <- "numeric"

na_col <- which(colnames(C) %in% c("NA", "NA.", "") | is.na(colnames(C)))
n_unassigned <- if (length(na_col)) C[, na_col, drop = FALSE] |> rowSums() else
                 rep(0, nrow(C))
if (length(na_col)) C <- C[, -na_col, drop = FALSE]

n_total    <- rowSums(C) + n_unassigned
n_assigned <- rowSums(C)

clu_levels <- natural_levels(rownames(C))
lab_levels <- natural_levels(colnames(C))
C <- C[clu_levels, lab_levels, drop = FALSE]
n_total    <- n_total[clu_levels]
n_assigned <- n_assigned[clu_levels]
n_unassigned <- n_unassigned[clu_levels]

cat(sprintf("%s: %d STx clusters x %d atlas labels | %s cells (%.1f%% unassigned)\n",
            CLASS, nrow(C), ncol(C), format(sum(n_total), big.mark = ","),
            100 * sum(n_unassigned) / sum(n_total)))

# ---------------------------------------------------------------------------
# 2. The two normalisations
# ---------------------------------------------------------------------------
M_fwd <- sweep(C, 1, pmax(n_assigned, 1), "/")            # cluster -> label
M_rev <- t(sweep(C, 2, pmax(colSums(C), 1), "/"))          # label -> cluster

write.csv(round(M_fwd, 6),
          file.path(OUTDIR, sprintf("correspondence_matrix_forward_STx_%s.csv", CLASS)))
write.csv(round(M_rev, 6),
          file.path(OUTDIR, sprintf("correspondence_matrix_reverse_STx_%s.csv", CLASS)))

best_fwd <- colnames(M_fwd)[apply(M_fwd, 1, which.max)]; names(best_fwd) <- rownames(M_fwd)
best_rev <- colnames(M_rev)[apply(M_rev, 1, which.max)]; names(best_rev) <- rownames(M_rev)

# ---------------------------------------------------------------------------
# 3. Per-cluster summary table
# ---------------------------------------------------------------------------
summ <- lapply(rownames(M_fwd), function(q) {
  v  <- sort(M_fwd[q, ], decreasing = TRUE)
  bl <- names(v)[1]
  partners <- names(v)[v >= TH_FRAC_SECOND]
  data.frame(
    comparison          = sprintf("STx(Figs7-9) vs cross-species(Figs1-3), %s", CLASS),
    stx_cluster         = q,
    n_cells             = unname(n_total[q]),
    n_assigned          = unname(n_assigned[q]),
    frac_unassigned     = unname(n_unassigned[q] / pmax(n_total[q], 1)),
    best_label          = bl,
    frac_best           = unname(v[1]),
    second_label        = names(v)[2],
    frac_second         = unname(v[2]),
    n_partners_ge20pct  = length(partners),
    partners_ge20pct    = paste(partners, collapse = ";"),
    reverse_best_of_best_label = unname(best_rev[bl]),
    reciprocal          = !is.na(best_rev[bl]) && best_rev[bl] == q,
    many_to_one         = sum(best_fwd == bl) > 1,
    stringsAsFactors = FALSE)
}) |> dplyr::bind_rows()

summ$frac_assigned <- 1 - summ$frac_unassigned
summ$relation <- with(summ, dplyr::case_when(
  frac_best < TH_FRAC_MIN                                   ~ "unresolved",
  reciprocal & frac_best >= TH_FRAC_PRIMARY &
    n_partners_ge20pct == 1                                 ~ "1:1",
  n_partners_ge20pct >= 2                                   ~ "1:many",
  many_to_one                                               ~ "many:1",
  TRUE                                                      ~ "1:1 (weak)"))
summ$low_assignment_flag <- summ$frac_assigned < TH_ASSIGNED_MIN

summ <- summ[order(match(summ$stx_cluster, clu_levels)), ]
write.csv(summ, file.path(OUTDIR, sprintf("correspondence_table_STx_%s.csv", CLASS)),
          row.names = FALSE)

cat("\n---- relation summary ----\n"); print(table(summ$relation))
cat("reciprocal best hits: ", sum(summ$reciprocal), " / ", nrow(summ), "\n", sep = "")
cat("clusters flagged low-assignment (<50% assigned): ",
    sum(summ$low_assignment_flag), "\n", sep = "")
cat("atlas labels never a best match: ",
    paste(setdiff(lab_levels, summ$best_label), collapse = ", "), "\n", sep = "")
cat("median frac_best: ", round(median(summ$frac_best), 3),
    " | median frac_assigned: ", round(median(summ$frac_assigned), 3), "\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Figure + source data
# ---------------------------------------------------------------------------
df0 <- as.data.frame(as.table(M_fwd))
colnames(df0) <- c("stx_cluster", "atlas_label", "frac")
df0$stx_cluster <- as.character(df0$stx_cluster)
df0$atlas_label <- as.character(df0$atlas_label)

key <- paste(df0$stx_cluster, df0$atlas_label)
df0$mark <- dplyr::case_when(
  key %in% paste(rownames(M_fwd), best_fwd) &
    key %in% paste(best_rev, rownames(M_rev)) ~ "reciprocal best",
  key %in% paste(rownames(M_fwd), best_fwd) ~ "one-directional best",
  key %in% paste(best_rev, rownames(M_rev)) ~ "one-directional best",
  TRUE                                        ~ NA_character_)
df0$mark <- factor(df0$mark, levels = c("reciprocal best", "one-directional best"))

bar_df <- data.frame(stx_cluster = summ$stx_cluster,
                     frac_assigned = summ$frac_assigned,
                     frac_best = summ$frac_best,
                     n_cells = summ$n_cells, stringsAsFactors = FALSE)

sourcedata <- rbind(
  data.frame(panel = "A (heatmap)", stx_cluster = df0$stx_cluster,
             atlas_label = df0$atlas_label,
             value_type = "fraction_of_singlet_cells_in_stx_cluster",
             value = df0$frac, mark = df0$mark, stringsAsFactors = FALSE),
  data.frame(panel = "B (bar)", stx_cluster = bar_df$stx_cluster,
             atlas_label = NA_character_, value_type = "fraction_singlet_by_RCTD",
             value = bar_df$frac_assigned, mark = NA_character_),
  data.frame(panel = "C (bar)", stx_cluster = bar_df$stx_cluster,
             atlas_label = NA_character_, value_type = "fraction_in_best_matching_label",
             value = bar_df$frac_best, mark = NA_character_),
  data.frame(panel = "n", stx_cluster = bar_df$stx_cluster,
             atlas_label = NA_character_, value_type = "n_cells_in_stx_cluster",
             value = as.numeric(bar_df$n_cells), mark = NA_character_))
write.csv(sourcedata,
          file.path(OUTDIR, sprintf("sourcedata_fig_correspondence_STx_%s.csv", CLASS)),
          row.names = FALSE)

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
.cal <- list(EX = c(31.9, 11.6), IN = c(32.9, 11.4))[[CLASS]]
MARGIN_W_MM <- as.numeric(Sys.getenv("MARGIN_W_MM", .cal[1]))
MARGIN_H_MM <- as.numeric(Sys.getenv("MARGIN_H_MM", .cal[2]))
BAR_FRAC    <- as.numeric(Sys.getenv("BAR_FRAC", "0.22"))  # bar height / heatmap

base_theme <- theme_bw(base_size = BASE_PT) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length = unit(0.6, "mm"),
        legend.background = element_blank(),
        legend.key.height = unit(3.2, "mm"),
        legend.key.width  = unit(3.2, "mm"),
        legend.spacing.y  = unit(1, "mm"),
        legend.text  = element_text(size = BASE_PT),
        legend.title = element_text(size = BASE_PT),
        legend.margin = margin(0, 0, 0, 1, "mm"),
        plot.margin = margin(1, 1, 1, 1, "mm"))

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
  d$stx_cluster <- factor(d$stx_cluster, levels = x_ord)
  d$atlas_label     <- factor(d$atlas_label, levels = rev(y_ord))
  b <- bar_df; b$stx_cluster <- factor(b$stx_cluster, levels = x_ord)

  xbreaks <- x_ord[suppressWarnings(as.integer(x_ord)) %% TICK_EVERY == 0]
  if (!length(xbreaks)) xbreaks <- x_ord[seq(1, length(x_ord), by = TICK_EVERY)]

  pA <- ggplot(d, aes(stx_cluster, atlas_label, fill = frac)) +
    geom_tile() +
    # best matches are shown as a red frame around the tile: full red for
    # reciprocal best, pale red for one-directional (forward or reverse) best
    geom_tile(data = subset(d, !is.na(mark)),
              aes(colour = mark), fill = NA, linewidth = 0.3,
              width = 0.88, height = 0.88) +
    scale_fill_gradientn(
      colours = c("white", "#d9e7f5", "#6baed6", "#2171b5", "#08306b"),
      limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1.0"),
      name = "Fraction of\nsinglet cells") +
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

  pTop <- ggplot(b, aes(stx_cluster, .data[[conf_var]])) +
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

save_fig(make_fig(clu_levels, cols_for_fixed_rows(M_fwd, clu_levels),
                  conf_var = "frac_assigned", conf_lab = "frac. singlet",
                  y_title  = "Cross-species cluster No."),
         "fig_correspondence_STx_rowNumberColSorted")

message("done -> ", normalizePath(OUTDIR))
