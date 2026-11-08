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
  key %in% fwd_hit & key %in% rev_hit ~ "reciprocal best hit",
  key %in% fwd_hit                    ~ "forward best only",
  key %in% rev_hit                    ~ "reverse best only",
  TRUE                                ~ NA_character_)

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
base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        legend.key.height = unit(0.9, "lines"),
        legend.background = element_blank(),
        plot.title = element_text(size = 12, face = "plain"),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        axis.title.x = element_text(size = 9))

make_fig <- function(row_ord, col_ord, subtitle) {
  d <- df0
  d$query_cluster <- factor(d$query_cluster, levels = rev(row_ord))
  d$ref_label     <- factor(d$ref_label, levels = col_ord)
  b <- bar_df; b$query_cluster <- factor(b$query_cluster, levels = rev(row_ord))

  pA <- ggplot(d, aes(ref_label, query_cluster, fill = frac)) +
    geom_tile(colour = "grey92", linewidth = 0.12) +
    geom_point(data = subset(d, !is.na(mark)),
               aes(shape = mark, colour = mark), size = 1.7, stroke = 0.7) +
    scale_fill_gradientn(
      colours = c("white", "#d9e7f5", "#6baed6", "#2171b5", "#08306b"),
      limits = c(0, 1), name = "fraction of\nquery cluster") +
    scale_shape_manual(values = c("reciprocal best hit" = 21,
                                  "forward best only" = 4,
                                  "reverse best only" = 3), name = NULL) +
    scale_colour_manual(values = c("reciprocal best hit" = "#d62728",
                                   "forward best only" = "#d62728",
                                   "reverse best only" = "grey25"), name = NULL) +
    labs(x = "cross-species atlas subtype (celltype_display, Figs. 1-3)",
         y = "cross-condition cluster (Figs. 4-6)",
         title = sprintf("Label transfer, %s neurons", CLASS_LONG),
         subtitle = subtitle) +
    base_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))

  bar <- function(var, xint, xlab) {
    ggplot(b, aes(.data[[var]], query_cluster)) +
      geom_col(fill = "grey70", width = 0.75) +
      geom_vline(xintercept = xint, linetype = "dashed",
                 colour = "#d62728", linewidth = 0.4) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                         labels = c("0", "0.5", "1")) +
      labs(x = xlab, y = NULL) + base_theme +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  pA +
    bar("med_score", TH_SCORE_MIN,    "median score") +
    bar("frac_best", TH_FRAC_PRIMARY, "best-match frac.") +
    plot_layout(widths = c(6, 1.3, 1.3), guides = "collect")
}

# ---------------------------------------------------------------------------
# Seriation by greedy maximum matching.
#
# Repeatedly takes the largest remaining entry of the matrix and assigns its
# row and column to the next diagonal slot, removing both from further
# consideration. The strongest correspondence of each cluster therefore lands
# exactly on the diagonal, which a weighted-mean ordering cannot guarantee
# when many clusters are unresolved. Rows/columns left over once one axis is
# exhausted are inserted next to whichever partner they match best.
# ---------------------------------------------------------------------------
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

diagonal_order <- function(M) {
  rem_r <- rownames(M); rem_c <- colnames(M)
  pair_r <- character(0); pair_c <- character(0)
  while (length(rem_r) && length(rem_c)) {
    sub <- M[rem_r, rem_c, drop = FALSE]
    if (max(sub) <= 0) break
    ij <- which.max(sub)
    i <- ((ij - 1) %% nrow(sub)) + 1
    j <- ((ij - 1) %/% nrow(sub)) + 1
    pair_r <- c(pair_r, rem_r[i]); pair_c <- c(pair_c, rem_c[j])
    rem_r <- rem_r[-i]; rem_c <- rem_c[-j]
  }
  rkey <- stats::setNames(rep(NA_real_, nrow(M)), rownames(M))
  ckey <- stats::setNames(rep(NA_real_, ncol(M)), colnames(M))
  rkey[pair_r] <- seq_along(pair_r)
  ckey[pair_c] <- seq_along(pair_c)
  lo_r <- names(rkey)[is.na(rkey)]
  if (length(lo_r)) {
    bc <- colnames(M)[apply(M[lo_r, , drop = FALSE], 1, which.max)]
    k  <- ckey[bc]; k[is.na(k)] <- length(pair_r) + 1
    rkey[lo_r] <- k + 0.5
  }
  lo_c <- names(ckey)[is.na(ckey)]
  if (length(lo_c)) {
    br <- rownames(M)[apply(M[, lo_c, drop = FALSE], 2, which.max)]
    k  <- rkey[br]; k[is.na(k)] <- length(pair_c) + 1
    ckey[lo_c] <- k + 0.5
  }
  list(rows = names(sort(rkey)), cols = names(sort(ckey)))
}

save_fig <- function(p, stem) {
  h <- max(6, 0.17 * nrow(M_fwd) + 2.2)
  w <- max(8, 0.17 * ncol(M_fwd) + 5.6)
  ggsave(file.path(OUTDIR, sprintf("%s_%s.png", stem, CLASS)), p,
         width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(OUTDIR, sprintf("%s_%s.pdf", stem, CLASS)), p,
         width = w, height = h, bg = "white", limitsize = FALSE)
  # EPS for journal submission. The PostScript device has no alpha channel,
  # so the figure is built from fully opaque fills and strokes only; cairo_ps
  # is used when available because it embeds fonts more reliably.
  eps_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_ps else "eps"
  ggsave(file.path(OUTDIR, sprintf("%s_%s.eps", stem, CLASS)), p,
         device = eps_dev, width = w, height = h, bg = "white",
         fallback_resolution = 600, limitsize = FALSE)
}

row_native <- rownames(M_fwd)        # already in natural cluster order
col_native <- colnames(M_fwd)

# Column seriation against the NATIVE row order: each column is placed at the
# weighted mean row index of its own mass, so the best matches run top-left to
# bottom-right while the y axis keeps the manuscript's cluster numbering.
# Columns carrying no mass at all are parked at the right-hand edge.
w_mean <- apply(M_fwd, 2, function(x) {
  s <- sum(x); if (s <= 0) return(Inf)
  sum(x * seq_len(nrow(M_fwd))) / s
})
col_sorted <- names(sort(w_mean))

dg <- diagonal_order(M_fwd)
save_fig(make_fig(dg$rows, dg$cols,
  "rows and columns seriated by greedy maximum matching (strongest match on the diagonal)"),
  "fig_correspondence")

save_fig(make_fig(row_native, col_native,
  "rows and columns in native cluster / subtype numbering"),
  "fig_correspondence_byClusterNumber")

save_fig(make_fig(row_native, cols_for_fixed_rows(M_fwd, row_native),
  "rows in native cross-condition cluster order; columns sorted to the diagonal"),
  "fig_correspondence_rowNumberColSorted")

message("done -> ", normalizePath(OUTDIR))
