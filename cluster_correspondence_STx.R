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
TH_ASSIGNED_MIN <- 0.50   # below this the cluster is flagged low-assignment

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

write.csv(round(M_fwd, 4),
          file.path(OUTDIR, sprintf("correspondence_matrix_forward_STx_%s.csv", CLASS)))
write.csv(round(M_rev, 4),
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
    key %in% paste(best_rev, rownames(M_rev)) ~ "reciprocal best hit",
  key %in% paste(rownames(M_fwd), best_fwd)   ~ "forward best only",
  key %in% paste(best_rev, rownames(M_rev))   ~ "reverse best only",
  TRUE                                        ~ NA_character_)

bar_df <- data.frame(stx_cluster = summ$stx_cluster,
                     frac_assigned = summ$frac_assigned,
                     frac_best = summ$frac_best,
                     n_cells = summ$n_cells, stringsAsFactors = FALSE)

sourcedata <- rbind(
  data.frame(panel = "A (heatmap)", stx_cluster = df0$stx_cluster,
             atlas_label = df0$atlas_label,
             value_type = "fraction_of_assigned_cells_in_stx_cluster",
             value = df0$frac, mark = df0$mark, stringsAsFactors = FALSE),
  data.frame(panel = "B (bar)", stx_cluster = bar_df$stx_cluster,
             atlas_label = NA_character_, value_type = "fraction_assigned_by_RCTD",
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
  d$stx_cluster <- factor(d$stx_cluster, levels = rev(row_ord))
  d$atlas_label <- factor(d$atlas_label, levels = col_ord)
  b <- bar_df; b$stx_cluster <- factor(b$stx_cluster, levels = rev(row_ord))

  pA <- ggplot(d, aes(atlas_label, stx_cluster, fill = frac)) +
    geom_tile(colour = "grey92", linewidth = 0.12) +
    geom_point(data = subset(d, !is.na(mark)),
               aes(shape = mark, colour = mark), size = 1.7, stroke = 0.7) +
    scale_fill_gradientn(
      colours = c("white", "#d9e7f5", "#6baed6", "#2171b5", "#08306b"),
      limits = c(0, 1), name = "fraction of\nassigned cells") +
    scale_shape_manual(values = c("reciprocal best hit" = 21,
                                  "forward best only" = 4,
                                  "reverse best only" = 3), name = NULL) +
    scale_colour_manual(values = c("reciprocal best hit" = "#d62728",
                                   "forward best only" = "#d62728",
                                   "reverse best only" = "grey25"), name = NULL) +
    labs(x = "cross-species atlas subtype, RCTD assignment (Figs. 1-3)",
         y = "Xenium spatial cluster (Figs. 7-9)",
         title = sprintf("Xenium RCTD assignment vs cross-species atlas, %s neurons",
                         CLASS_LONG),
         subtitle = subtitle) +
    base_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))

  bar <- function(var, xint, xlab) {
    ggplot(b, aes(.data[[var]], stx_cluster)) +
      geom_col(fill = "grey70", width = 0.75) +
      geom_vline(xintercept = xint, linetype = "dashed",
                 colour = "#d62728", linewidth = 0.4) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                         labels = c("0", "0.5", "1")) +
      labs(x = xlab, y = NULL) + base_theme +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  pA +
    bar("frac_assigned", TH_ASSIGNED_MIN, "frac. assigned\nby RCTD") +
    bar("frac_best",     TH_FRAC_PRIMARY, "best-match frac.") +
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
  f <- function(ext) file.path(OUTDIR, sprintf("%s_%s.%s", stem, CLASS, ext))
  ggsave(f("png"), p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(f("pdf"), p, width = w, height = h, bg = "white", limitsize = FALSE)
  eps_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_ps else "eps"
  ggsave(f("eps"), p, device = eps_dev, width = w, height = h, bg = "white",
         fallback_resolution = 600, limitsize = FALSE)
}

w_mean <- apply(M_fwd, 2, function(x) {
  s <- sum(x); if (s <= 0) return(Inf); sum(x * seq_len(nrow(M_fwd))) / s })
col_sorted <- names(sort(w_mean))

dg <- diagonal_order(M_fwd)
save_fig(make_fig(dg$rows, dg$cols,
  "rows and columns seriated by greedy maximum matching (strongest match on the diagonal)"),
  "fig_correspondence_STx")

save_fig(make_fig(clu_levels, lab_levels,
  "rows and columns in native cluster / subtype numbering"),
  "fig_correspondence_STx_byClusterNumber")

save_fig(make_fig(clu_levels, cols_for_fixed_rows(M_fwd, clu_levels),
  "rows in native Xenium cluster order; columns sorted to the diagonal"),
  "fig_correspondence_STx_rowNumberColSorted")

message("done -> ", normalizePath(OUTDIR))
