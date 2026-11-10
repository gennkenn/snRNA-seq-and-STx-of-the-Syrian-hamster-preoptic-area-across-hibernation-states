# =============================================================================
# Independent check of the label-transfer correspondence by marker overlap.
#
# For every subtype of the cross-species dataset (REF) and every subtype of the
# cross-condition dataset (QRY), marker genes are identified separately within
# each dataset with FindMarkers (one subtype vs all other subtypes of the same
# dataset; Wilcoxon, only.pos, min.pct 0.1, log2FC > 0.25, Bonferroni-adjusted
# p < 0.05), restricted to the gene space shared by the two datasets. The
# overlap of the two marker sets is then scored for every QRY x REF pair by the
# Jaccard index and a hypergeometric test over the shared gene space. Because
# the marker sets are derived without any reference to the other dataset, the
# overlap is independent of the label transfer and can be used to check it:
# for each QRY subtype we report whether the label-transfer best match is also
# the REF subtype with the highest marker overlap, and its rank.
#
# Inputs (env): REF_RDS, QRY_RDS, OUTDIR, CLASS (Ex|Inh), LT_TABLE
#   (correspondence_table_<CLASS>.csv from cluster_correspondence_labeltransfer.R),
#   optional QRY_LABEL_COL, N_SHOW (pairs listed in detail, default 5),
#   TOP_N (cap on markers per subtype, by avg_log2FC; default 100; 0 = no cap)
# Outputs: marker_overlap_jaccard_<CLASS>.csv      QRY x REF Jaccard matrix
#          marker_overlap_summary_<CLASS>.csv      per QRY subtype
#          marker_overlap_selected_<CLASS>.csv     shared markers of N_SHOW pairs
#          markers_ref_<CLASS>.csv, markers_qry_<CLASS>.csv
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

REF_PATH <- Sys.getenv("REF_RDS"); QRY_PATH <- Sys.getenv("QRY_RDS")
OUTDIR   <- Sys.getenv("OUTDIR", "."); CLASS <- Sys.getenv("CLASS", "Ex")
LT_TABLE <- Sys.getenv("LT_TABLE", "")
QRY_LABEL_COL <- Sys.getenv("QRY_LABEL_COL", "")
N_SHOW <- as.integer(Sys.getenv("N_SHOW", "5"))
TOP_N  <- as.integer(Sys.getenv("TOP_N", "100"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
out <- function(stem) file.path(OUTDIR, sprintf("%s_%s.csv", stem, CLASS))

natural <- function(x) { x <- unique(x); x[order(suppressWarnings(as.numeric(gsub("\\D", "", x))), x)] }

prep <- function(obj, genes) {
  DefaultAssay(obj) <- "RNA"
  if (length(Layers(obj[["RNA"]], search = "counts")) > 1) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- subset(obj, features = genes)
  NormalizeData(obj, verbose = FALSE)
}

ref <- readRDS(REF_PATH); qry <- readRDS(QRY_PATH)
ref$.label <- as.character(ref$celltype_display)
qry$.label <- if (nzchar(QRY_LABEL_COL)) as.character(qry[[QRY_LABEL_COL]][, 1]) else as.character(Idents(qry))
shared <- intersect(rownames(ref[["RNA"]]), rownames(qry[["RNA"]]))
cat(sprintf("shared genes: %d | ref %d cells / %d subtypes | qry %d cells / %d subtypes\n",
            length(shared), ncol(ref), length(unique(ref$.label)), ncol(qry), length(unique(qry$.label))))
ref <- prep(ref, shared); qry <- prep(qry, shared)
Idents(ref) <- ref$.label; Idents(qry) <- qry$.label

find_all <- function(obj, tag) {
  t0 <- Sys.time()
  m <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                      test.use = "wilcox", verbose = FALSE)
  m <- m[m$p_val_adj < 0.05, ]
  cat(sprintf("%s: %d markers across %d subtypes (%.1f min)\n", tag, nrow(m),
              length(unique(m$cluster)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  m
}
m_ref <- find_all(ref, "REF"); m_qry <- find_all(qry, "QRY")
write.csv(m_ref, out("markers_ref"), row.names = FALSE)
write.csv(m_qry, out("markers_qry"), row.names = FALSE)

to_sets <- function(m) {
  m <- m[order(m$cluster, -m$avg_log2FC), ]
  if (TOP_N > 0) m <- m %>% group_by(cluster) %>% slice_head(n = TOP_N) %>% ungroup()
  split(m$gene, as.character(m$cluster))
}
S_ref <- to_sets(m_ref); S_qry <- to_sets(m_qry)
r_lev <- natural(names(S_ref)); q_lev <- natural(names(S_qry))
N <- length(shared)

J <- matrix(0, length(q_lev), length(r_lev), dimnames = list(q_lev, r_lev))
P <- J; OV <- J
for (q in q_lev) for (r in r_lev) {
  a <- S_qry[[q]]; b <- S_ref[[r]]; k <- length(intersect(a, b))
  OV[q, r] <- k
  J[q, r]  <- k / length(union(a, b))
  P[q, r]  <- phyper(k - 1, length(b), N - length(b), length(a), lower.tail = FALSE)
}
write.csv(round(J, 4), out("marker_overlap_jaccard"))
write.csv(OV, out("marker_overlap_count"))

# ---- compare with the label-transfer best match ------------------------------
lt <- read.csv(LT_TABLE, stringsAsFactors = FALSE)
lt$query_cluster <- as.character(lt$query_cluster)
summ <- lapply(q_lev, function(q) {
  jr <- sort(J[q, ], decreasing = TRUE)
  ltb <- lt$best_label[lt$query_cluster == q]
  if (!length(ltb)) ltb <- NA
  data.frame(
    query_cluster = q, n_markers_qry = length(S_qry[[q]]),
    lt_best = ltb, lt_frac_best = lt$frac_best[lt$query_cluster == q][1],
    lt_reciprocal = lt$reciprocal[lt$query_cluster == q][1],
    marker_best = names(jr)[1], marker_best_jaccard = unname(jr[1]),
    marker_second = names(jr)[2], marker_second_jaccard = unname(jr[2]),
    lt_best_jaccard = if (!is.na(ltb) && ltb %in% r_lev) J[q, ltb] else NA,
    lt_best_overlap = if (!is.na(ltb) && ltb %in% r_lev) OV[q, ltb] else NA,
    lt_best_n_markers_ref = if (!is.na(ltb) && ltb %in% names(S_ref)) length(S_ref[[ltb]]) else NA,
    lt_best_hyper_p = if (!is.na(ltb) && ltb %in% r_lev) P[q, ltb] else NA,
    rank_of_lt_best = if (!is.na(ltb) && ltb %in% r_lev) which(names(jr) == ltb) else NA,
    agree = !is.na(ltb) && names(jr)[1] == ltb,
    stringsAsFactors = FALSE)
}) %>% bind_rows()
write.csv(summ, out("marker_overlap_summary"), row.names = FALSE)

cat("\n---- label-transfer best match vs marker-overlap best match ----\n")
cat(sprintf("agree (rank 1): %d / %d | within top 3: %d / %d | median rank: %s\n",
            sum(summ$agree), nrow(summ), sum(summ$rank_of_lt_best <= 3, na.rm = TRUE), nrow(summ),
            median(summ$rank_of_lt_best, na.rm = TRUE)))
cat(sprintf("reciprocal pairs: agree %d / %d\n",
            sum(summ$agree[summ$lt_reciprocal %in% c(TRUE, "TRUE")]),
            sum(summ$lt_reciprocal %in% c(TRUE, "TRUE"))))

# ---- N_SHOW reciprocal pairs spanning the frac_best range, with shared genes --
rec <- summ[summ$lt_reciprocal %in% c(TRUE, "TRUE") & !is.na(summ$lt_frac_best), ]
rec <- rec[order(-rec$lt_frac_best), ]
pick <- rec[unique(round(seq(1, nrow(rec), length.out = N_SHOW))), ]
sel <- lapply(seq_len(nrow(pick)), function(i) {
  q <- pick$query_cluster[i]; r <- pick$lt_best[i]
  g <- intersect(S_qry[[q]], S_ref[[r]])
  g <- g[order(-(m_qry$avg_log2FC[match(g, m_qry$gene[m_qry$cluster == q])]))]
  data.frame(query_cluster = q, ref_subtype = r, lt_frac_best = pick$lt_frac_best[i],
             n_markers_qry = length(S_qry[[q]]), n_markers_ref = length(S_ref[[r]]),
             n_shared = length(g), jaccard = round(J[q, r], 3),
             hypergeom_p = signif(P[q, r], 3), rank_of_lt_best = pick$rank_of_lt_best[i],
             shared_markers = paste(head(g, 40), collapse = " "), stringsAsFactors = FALSE)
}) %>% bind_rows()
write.csv(sel, out("marker_overlap_selected"), row.names = FALSE)
cat("\n---- selected pairs ----\n")
print(sel[, c("query_cluster", "ref_subtype", "lt_frac_best", "n_markers_qry", "n_markers_ref",
              "n_shared", "jaccard", "hypergeom_p", "rank_of_lt_best")], row.names = FALSE)
cat("\nDONE\n")
