#!/usr/bin/env Rscript
# =============================================================================
# Cluster correspondence by bidirectional label transfer
# Syrian hamster preoptic area (POA) atlas -- excitatory or inhibitory neurons
#
# Establishes the correspondence between:
#   REF  : cross-species snRNA-seq dataset (Figs. 1-3)
#          object MPN_merged.neuron_vGluT2_QC.rds   (excitatory)
#                 MPN_merged.neuron_vGat_QC.rds     (inhibitory)
#          labels = meta.data$celltype_display
#          NOTE: contains BOTH mouse and hamster cells, in a 1:1-ortholog
#                gene space (mouse gene symbols).
#   QRY  : cross-condition (hibernating vs non-hibernating) snRNA-seq
#          dataset (Figs. 4-6)
#          object sumwinHam_merged.neuron_vGluT2.rds        (excitatory)
#                 sumwinHam_merged.neuron_vGat_afterQC.rds  (inhibitory)
#          labels = Idents()  (Seurat cluster ids)
#
# The same script serves both neuronal classes; set CLASS=Ex or CLASS=Inh.
#
# The two datasets share NO cells (independent animal cohorts), so the
# correspondence is established computationally.
#
# Method
#   1. Restrict both objects to the shared gene space.
#   2. Re-normalise both from raw counts with a single, common pipeline
#      (LogNormalize) so that reference and query are directly comparable.
#      The query's SCT assay is deliberately NOT used, because the reference
#      was log-normalised; mixing normalisations biases anchor finding.
#   3. Forward transfer : REF -> QRY   (each query cell gets a predicted e-label)
#   4. Reverse transfer : QRY -> REF   (each reference hamster cell gets a
#      predicted cross-condition cluster id)
#   5. Cross-tabulate, take reciprocal best hits, and classify each pair as
#      1:1 / 1:many / many:1 / unresolved using explicit thresholds.
#
#   The primary reference is the HAMSTER-ONLY subset of REF (species-matched
#   to the query). A sensitivity run using all REF cells (mouse + hamster) is
#   also performed, because several e-labels contain few or no hamster cells
#   and are therefore not representable in a hamster-only reference.
#
# Outputs (written to OUTDIR; <CLASS> is Ex or Inh)
#   correspondence_table_<CLASS>.csv              per-query-cluster summary
#   correspondence_matrix_forward_<CLASS>.csv     query cluster x ref label
#   correspondence_matrix_reverse_<CLASS>.csv     ref label x query cluster
#   percell_predictions_query_<CLASS>.csv         per-cell predicted label + score
#   percell_predictions_reference_<CLASS>.csv     per-cell predicted cluster + score
#   reference_label_representation_<CLASS>.csv    hamster cell count per ref label
#
# The figure and its source-data file are produced separately by
# replot_correspondence_from_csv.R, which reads the CSVs above.
#
# Seurat 5.5.1 / SeuratObject 5.4.0
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

set.seed(1234)

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------
REF_PATH <- Sys.getenv("REF_RDS", "MPN_merged.neuron_vGluT2_QC.rds")
QRY_PATH <- Sys.getenv("QRY_RDS", "sumwinHam_merged.neuron_vGluT2.rds")
OUTDIR   <- Sys.getenv("OUTDIR",  "correspondence_out")
CLASS    <- Sys.getenv("CLASS",   "Ex")   # "Ex" (vGluT2+) or "Inh" (vGat+)

CLASS_LONG <- c(Ex = "excitatory (vGluT2+)",
                Inh = "inhibitory (vGat+)")[CLASS]
if (is.na(CLASS_LONG)) CLASS_LONG <- CLASS
out <- function(stem, ext) file.path(OUTDIR, sprintf("%s_%s.%s", stem, CLASS, ext))

REF_LABEL_COL <- "celltype_display"   # e1 - e44
SPECIES_COL   <- "Species"
REF_SPECIES   <- "hamster"            # species-matched reference subset

N_HVG   <- 2000   # variable features used to build each reference
N_PCS   <- 30     # dimensions used for anchor finding and weighting
K_ANCHOR <- 5
MIN_REF_CELLS <- 10   # e-labels with fewer hamster cells are not assessable

# Correspondence-calling thresholds (stated explicitly so they can be reported)
TH_FRAC_PRIMARY  <- 0.50  # >= 50% of a query cluster assigned to one label
TH_FRAC_SECOND   <- 0.20  # a label holding >= 20% counts as an extra partner
TH_FRAC_MIN      <- 0.30  # below this, the cluster is called unresolved
TH_SCORE_MIN     <- 0.40  # median prediction.score.max floor

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 1. Helper functions
# ----------------------------------------------------------------------------

#' Collapse split v5 layers and re-normalise from raw counts on a fixed
#' gene universe, then build an independent PCA for anchor projection.
prep_object <- function(obj, features_universe, assay = "RNA",
                        nhvg = N_HVG, npcs = N_PCS) {
  DefaultAssay(obj) <- assay
  # v5 objects may hold one counts layer per sample; anchor finding needs one.
  if (length(Layers(obj[[assay]], search = "counts")) > 1) {
    obj[[assay]] <- JoinLayers(obj[[assay]])
  }
  obj <- subset(obj, features = features_universe)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = nhvg, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs,
                verbose = FALSE)
  obj
}

#' One direction of anchor-based label transfer.
#' Returns the TransferData prediction data.frame (rows = query cells).
run_transfer <- function(reference, query, refdata,
                         npcs = N_PCS, k.anchor = K_ANCHOR) {
  feats <- intersect(VariableFeatures(reference), rownames(query))
  message(sprintf("  anchor features: %d", length(feats)))
  anchors <- FindTransferAnchors(
    reference            = reference,
    query                = query,
    normalization.method = "LogNormalize",
    reduction            = "pcaproject",
    features             = feats,
    dims                 = 1:npcs,
    k.anchor             = k.anchor,
    k.filter             = NA,      # no spatial filtering: references are small
    verbose              = FALSE
  )
  message(sprintf("  anchors found: %d", nrow(anchors@anchors)))
  TransferData(
    anchorset        = anchors,
    refdata          = refdata,
    weight.reduction = "pcaproject",
    dims             = 1:npcs,
    verbose          = FALSE
  )
}

#' Natural ordering of label vectors: "e1, e2, ... e44" and "0, 1, ... 50"
#' rather than the lexicographic "e1, e10, e11" / "0, 1, 10".
natural_levels <- function(x) {
  u <- unique(as.character(x))
  m <- regmatches(u, regexpr("[0-9]+$", u))
  if (length(m) == length(u) && all(nzchar(m))) {
    pre <- sub("[0-9]+$", "", u)
    u[order(pre, as.integer(m))]
  } else {
    sort(u)
  }
}

#' Row-normalised contingency matrix: rows = groups, cols = predicted labels.
#' Columns are forced to `all_labels` so that never-predicted labels are kept.
frac_matrix <- function(group, predicted, all_labels) {
  tb <- table(factor(group), factor(predicted, levels = all_labels))
  sweep(as.matrix(tb), 1, pmax(rowSums(tb), 1), "/")
}

# ----------------------------------------------------------------------------
# 2. Load
# ----------------------------------------------------------------------------
message("== loading objects ==")
ref_full <- readRDS(REF_PATH)
qry_raw  <- readRDS(QRY_PATH)

for (nm in c("reference" = "ref_full", "query" = "qry_raw")) {
  o <- get(nm)
  if (!"RNA" %in% names(o@assays)) {
    stop(sprintf("Object '%s' has no RNA assay. Assays present: %s",
                 nm, paste(names(o@assays), collapse = ", ")))
  }
}

if (!REF_LABEL_COL %in% colnames(ref_full@meta.data)) {
  stop(sprintf("Reference label column '%s' not found. Available: %s",
               REF_LABEL_COL, paste(colnames(ref_full@meta.data), collapse = ", ")))
}
if (!SPECIES_COL %in% colnames(ref_full@meta.data)) {
  stop(sprintf("Species column '%s' not found in the reference.", SPECIES_COL))
}

ref_full$.label <- as.character(ref_full[[REF_LABEL_COL]][, 1])

# Query labels: Idents() by default, or a named metadata column via QRY_LABEL_COL.
QRY_LABEL_COL <- Sys.getenv("QRY_LABEL_COL", "")
if (nzchar(QRY_LABEL_COL)) {
  if (!QRY_LABEL_COL %in% colnames(qry_raw@meta.data)) {
    stop(sprintf("Query label column '%s' not found. Available: %s",
                 QRY_LABEL_COL, paste(colnames(qry_raw@meta.data), collapse = ", ")))
  }
  qry_raw$.label <- as.character(qry_raw[[QRY_LABEL_COL]][, 1])
} else {
  qry_raw$.label <- as.character(Idents(qry_raw))
}

cat(sprintf("reference labels: %d | query labels: %d\n",
            length(unique(ref_full$.label)), length(unique(qry_raw$.label))))

e_levels <- natural_levels(ref_full$.label)   # reference subtype labels
q_levels <- natural_levels(qry_raw$.label)    # query cluster ids

shared_genes <- intersect(rownames(ref_full[["RNA"]]), rownames(qry_raw[["RNA"]]))
message(sprintf("reference: %d cells / query: %d cells / shared genes: %d",
                ncol(ref_full), ncol(qry_raw), length(shared_genes)))

# How well is each e-label represented by hamster cells? This bounds what a
# species-matched reference can possibly recover.
rep_tab <- as.data.frame(table(ref_full$.label,
                               ref_full[[SPECIES_COL]][, 1])) |>
  tidyr::pivot_wider(names_from = Var2, values_from = Freq) |>
  dplyr::rename(label = Var1) |>
  dplyr::mutate(label = factor(label, levels = e_levels)) |>
  dplyr::arrange(label)
rep_tab$assessable_hamster_ref <- rep_tab[[REF_SPECIES]] >= MIN_REF_CELLS
write.csv(rep_tab, out("reference_label_representation", "csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# 3. Prepare the three working objects
# ----------------------------------------------------------------------------
message("== preparing objects ==")
qry <- prep_object(qry_raw, shared_genes)

ref_ham <- subset(ref_full, cells = colnames(ref_full)[
  ref_full[[SPECIES_COL]][, 1] == REF_SPECIES])
keep_lab <- as.character(rep_tab$label[rep_tab$assessable_hamster_ref])
ref_ham <- subset(ref_ham, cells = colnames(ref_ham)[ref_ham$.label %in% keep_lab])
ref_ham <- prep_object(ref_ham, shared_genes)
message(sprintf("hamster-only reference: %d cells, %d labels retained (>= %d cells)",
                ncol(ref_ham), length(keep_lab), MIN_REF_CELLS))

ref_all <- prep_object(ref_full, shared_genes)   # sensitivity reference

# ----------------------------------------------------------------------------
# 4. Forward transfer: reference e-labels -> query cells
# ----------------------------------------------------------------------------
message("== forward transfer (hamster-only reference) ==")
pred_fwd <- run_transfer(ref_ham, qry, refdata = ref_ham$.label)

message("== forward transfer (all-species reference, sensitivity) ==")
pred_fwd_all <- run_transfer(ref_all, qry, refdata = ref_all$.label)

qry$pred_label       <- pred_fwd$predicted.id
qry$pred_score       <- pred_fwd$prediction.score.max
qry$pred_label_allsp <- pred_fwd_all$predicted.id
qry$pred_score_allsp <- pred_fwd_all$prediction.score.max

percell_q <- data.frame(
  cell              = colnames(qry),
  query_cluster     = qry$.label,
  condition         = if ("Mode" %in% colnames(qry@meta.data)) qry$Mode else NA,
  predicted_label   = qry$pred_label,
  prediction_score  = qry$pred_score,
  predicted_label_allspecies_ref = qry$pred_label_allsp,
  prediction_score_allspecies_ref = qry$pred_score_allsp,
  stringsAsFactors  = FALSE
)
write.csv(percell_q, out("percell_predictions_query", "csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# 5. Reverse transfer: query cluster ids -> reference hamster cells
# ----------------------------------------------------------------------------
message("== reverse transfer (query as reference) ==")
ref_ham_q <- ref_ham   # query role in this direction
pred_rev <- run_transfer(qry, ref_ham_q, refdata = qry$.label)

ref_ham$pred_cluster <- pred_rev$predicted.id
ref_ham$pred_score   <- pred_rev$prediction.score.max

percell_r <- data.frame(
  cell             = colnames(ref_ham),
  reference_label  = ref_ham$.label,
  predicted_cluster = ref_ham$pred_cluster,
  prediction_score = ref_ham$pred_score,
  stringsAsFactors = FALSE
)
write.csv(percell_r, out("percell_predictions_reference", "csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# 6. Contingency matrices and reciprocal best hits
# ----------------------------------------------------------------------------
message("== building correspondence ==")
lab_used <- intersect(e_levels, keep_lab)

M_fwd <- frac_matrix(factor(qry$.label, levels = q_levels),
                     qry$pred_label, lab_used)                 # cluster x label
M_rev <- frac_matrix(factor(ref_ham$.label, levels = lab_used),
                     ref_ham$pred_cluster, q_levels)           # label x cluster

write.csv(round(M_fwd, 6), out("correspondence_matrix_forward", "csv"))
write.csv(round(M_rev, 6), out("correspondence_matrix_reverse", "csv"))

best_fwd <- colnames(M_fwd)[apply(M_fwd, 1, which.max)]   # per query cluster
names(best_fwd) <- rownames(M_fwd)
best_rev <- colnames(M_rev)[apply(M_rev, 1, which.max)]   # per e-label
names(best_rev) <- rownames(M_rev)

med_score <- tapply(qry$pred_score, factor(qry$.label, levels = q_levels),
                    median)

summ <- lapply(rownames(M_fwd), function(q) {
  v  <- sort(M_fwd[q, ], decreasing = TRUE)
  bl <- names(v)[1]
  partners <- names(v)[v >= TH_FRAC_SECOND]
  recip <- !is.na(best_rev[bl]) && best_rev[bl] == q
  many_to_one <- sum(best_fwd == bl) > 1
  data.frame(
    query_cluster       = q,
    n_cells             = sum(qry$.label == q),
    best_label          = bl,
    frac_best           = unname(v[1]),
    second_label        = names(v)[2],
    frac_second         = unname(v[2]),
    n_partners_ge20pct  = length(partners),
    partners_ge20pct    = paste(partners, collapse = ";"),
    median_pred_score   = unname(med_score[q]),
    reverse_best_of_best_label = unname(best_rev[bl]),
    reciprocal          = recip,
    many_to_one         = many_to_one,
    best_label_allspecies_ref = names(sort(table(
      qry$pred_label_allsp[qry$.label == q]), decreasing = TRUE))[1],
    stringsAsFactors = FALSE
  )
}) |> dplyr::bind_rows()

summ$relation <- with(summ, dplyr::case_when(
  frac_best < TH_FRAC_MIN | median_pred_score < TH_SCORE_MIN ~ "unresolved",
  reciprocal & frac_best >= TH_FRAC_PRIMARY & n_partners_ge20pct == 1 ~ "1:1",
  n_partners_ge20pct >= 2 ~ "1:many",
  many_to_one ~ "many:1",
  TRUE ~ "1:1 (weak)"
))
summ$agree_with_allspecies_ref <- summ$best_label == summ$best_label_allspecies_ref

summ <- summ[order(match(summ$query_cluster, q_levels)), ]
write.csv(summ, out("correspondence_table", "csv"), row.names = FALSE)

cat("\n---- relation summary ----\n"); print(table(summ$relation))
cat("---- reciprocal best hits:", sum(summ$reciprocal), "/", nrow(summ), "----\n")
cat("---- agreement with all-species reference:",
    sum(summ$agree_with_allspecies_ref), "/", nrow(summ), "----\n")
cat("---- e-labels never used as a best match:",
    paste(setdiff(lab_used, summ$best_label), collapse = ", "), "----\n\n")

# ----------------------------------------------------------------------------
# 7. Figure
#
# The published panel (Fig. 5i, excitatory; Fig. 6f, inhibitory) is drawn by
# replot_correspondence_from_csv.R from the CSVs written above, so that
# changes to the figure never require re-running the label transfer.
# ----------------------------------------------------------------------------
message("tables written -> ", normalizePath(OUTDIR),
        "  (draw the figure with replot_correspondence_from_csv.R)")
