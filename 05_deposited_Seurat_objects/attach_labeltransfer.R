# =============================================================================
# Attach the per-nucleus label-transfer result to the deposited cross-condition
# (GSE346081) Seurat objects.
#
# Adds two metadata columns from percell_predictions_query_<CLASS>.csv, the
# output of cluster_correspondence_labeltransfer.R (REF_MODE=all):
#   crossspecies_subtype_transferred  cross-species subtype (e1-e44 / i1-i48)
#                                     transferred to the nucleus
#   crossspecies_transfer_score       Seurat prediction.score.max (0-1)
# Only the deposit copies are modified; the source objects are untouched.
#
# Usage:
#   OUTDIR=~/claude/POA_deposit/objects Rscript attach_labeltransfer.R
#   (the two CSVs are expected in the working directory, possibly with a
#    storage prefix in the file name)
# =============================================================================
suppressPackageStartupMessages({ library(Seurat) })
OUTDIR <- path.expand(Sys.getenv("OUTDIR", "~/claude/POA_deposit/objects"))

find_csv <- function(cls) {
  hits <- Sys.glob(sprintf("*percell_predictions_query_%s.csv", cls))
  if (!length(hits)) stop("no per-cell prediction table for ", cls)
  hits[1]
}

JOBS <- list(
  list(cls = "Ex",  obj = "GSE346081_sumwin_vGluT2_deposit.rds"),
  list(cls = "Inh", obj = "GSE346081_sumwin_vGat_deposit.rds"))

for (J in JOBS) {
  p <- file.path(OUTDIR, J$obj)
  cat("\n==== ", J$obj, " ====\n", sep = "")
  obj <- readRDS(p); n0 <- ncol(obj)
  pr  <- read.csv(find_csv(J$cls), stringsAsFactors = FALSE)
  rownames(pr) <- pr$cell
  matched <- sum(colnames(obj) %in% rownames(pr))
  cat(sprintf("cells: %d | prediction rows: %d | matched: %d\n", n0, nrow(pr), matched))
  stopifnot(matched == n0)                     # every nucleus must have a prediction
  # the query cluster recorded by the transfer must equal the deposited clusters
  stopifnot(identical(as.character(pr[colnames(obj), "query_cluster"]),
                      as.character(obj$seurat_clusters)))

  obj$crossspecies_subtype_transferred <- pr[colnames(obj), "predicted_label"]
  obj$crossspecies_transfer_score      <- pr[colnames(obj), "prediction_score"]
  stopifnot(ncol(obj) == n0, !anyNA(obj$crossspecies_subtype_transferred))

  saveRDS(obj, p)
  cat("rewritten: ", p, "\n", sep = "")
  cat("columns: ", paste(colnames(obj@meta.data), collapse = ", "), "\n", sep = "")
  cat("transferred subtypes (top 5):\n")
  print(head(sort(table(obj$crossspecies_subtype_transferred), decreasing = TRUE), 5))
  cat(sprintf("transfer score: median %.3f, IQR %.3f-%.3f\n",
              median(obj$crossspecies_transfer_score),
              quantile(obj$crossspecies_transfer_score, 0.25),
              quantile(obj$crossspecies_transfer_score, 0.75)))
}
cat("\nDONE\n")
