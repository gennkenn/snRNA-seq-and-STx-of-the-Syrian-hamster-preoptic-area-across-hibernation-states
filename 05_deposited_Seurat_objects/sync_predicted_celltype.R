# =============================================================================
# Make predicted.celltype identical to the RCTD first_type column in the two
# deposited Xenium objects.
#
# The deposited objects carry the complete RCTD output (spot_class, first_type,
# second_type, min_score, singlet_score, conv_all, conv_doublet). Before this
# step, predicted.celltype held the first_type call of the original RCTD run
# and the other columns came from a verification re-run with identical
# settings; the two agreed for >99.5% of cells. So that the deposited columns
# are mutually consistent, predicted.celltype is set equal to first_type here.
# Only the deposit copies are modified; the source objects are untouched.
#
# Usage:  OUTDIR=~/claude/POA_deposit/objects Rscript sync_predicted_celltype.R
# =============================================================================
suppressPackageStartupMessages({ library(Seurat) })
OUTDIR <- path.expand(Sys.getenv("OUTDIR", "~/claude/POA_deposit/objects"))

for (f in c("GSE345359_Xenium_EX_deposit.rds", "GSE345359_Xenium_IN_deposit.rds")) {
  p <- file.path(OUTDIR, f)
  cat("\n==== ", f, " ====\n", sep = "")
  obj <- readRDS(p)
  n0  <- ncol(obj)
  old <- as.character(obj$predicted.celltype)
  new <- as.character(obj$first_type)

  # both columns must label exactly the same set of cells
  stopifnot(identical(is.na(old), is.na(new)))
  n_lab   <- sum(!is.na(new))
  n_diff  <- sum(old != new, na.rm = TRUE)
  cat(sprintf("cells: %d | labelled: %d | predicted.celltype != first_type: %d (%.3f%%)\n",
              n0, n_lab, n_diff, 100 * n_diff / n_lab))

  obj$predicted.celltype <- obj$first_type
  stopifnot(ncol(obj) == n0,
            identical(as.character(obj$predicted.celltype), as.character(obj$first_type)))
  saveRDS(obj, p)
  cat("rewritten: ", p, "\n", sep = "")
  cat("predicted.celltype table (top 5):\n")
  print(head(sort(table(obj$predicted.celltype), decreasing = TRUE), 5))
}
cat("\nDONE\n")
