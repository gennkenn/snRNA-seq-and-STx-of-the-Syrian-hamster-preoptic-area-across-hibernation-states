# =============================================================================
# Attach the RCTD results_df columns to the excitatory Xenium deposit object.
#
# make_deposit_objects.R attached these columns to the inhibitory object but
# not to the excitatory one, because rctd_results_EX.csv.gz was no longer
# present on the host when that job ran (the excitatory RCTD had been produced
# by an earlier job whose scratch directory was cleaned). This script repairs
# only that one object, using exactly the merge the original script used, and
# leaves the other five deposit objects untouched.
#
# Usage:
#   OUTDIR=~/claude/POA_deposit/objects \
#   RCTD=~/claude/POA_cluster_correspondence/rctd/rctd_results_EX.csv.gz \
#     Rscript attach_rctd_EX.R
# =============================================================================
suppressPackageStartupMessages({ library(Seurat) })

OUTDIR <- path.expand(Sys.getenv("OUTDIR", "~/claude/POA_deposit/objects"))
RCTD   <- path.expand(Sys.getenv("RCTD", ""))
OBJ    <- file.path(OUTDIR, "GSE345359_Xenium_EX_deposit.rds")

RCTD_COLS <- c("spot_class", "first_type", "second_type",
               "min_score", "singlet_score", "conv_all", "conv_doublet")

# the table may be staged under a prefixed filename, so fall back to a glob
if (!nzchar(RCTD) || !file.exists(RCTD)) {
  hits <- Sys.glob(c("*rctd_results_EX.csv.gz", "./*rctd_results_EX.csv.gz"))
  if (length(hits)) RCTD <- hits[1]
}
stopifnot(file.exists(OBJ), file.exists(RCTD))
cat("RCTD table: ", RCTD, "\n", sep = "")
obj <- readRDS(OBJ)
n0  <- ncol(obj)
cat("object: ", OBJ, "\ncells: ", n0, "\n", sep = "")
cat("columns before: ", paste(colnames(obj@meta.data), collapse = ", "), "\n", sep = "")

rr <- read.csv(gzfile(RCTD), stringsAsFactors = FALSE)
rr <- rr[!duplicated(rr$cell), ]
rownames(rr) <- rr$cell
cat("RCTD rows: ", nrow(rr), "\n", sep = "")

present <- intersect(RCTD_COLS, colnames(rr))
for (cl in present) obj[[cl]] <- rr[colnames(obj), cl]
cat("attached: ", paste(present, collapse = ", "), "\n", sep = "")

# verification --------------------------------------------------------------
stopifnot(ncol(obj) == n0)                       # no cells gained or lost
matched <- sum(colnames(obj) %in% rownames(rr))
cat("cells matched to an RCTD row: ", matched, " (",
    sprintf("%.1f%%", 100 * matched / n0), ")\n", sep = "")
cat("spot_class table:\n"); print(table(obj$spot_class, useNA = "ifany"))
# first_type must reproduce the predicted.celltype the pipeline stored
agree <- mean(as.character(obj$first_type) ==
              as.character(obj$predicted.celltype), na.rm = TRUE)
cat(sprintf("first_type vs stored predicted.celltype: %.3f%% agreement\n",
            100 * agree))

saveRDS(obj, OBJ)
cat("rewritten: ", OBJ, "\n", sep = "")
cat("columns after: ", paste(colnames(obj@meta.data), collapse = ", "), "\n", sep = "")
cat("DONE\n")
