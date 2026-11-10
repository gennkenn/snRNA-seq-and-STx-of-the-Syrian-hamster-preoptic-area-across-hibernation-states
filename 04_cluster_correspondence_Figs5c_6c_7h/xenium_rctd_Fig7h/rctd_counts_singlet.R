#!/usr/bin/env Rscript
# =============================================================================
# Build the Xenium cluster x atlas-subtype contingency table from the RCTD
# re-run, restricted to cells RCTD called `singlet`.
#
# Why restrict. RCTD's results_df carries a first_type for every cell it
# processed, whatever the spot_class, so the predicted.celltype column stored
# in the original objects mixes confident singlet calls with cells RCTD
# explicitly rejected. In the excitatory data that is 73,500 singlets against
# 82,626 rejects: more than half of the stored labels are calls RCTD itself
# did not trust. Correspondence built on all of them understates how cleanly
# the Xenium clusters map onto the cross-species atlas.
#
# Output format matches extract_stx_contingency.R so that
# cluster_correspondence_STx.R consumes it unchanged: a cluster x label count
# matrix whose `NA` column holds every cell that is NOT a confident singlet.
# Downstream, `frac_unassigned` therefore reads as "fraction of the cluster
# RCTD did not call a singlet", and `frac_assigned` as the singlet fraction.
#
# Input : rctd_vs_stored_<CLASS>.csv  (cell, seurat_clusters, stored, rerun,
#                                      spot_class, agree) from rerun_rctd_stx.R
# Output: counts_<CLASS>.csv
#         spotclass_by_cluster_<CLASS>.csv   spot_class composition per cluster
#
# Usage: INDIR=rctd_out OUTDIR=rctd_out CLASS=EX Rscript rctd_counts_singlet.R
# =============================================================================

INDIR  <- Sys.getenv("INDIR",  "rctd_out")
OUTDIR <- Sys.getenv("OUTDIR", "rctd_out")
CLASS  <- Sys.getenv("CLASS",  "EX")
KEEP   <- strsplit(Sys.getenv("SPOT_CLASSES", "singlet"), ",")[[1]]
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cmp <- read.csv(file.path(INDIR, sprintf("rctd_vs_stored_%s.csv", CLASS)),
                stringsAsFactors = FALSE,
                colClasses = c(seurat_clusters = "character"))
cat(sprintf("%s: %d cells in object\n", CLASS, nrow(cmp)))
print(table(cmp$spot_class, useNA = "ifany"))

keep <- !is.na(cmp$spot_class) & cmp$spot_class %in% KEEP & !is.na(cmp$rerun)
cat(sprintf("kept as confident (%s): %d (%.1f%%)\n",
            paste(KEEP, collapse = "/"), sum(keep), 100 * mean(keep)))

clu_all <- sort(unique(as.integer(cmp$seurat_clusters)))
clu_all <- as.character(clu_all)
labs    <- sort(unique(cmp$rerun[keep]))

tb <- table(factor(cmp$seurat_clusters[keep], levels = clu_all),
            factor(cmp$rerun[keep],           levels = labs))
M <- as.matrix(tb)

# every cell that is not a confident singlet goes into the NA column, so the
# downstream script's "unassigned" bookkeeping stays meaningful
n_tot  <- as.integer(table(factor(cmp$seurat_clusters, levels = clu_all)))
M <- cbind(M, `NA` = n_tot - rowSums(M))
stopifnot(all(M[, "NA"] >= 0), sum(M) == nrow(cmp))

write.csv(M, file.path(OUTDIR, sprintf("counts_%s.csv", CLASS)))
cat(sprintf("counts written: %d clusters x %d labels (+NA column)\n",
            nrow(M), ncol(M) - 1))

sc <- as.data.frame.matrix(table(factor(cmp$seurat_clusters, levels = clu_all),
                                 factor(cmp$spot_class)))
sc$n_cells <- rowSums(sc)
sc$frac_singlet <- round(sc[["singlet"]] / pmax(sc$n_cells, 1), 4)
sc <- cbind(seurat_clusters = rownames(sc), sc)
write.csv(sc, file.path(OUTDIR, sprintf("spotclass_by_cluster_%s.csv", CLASS)),
          row.names = FALSE)
cat(sprintf("singlet fraction per cluster: median %.3f, range %.3f - %.3f\n",
            median(sc$frac_singlet), min(sc$frac_singlet), max(sc$frac_singlet)))
