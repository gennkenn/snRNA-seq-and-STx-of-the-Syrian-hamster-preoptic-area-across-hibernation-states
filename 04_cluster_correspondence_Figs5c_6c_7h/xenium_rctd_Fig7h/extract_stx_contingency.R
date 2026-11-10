#!/usr/bin/env Rscript
# =============================================================================
# Extract the cluster x RCTD-label contingency tables from the Xenium (STx)
# Seurat objects, so that all downstream table/figure work can be done from
# small CSVs without re-loading the multi-GB objects.
#
# Inputs  : merged.obj_EX_RCTD.rds / merged.obj_IN_RCTD.rds  (~2.0 / 2.5 GB)
# Outputs : counts_<CLASS>.csv          seurat_clusters x predicted.celltype counts
#           clustermeta_<CLASS>.csv     per-cluster n cells (+ extra columns found)
#           metainventory_<CLASS>.csv   every metadata column with its cardinality
#
# Objects are loaded one at a time and released, to keep peak memory low on a
# shared machine.
# =============================================================================

suppressPackageStartupMessages(library(Seurat))

STXDIR <- Sys.getenv("STXDIR", ".")
OUTDIR <- Sys.getenv("OUTDIR", ".")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

CLUSTER_COL <- "seurat_clusters"
LABEL_COL   <- "predicted.celltype"

for (cls in c("EX", "IN")) {
  f <- file.path(STXDIR, sprintf("merged.obj_%s_RCTD.rds", cls))
  cat("\n==================== ", cls, " ====================\n")
  cat("file: ", f, " (", round(file.size(f) / 1e9, 2), " GB)\n", sep = "")
  o <- readRDS(f)
  cat("cells: ", ncol(o), " | features: ", nrow(o), "\n", sep = "")
  cat("assays: ", paste(names(o@assays), collapse = ", "), "\n", sep = "")

  md <- o@meta.data

  # ---- inventory of every metadata column --------------------------------
  inv <- data.frame(
    column    = colnames(md),
    class     = vapply(md, function(x) class(x)[1], ""),
    n_unique  = vapply(md, function(x) length(unique(x)), 0L),
    n_missing = vapply(md, function(x) sum(is.na(x)), 0L),
    example   = vapply(md, function(x) paste(head(unique(as.character(x)), 3),
                                             collapse = " | "), ""),
    stringsAsFactors = FALSE)
  write.csv(inv, file.path(OUTDIR, sprintf("metainventory_%s.csv", cls)),
            row.names = FALSE)
  cat("\n-- metadata inventory --\n"); print(inv, right = FALSE)

  if (!all(c(CLUSTER_COL, LABEL_COL) %in% colnames(md))) {
    cat("\n!! required columns missing; skipping contingency for ", cls, "\n")
    rm(o, md); invisible(gc()); next
  }

  clu <- as.character(md[[CLUSTER_COL]])
  lab <- as.character(md[[LABEL_COL]])
  cat(sprintf("\n%s: %d levels | %s: %d levels | cells with NA label: %d\n",
              CLUSTER_COL, length(unique(clu)),
              LABEL_COL, length(unique(na.omit(lab))), sum(is.na(lab))))

  # ---- the contingency table: the sufficient statistic for everything -----
  tb <- table(cluster = clu, label = lab, useNA = "ifany")
  write.csv(as.matrix(tb), file.path(OUTDIR, sprintf("counts_%s.csv", cls)))
  cat("counts table written: ", nrow(tb), " clusters x ", ncol(tb), " labels\n",
      sep = "")

  # ---- per-cluster extras: any other low-cardinality metadata ------------
  extra <- setdiff(inv$column[inv$n_unique <= 20], c(CLUSTER_COL, LABEL_COL))
  cm <- data.frame(cluster = names(table(clu)),
                   n_cells = as.integer(table(clu)),
                   stringsAsFactors = FALSE)
  for (e in extra) {
    sub <- table(clu, as.character(md[[e]]))
    frac <- sweep(sub, 1, pmax(rowSums(sub), 1), "/")
    colnames(frac) <- paste0(e, "::", colnames(frac))
    cm <- cbind(cm, round(as.data.frame.matrix(frac)[cm$cluster, , drop = FALSE], 4))
  }
  write.csv(cm, file.path(OUTDIR, sprintf("clustermeta_%s.csv", cls)),
            row.names = FALSE)
  cat("extra per-cluster columns: ", paste(extra, collapse = ", "), "\n", sep = "")

  rm(o, md, clu, lab, tb); invisible(gc())
}

cat("\nDONE\n")
sessionInfo()
