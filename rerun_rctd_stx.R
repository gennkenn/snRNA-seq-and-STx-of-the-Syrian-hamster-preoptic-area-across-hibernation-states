#!/usr/bin/env Rscript
# =============================================================================
# Re-run RCTD class-matched on the Xenium POA data and verify the result
# against the predicted.celltype already stored in the objects.
#
#   excitatory : merged.obj_EX_RCTD.rds  <- MPN_merged.neuron_vGluT2_QC.rds
#   inhibitory : merged.obj_IN_RCTD.rds  <- MPN_merged.neuron_vGat_QC.rds
#
# The procedure replicates Fig7h_RCTD_hamster.R from the project repository:
#   * the Xenium assay is joined and its full counts matrix is passed to
#     SpatialRNA together with the centroids of ONE FOV at a time;
#   * one RCTD run per FOV (fov, fov.Ctrl2, fov.Winter1, fov.Winter2);
#   * create.RCTD(..., max_cores) then run.RCTD(doublet_mode = "doublet");
#   * the per-FOV results_df are concatenated.
# Only max_cores differs from the original (8); RCTD fits each cell
# independently, so this affects runtime and not the fitted values.
#
# Two deliberate deviations from the original script, both recorded here:
#   1. Reference files. The original read
#      integration_CCA_LogNormalize/MPN_merged.neuron_vGluT2_afterQC.rds and
#      ..._vGat_afterQC.rds. This run uses the *_QC.rds files supplied in
#      STXDIR (the '_official' variants carrying the published
#      celltype_display labels). Any disagreement with the stored
#      predicted.celltype is expected to come from here first.
#   2. The original assigned merged.obj$predicted.celltype once per FOV,
#      each assignment overwriting the previous one, before finally
#      concatenating all four. Only the final concatenation is reproduced.
#
# Unlike the stored objects, the full results_df is kept, so spot_class,
# second_type, min_score, singlet_score and the convergence flags become
# available for downstream confidence filtering.
#
# Outputs (OUTDIR)
#   rctd_results_<CLS>.csv.gz        per-cell results_df + fov, all FOVs
#   rctd_results_<CLS>_<FOV>.csv.gz  per-FOV intermediate (written as it goes)
#   rctd_vs_stored_<CLS>.csv         per-cell comparison with predicted.celltype
#   rctd_agreement_<CLS>.txt         agreement summary
#   rctd_spotclass_<CLS>.csv         spot_class x seurat_clusters counts
#   sessioninfo.txt
#
# Usage:
#   STXDIR=~/claude/STx_hamster OUTDIR=./out MAXCORES=16 Rscript rerun_rctd_stx.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(Matrix)
})

STXDIR   <- Sys.getenv("STXDIR", ".")
OUTDIR   <- Sys.getenv("OUTDIR", "./out")
MAXCORES <- as.integer(Sys.getenv("MAXCORES", "16"))
# Which neuronal classes to run, and where to look for per-FOV results that a
# previous run already produced. A FOV whose CSV is found in PRECOMPUTED_DIR is
# loaded instead of re-fitted, so an interrupted run resumes where it stopped
# rather than repeating hours of identical computation.
CLASSES  <- strsplit(Sys.getenv("CLASSES", "EX,IN"), ",")[[1]]
PRECOMP  <- Sys.getenv("PRECOMPUTED_DIR", "")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat(R.version.string, "\n")
cat(sprintf("Seurat %s | SeuratObject %s | spacexr %s | max_cores %d\n",
            packageVersion("Seurat"), packageVersion("SeuratObject"),
            packageVersion("spacexr"), MAXCORES))

gzw <- function(df, path) {
  con <- gzfile(path, "wt"); write.csv(df, con, row.names = FALSE); close(con)
}

JOBS <- list(
  list(cls = "EX", obj = "merged.obj_EX_RCTD.rds", ref = "MPN_merged.neuron_vGluT2_QC.rds"),
  list(cls = "IN", obj = "merged.obj_IN_RCTD.rds", ref = "MPN_merged.neuron_vGat_QC.rds")
)

JOBS <- Filter(function(J) J$cls %in% CLASSES, JOBS)
cat("classes to run: ", paste(vapply(JOBS, function(J) J$cls, ""), collapse = ", "),
    if (nzchar(PRECOMP)) sprintf(" | reusing per-FOV results from %s", PRECOMP) else "",
    "\n", sep = "")

for (J in JOBS) {
  cls <- J$cls
  cat("\n", strrep("=", 30), " ", cls, " ", strrep("=", 30), "\n", sep = "")
  t0 <- Sys.time()

  # ---- query ---------------------------------------------------------------
  obj <- readRDS(file.path(STXDIR, J$obj))
  obj@active.assay <- "Xenium"
  obj <- JoinLayers(obj)
  query.counts <- GetAssayData(obj, assay = "Xenium", layer = "counts")
  cat(sprintf("query object: %d cells x %d genes | FOVs: %s\n",
              ncol(obj), nrow(query.counts), paste(Images(obj), collapse = ", ")))

  stored <- as.character(obj$predicted.celltype)
  names(stored) <- colnames(obj)
  stored_clusters <- as.character(obj$seurat_clusters)
  names(stored_clusters) <- colnames(obj)

  # ---- reference -----------------------------------------------------------
  ref <- readRDS(file.path(STXDIR, J$ref))
  Idents(ref) <- "celltype_display"
  ref@active.assay <- "RNA"
  ref.counts <- GetAssayData(ref, assay = "RNA", layer = "counts")
  cluster <- as.factor(ref$celltype_display)
  names(cluster) <- colnames(ref)
  nUMI <- colSums(ref.counts)
  names(nUMI) <- colnames(ref)
  levels(cluster) <- gsub("/", "-", levels(cluster))
  reference <- Reference(ref.counts, cluster, nUMI)
  cat(sprintf("reference: %d cells, %d cell types, %d genes\n",
              ncol(ref.counts), nlevels(cluster), nrow(ref.counts)))
  cat("reference species composition: ")
  print(table(ref$Species))
  rm(ref, ref.counts); invisible(gc())

  # ---- one RCTD run per FOV ------------------------------------------------
  # The objects carry extra FOVs beyond the four slides -- a cropped `zoom`
  # region and a `combined_*` view that both re-present cells already covered
  # by the slide FOVs. The original script used only fov, fov.Ctrl2,
  # fov.Winter1 and fov.Winter2, so the loop is restricted to those; including
  # the others would roughly double the runtime and contribute only duplicate
  # cells that the de-duplication below would discard anyway.
  fovs <- grep("^fov", Images(obj), value = TRUE)
  cat("FOVs used: ", paste(fovs, collapse = ", "),
      " (skipped: ", paste(setdiff(Images(obj), fovs), collapse = ", "), ")\n", sep = "")

  res_list <- list()
  for (f in fovs) {
    cat("\n-- FOV", f, "--\n")

    pre <- if (nzchar(PRECOMP))
      file.path(PRECOMP, sprintf("rctd_results_%s_%s.csv.gz", cls, f)) else ""
    if (nzchar(pre) && file.exists(pre)) {
      df <- read.csv(gzfile(pre), stringsAsFactors = FALSE)
      cat(sprintf("  reused precomputed result: %d cells\n", nrow(df)))
      print(table(df$spot_class))
      gzw(df, file.path(OUTDIR, sprintf("rctd_results_%s_%s.csv.gz", cls, f)))
      res_list[[f]] <- df
      next
    }

    coords <- GetTissueCoordinates(obj[[f]], which = "centroids")
    rownames(coords) <- coords$cell
    coords$cell <- NULL
    query <- SpatialRNA(coords, query.counts, colSums(query.counts))
    cat(sprintf("  SpatialRNA: %d cells\n", length(query@nUMI)))

    rctd <- create.RCTD(query, reference, max_cores = MAXCORES)
    rctd <- run.RCTD(rctd, doublet_mode = "doublet")

    df <- rctd@results$results_df
    df$cell <- rownames(df)
    df$fov  <- f
    cat(sprintf("  results_df: %d cells\n", nrow(df)))
    print(table(df$spot_class))
    gzw(df, file.path(OUTDIR, sprintf("rctd_results_%s_%s.csv.gz", cls, f)))
    res_list[[f]] <- df
    rm(rctd, query); invisible(gc())
  }

  allres <- do.call(rbind, res_list)
  rownames(allres) <- NULL
  cat(sprintf("FOV coverage: %d of %d object cells appear in at least one FOV run (%.1f%%)\n",
              length(intersect(unique(allres$cell), colnames(obj))), ncol(obj),
              100 * length(intersect(unique(allres$cell), colnames(obj))) / ncol(obj)))
  gzw(allres, file.path(OUTDIR, sprintf("rctd_results_%s.csv.gz", cls)))
  cat(sprintf("\ncombined results_df: %d rows | %d unique cells\n",
              nrow(allres), length(unique(allres$cell))))

  # ---- comparison with the stored predicted.celltype -----------------------
  # The original concatenated the four FOV vectors; with duplicated cells the
  # first occurrence wins, so the same rule is applied here.
  new <- allres$first_type[!duplicated(allres$cell)]
  new <- as.character(new)
  names(new) <- allres$cell[!duplicated(allres$cell)]

  cells <- colnames(obj)
  cmp <- data.frame(
    cell            = cells,
    seurat_clusters = unname(stored_clusters[cells]),
    stored          = unname(stored[cells]),
    rerun           = unname(new[cells]),
    stringsAsFactors = FALSE)
  sc <- allres$spot_class[!duplicated(allres$cell)]
  names(sc) <- allres$cell[!duplicated(allres$cell)]
  cmp$spot_class <- as.character(unname(sc[cells]))
  cmp$agree <- !is.na(cmp$stored) & !is.na(cmp$rerun) & cmp$stored == cmp$rerun
  write.csv(cmp, file.path(OUTDIR, sprintf("rctd_vs_stored_%s.csv", cls)),
            row.names = FALSE)

  both <- !is.na(cmp$stored) & !is.na(cmp$rerun)
  lines <- c(
    sprintf("class: %s", cls),
    sprintf("reference used: %s", J$ref),
    sprintf("cells in object: %d", nrow(cmp)),
    sprintf("stored predicted.celltype non-NA: %d (%.1f%%)",
            sum(!is.na(cmp$stored)), 100 * mean(!is.na(cmp$stored))),
    sprintf("rerun first_type non-NA:          %d (%.1f%%)",
            sum(!is.na(cmp$rerun)), 100 * mean(!is.na(cmp$rerun))),
    sprintf("both non-NA: %d", sum(both)),
    sprintf("IDENTICAL on both-non-NA: %d (%.3f%%)",
            sum(cmp$agree), 100 * sum(cmp$agree) / max(sum(both), 1)),
    sprintf("stored NA but rerun assigned: %d", sum(is.na(cmp$stored) & !is.na(cmp$rerun))),
    sprintf("rerun NA but stored assigned: %d", sum(!is.na(cmp$stored) & is.na(cmp$rerun))),
    "",
    "agreement by spot_class:")
  tb <- table(cmp$spot_class[both], cmp$agree[both])
  lines <- c(lines, capture.output(print(tb)), "",
             "top 15 disagreeing stored -> rerun pairs:")
  dis <- cmp[both & !cmp$agree, ]
  if (nrow(dis)) {
    tt <- sort(table(paste(dis$stored, "->", dis$rerun)), decreasing = TRUE)
    lines <- c(lines, capture.output(print(head(tt, 15))))
  } else lines <- c(lines, "(none)")
  writeLines(lines, file.path(OUTDIR, sprintf("rctd_agreement_%s.txt", cls)))
  cat(paste(lines, collapse = "\n"), "\n")

  sct <- as.data.frame(table(spot_class = cmp$spot_class,
                             seurat_clusters = cmp$seurat_clusters))
  write.csv(sct, file.path(OUTDIR, sprintf("rctd_spotclass_%s.csv", cls)),
            row.names = FALSE)

  cat(sprintf("\n%s done in %.1f min\n", cls,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(obj, query.counts, reference, allres, res_list, cmp); invisible(gc())
}

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessioninfo.txt"))
cat("\nALL DONE\n")
