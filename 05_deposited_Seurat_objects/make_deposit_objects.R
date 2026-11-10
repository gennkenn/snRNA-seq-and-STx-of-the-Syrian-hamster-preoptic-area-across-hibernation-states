#!/usr/bin/env Rscript
# =============================================================================
# Prepare the Seurat objects for public deposition.
#
# For each object this script
#   * keeps only the metadata columns a re-user needs (see KEEP below),
#     dropping intermediate clustering resolutions that were never used for
#     the published annotation;
#   * keeps ALL dimensional reductions (pca, integrated.dr, umap) -- the UMAP
#     coordinates are the only way to reproduce the published embeddings, and
#     integrated.dr is the CCA space the cross-species comparison was built in;
#   * keeps the Xenium @images (spatial coordinates), without which the
#     spatial data is unusable;
#   * for the Xenium objects, adds the full RCTD results_df columns
#     (spot_class, second_type, min_score, singlet_score, convergence flags)
#     that the original pipeline discarded when it kept only first_type;
#   * writes the object under a NEW name in OUTDIR. Input files are opened
#     read-only and never modified.
#
# It also emits deposit_metadata_dictionary.csv describing every retained
# column, for the Data Records section.
#
# Usage:
#   OUTDIR=~/claude/POA_deposit RCTDDIR=~/claude/POA_cluster_correspondence/rctd \
#     Rscript make_deposit_objects.R
# =============================================================================

suppressPackageStartupMessages({ library(Seurat) })

OUTDIR  <- Sys.getenv("OUTDIR",  "~/claude/POA_deposit")
RCTDDIR <- Sys.getenv("RCTDDIR", "")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# What to keep, per dataset family
# ---------------------------------------------------------------------------
KEEP <- list(
  crossspecies = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
                   "Species", "Sex", "celltype_display", "NeuronType"),
  crosscondition = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
                     "nCount_SCT", "nFeature_SCT", "Mode", "Sex",
                     "seurat_clusters"),
  xenium = c("orig.ident", "nCount_Xenium", "nFeature_Xenium",
             "nCount_BlankCodeword", "nCount_ControlCodeword",
             "nCount_ControlProbe", "Condition", "animal",
             "seurat_clusters", "predicted.celltype", "AVPe_MnPO_MPN")
)

DESCRIBE <- c(
  orig.ident             = "Sequencing/imaging run identifier",
  nCount_RNA             = "Total UMI count per nucleus (snRNA-seq)",
  nFeature_RNA           = "Number of genes detected per nucleus",
  percent.mt             = "Percentage of mitochondrial transcripts (QC)",
  nCount_SCT             = "Total count after SCTransform normalisation",
  nFeature_SCT           = "Number of genes after SCTransform normalisation",
  Species                = "mouse or hamster (cross-species dataset only)",
  Sex                    = "Sex of the donor animal",
  Mode                   = "summer (non-hibernating) or winter (hibernating)",
  Condition              = "Ctrl (non-hibernating) or Winter (hibernating)",
  animal                 = "Donor animal identifier",
  celltype_display       = "Published neuronal subtype label (e1-e44 / i1-i48)",
  NeuronType             = "Excitatory (Ex) or inhibitory (In) assignment",
  seurat_clusters        = "Published cluster identity used in the figures",
  predicted.celltype     = "RCTD first_type: cross-species subtype assigned to this cell",
  nCount_Xenium          = "Total transcript count per cell (Xenium)",
  nFeature_Xenium        = "Number of panel genes detected per cell",
  nCount_BlankCodeword   = "Negative control: blank codeword counts",
  nCount_ControlCodeword = "Negative control: control codeword counts",
  nCount_ControlProbe    = "Negative control: control probe counts",
  AVPe_MnPO_MPN          = "TRUE if the cell lies in the AVPe/MnPO/MPN region",
  spot_class             = "RCTD call: singlet, doublet_certain, doublet_uncertain or reject",
  first_type             = "RCTD best-fitting cross-species subtype",
  second_type            = "RCTD second cell type in the doublet fit",
  min_score              = "RCTD negative log-likelihood of the best fit",
  singlet_score          = "RCTD fit score under the singlet assumption",
  conv_all               = "RCTD convergence flag, all-type fit",
  conv_doublet           = "RCTD convergence flag, doublet fit"
)

JOBS <- list(
  list(out = "GSE345405_MPN_crossspecies_vGluT2_deposit.rds",
       family = "crossspecies", rctd = NA,
       src = "~/claude/STx_hamster/MPN_merged.neuron_vGluT2_QC.rds"),
  list(out = "GSE345405_MPN_crossspecies_vGat_deposit.rds",
       family = "crossspecies", rctd = NA,
       src = "~/claude/STx_hamster/MPN_merged.neuron_vGat_QC.rds"),
  list(out = "GSE346081_sumwin_vGluT2_deposit.rds",
       family = "crosscondition", rctd = NA,
       src = Sys.getenv("SUMWIN_EX", "~/claude/sumwinHam_merged.neuron_vGluT2.rds")),
  list(out = "GSE346081_sumwin_vGat_deposit.rds",
       family = "crosscondition", rctd = NA,
       src = Sys.getenv("SUMWIN_IN", "~/claude/sumwinHam_merged.neuron_vGat_afterQC.rds")),
  list(out = "GSE345359_Xenium_EX_deposit.rds",
       family = "xenium", rctd = "EX",
       src = "~/claude/STx_hamster/merged.obj_EX_RCTD.rds"),
  list(out = "GSE345359_Xenium_IN_deposit.rds",
       family = "xenium", rctd = "IN",
       src = "~/claude/STx_hamster/merged.obj_IN_RCTD.rds")
)

RCTD_COLS <- c("spot_class", "first_type", "second_type",
               "min_score", "singlet_score", "conv_all", "conv_doublet")

dict <- list()

for (J in JOBS) {
  src <- path.expand(J$src)
  if (!file.exists(src)) { cat("\n!! MISSING, skipped: ", src, "\n", sep = ""); next }
  cat("\n", strrep("=", 70), "\n", J$out, "\n", sep = "")
  obj <- readRDS(src)
  n0 <- ncol(obj)
  md0 <- colnames(obj@meta.data)

  # ---- attach the RCTD columns the original pipeline dropped --------------
  if (!is.na(J$rctd) && nzchar(RCTDDIR)) {
    f <- file.path(path.expand(RCTDDIR), sprintf("rctd_results_%s.csv.gz", J$rctd))
    if (file.exists(f)) {
      rr <- read.csv(gzfile(f), stringsAsFactors = FALSE)
      rr <- rr[!duplicated(rr$cell), ]
      rownames(rr) <- rr$cell
      for (cl in intersect(RCTD_COLS, colnames(rr))) {
        obj[[cl]] <- rr[colnames(obj), cl]
      }
      cat("RCTD columns attached: ",
          paste(intersect(RCTD_COLS, colnames(rr)), collapse = ", "), "\n", sep = "")
    } else cat("!! RCTD results not found at ", f, " -- columns not attached\n", sep = "")
  }

  keep <- intersect(c(KEEP[[J$family]], RCTD_COLS), colnames(obj@meta.data))
  dropped <- setdiff(colnames(obj@meta.data), keep)
  obj@meta.data <- obj@meta.data[, keep, drop = FALSE]

  cat("cells: ", n0, " -> ", ncol(obj), "\n", sep = "")
  cat("kept  (", length(keep), "): ", paste(keep, collapse = ", "), "\n", sep = "")
  cat("dropped (", length(dropped), "): ", paste(dropped, collapse = ", "), "\n", sep = "")
  cat("reductions kept: ", paste(Reductions(obj), collapse = ", "), "\n", sep = "")
  cat("images kept: ", paste(Images(obj), collapse = ", "), "\n", sep = "")
  stopifnot(ncol(obj) == n0)

  for (cl in keep) {
    v <- obj@meta.data[[cl]]
    dict[[length(dict) + 1]] <- data.frame(
      file = J$out, column = cl, type = class(v)[1],
      n_unique = length(unique(v)), n_missing = sum(is.na(v)),
      example = paste(head(unique(as.character(v)), 3), collapse = " | "),
      description = ifelse(cl %in% names(DESCRIBE), DESCRIBE[[cl]], ""),
      stringsAsFactors = FALSE)
  }

  saveRDS(obj, file.path(path.expand(OUTDIR), J$out))
  cat("written: ", file.path(OUTDIR, J$out), "\n", sep = "")
  rm(obj); invisible(gc())
}

dd <- do.call(rbind, dict)
write.csv(dd, file.path(path.expand(OUTDIR), "deposit_metadata_dictionary.csv"),
          row.names = FALSE)
cat("\ndictionary rows: ", nrow(dd), "\n", sep = "")
cat("DONE\n")
