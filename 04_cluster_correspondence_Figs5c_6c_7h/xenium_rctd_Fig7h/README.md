# Xenium clusters ↔ cross-species subtypes via RCTD (Fig. 7h)

Command-line scripts (`Rscript`). One run per neuronal class (`CLASS=EX` or `CLASS=IN`). Run in this order:

| # | Script | Input | Output |
|---|---|---|---|
| 1 | `rerun_rctd_stx.R` | `merged.obj_EX_RCTD.rds` / `merged.obj_IN_RCTD.rds` (folder 03); reference `MPN_merged.neuron_vGluT2_QC.rds` / `MPN_merged.neuron_vGat_QC.rds` (folder 01, all nuclei) | RCTD (`spacexr`, doublet mode, one run per FOV) → `rctd_results_<CLASS>.csv.gz` (per cell: `first_type`, `spot_class`, …) and `rctd_vs_stored_<CLASS>.csv` comparing with the `predicted.celltype` stored in the object. Resume-capable (per-FOV checkpoints). |
| 2 | `rctd_counts_singlet.R` | `rctd_vs_stored_<CLASS>.csv` | `counts_<CLASS>.csv` — Xenium cluster × atlas subtype contingency table restricted to RCTD **singlets** (non-singlets are counted in an `NA` column); `spotclass_by_cluster_<CLASS>.csv` |
| 3 | `cluster_correspondence_STx.R` | `counts_<CLASS>.csv` | forward/reverse matrices, `correspondence_table_STx_<CLASS>.csv`, `sourcedata_fig_correspondence_STx_<CLASS>.csv`, panel `fig_correspondence_STx_rowNumberColSorted_<CLASS>.eps/.png` (= Fig. 7h) |
| — | `extract_stx_contingency.R` | the Xenium objects | alternative to step 1–2: contingency table of `seurat_clusters` × stored `predicted.celltype` using **all** cells (no singlet filter); kept for reference, not used for Fig. 7h |

Why singlets only: RCTD returns a `first_type` for every cell it processes, including cells whose `spot_class` is `reject`; building the correspondence on confident singlets only gives an interpretable per-cluster confidence (`frac_singlet`, the bar above the Fig. 7h heatmaps).
