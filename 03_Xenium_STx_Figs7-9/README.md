# 03 — Xenium spatial transcriptomics (hamster POA, 100-gene panel), Figs. 7–9

GEO: GSE345359. Four tissue sections / FOVs: `fov` (Ctrl1), `fov.Ctrl2`, `fov.Winter1`, `fov.Winter2`.

Run in this order:

| # | Script | Figure panels | Input | Output |
|---|---|---|---|---|
| 1 | `20260710_ROIextraction.Rmd` | — (preprocessing) | Xenium output bundles; ROI polygons exported from Xenium Explorer (`ROI_geojson/<animal>.csv`) | Cells inside each POA ROI are tagged with the animal id (`meta.data$animal`) using `sf` point-in-polygon tests |
| 2 | `STx_hamster_clean.R` | Fig. 7d–g, Fig. 8a, Fig. 9a, b | Xenium objects from step 1; cross-species neuron objects from folder 01 (reference for RCTD) | Gene renaming to mouse symbols, SCTransform, PCA/UMAP, clustering at several resolutions; cell-class assignment (`Cell_class`: Ex, In, MO, IO, Microglia, Astrocyte, Fibroblast, Endothelial, Ependymal, Ambiguous); neuron sub-clustering → `merged.obj_EX_RCTD.rds` (excitatory clusters 0–43) and `merged.obj_IN_RCTD.rds` (inhibitory clusters 0–56); Fig. 8a pseudobulk (per animal) IEG comparison with two-sided Welch's t-test and Benjamini–Hochberg correction across all cluster × gene tests (`PseudoBulk_WelchT_BH_Results.csv`) |
| 3 | `Fig8b.R` | Fig. 8b | `merged.obj`, `merged.obj_EX` from step 2 | Spatial map of excitatory clusters within a defined window of the Ctrl2 section |
| 4 | `Fig9cd.R` | Fig. 9c, d | `merged.obj_EX`, `merged.obj_IN` from step 2 | Fraction of Fos-positive / IEG-score-positive cells per cluster and animal (threshold = mean + 2 SD); per-cluster two-sided Welch's t-test between conditions with Benjamini–Hochberg correction across clusters; bar graphs |

The RCTD assignment of Xenium cells to cross-species subtypes (`predicted.celltype`, e1–e44 / i1–i48) is what Fig. 7h is built on; the Fig. 7h analysis itself lives in `../04_cluster_correspondence_Figs5c_6c_7h/xenium_rctd_Fig7h/`.
