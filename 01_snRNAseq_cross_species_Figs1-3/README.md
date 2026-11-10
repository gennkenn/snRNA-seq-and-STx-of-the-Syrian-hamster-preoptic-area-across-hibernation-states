# 01 — Cross-species snRNA-seq atlas (hamster + mouse POA), Figs. 1–3

GEO: GSE345405. Interactive RStudio session scripts (run top to bottom).

| Script | Figure panels | Input | Output |
|---|---|---|---|
| `CrossSpecies_ham_clean.R` | Fig. 1c–j, Fig. 2a–c, Fig. 3a–c | Cell Ranger `filtered_feature_bc_matrix` of 2 hamster and 4 mouse POA samples; `one-to-one_ortholog_genes.csv` | QC, ortholog-restricted merge, CCA integration (LogNormalize), clustering; neuron objects `MPN_merged.neuron_vGluT2_QC.rds` (excitatory, subtypes **e1–e44**) and `MPN_merged.neuron_vGat_QC.rds` (inhibitory, subtypes **i1–i48**); subtype labels in `meta.data$celltype_display` |
| `Fig2d.R` | Fig. 2d | `MPN_merged.neuron_vGluT2_QC` (and the vGat object) from the script above | Per-subtype hamster-vs-mouse differential expression (`FindMarkers`, Wilcoxon, subtypes with ≥ 3 cells of each species) and the summary plot |

The two neuron objects produced here are the **reference** used by the Xenium RCTD assignment (folder 03), by the cluster-correspondence analyses (folder 04) and as the basis of the deposited objects (folder 05).
