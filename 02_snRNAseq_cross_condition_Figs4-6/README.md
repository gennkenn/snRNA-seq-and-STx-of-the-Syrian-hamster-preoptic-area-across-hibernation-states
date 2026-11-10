# 02 — Cross-condition snRNA-seq (hibernating vs non-hibernating hamster POA), Figs. 4–6

GEO: GSE346081. Interactive RStudio session script (run top to bottom).

| Script | Figure panels | Input | Output |
|---|---|---|---|
| `hamster_sumwin_clean.R` | Fig. 4b–i; Fig. 5a, b, d–i; Fig. 6a, b, d–f | Cell Ranger `filtered_feature_bc_matrix` of 5 hamster POA samples (2 summer/non-hibernating, 3 winter/hibernating); `one-to-one_ortholog_genes.csv` | QC, SCTransform integration, clustering; neuron objects `sumwinHam_merged.neuron_vGluT2.rds` (excitatory clusters **0–50**) and `sumwinHam_merged.neuron_vGat_afterQC.rds` (inhibitory clusters **0–52**); winter-vs-summer differential expression (volcano plots), gene-set enrichment |

**Panel letters.** The `# Fig5x` / `# Fig6x` comments in this script predate the insertion of the correspondence panels Fig. 5c and Fig. 6c (folder 04). Read them as follows:

| comment in script | current panel |
|---|---|
| `# Fig5a`, `# Fig5b` | Fig. 5a, 5b |
| `# Fig5c` (winter vs summer volcano) | Fig. 5d |
| `# Fig5d` … `# Fig5f` | Fig. 5e … 5g |
| `# Fig5g,h` | Fig. 5h, 5i |
| `# Fig6a`, `# Fig6b` | Fig. 6a, 6b |
| `# Fig6c` (winter vs summer volcano) | Fig. 6d |
| `# Fig6d,e` / `# Fig6e` | Fig. 6e, 6f |

Fig. 5c and Fig. 6c themselves are produced in `../04_cluster_correspondence_Figs5c_6c_7h/label_transfer_Fig5c_6c/`.
