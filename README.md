# A single-nucleus and spatial transcriptomic atlas of the Syrian hamster preoptic area across hibernation states

Analysis code for Sato *et al.*, *Scientific Data* (Data Descriptor).
Syrian hamster (*Mesocricetus auratus*) preoptic area (POA); single-nucleus RNA-seq (snRNA-seq, 10x Genomics) and Xenium *in situ* spatial transcriptomics (STx).

## Data

| Dataset | Figures | GEO accession | Content |
|---|---|---|---|
| Cross-species snRNA-seq (hamster + mouse POA) | Figs. 1–3 | [GSE345405](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345405) | 13,233 hamster + 36,201 mouse nuclei; neuronal subtypes e1–e44 (excitatory) and i1–i48 (inhibitory) |
| Cross-condition snRNA-seq (hibernating vs non-hibernating hamster POA) | Figs. 4–6 | [GSE346081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346081) | 75,111 nuclei; excitatory clusters 0–50, inhibitory clusters 0–52 |
| Xenium spatial transcriptomics (hamster POA, 100-gene panel) | Figs. 7–9 | [GSE345359](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345359) | 224,350 excitatory / 507,222 inhibitory cells; excitatory clusters 0–43, inhibitory clusters 0–56 |

Processed Seurat objects (`.rds`) and a metadata dictionary are deposited alongside the raw data (see `05_deposited_Seurat_objects/`).
Hamster genes are named by their one-to-one mouse orthologs (mouse gene symbols) throughout; the ortholog table is `one-to-one_ortholog_genes.csv` (BCM_Maur_2.0 vs GRCm39).

## Repository layout

Folders are numbered in the order the analyses depend on each other. Each folder has its own README listing inputs and outputs of every script.

```
01_snRNAseq_cross_species_Figs1-3/       QC, integration and clustering of hamster + mouse POA nuclei
02_snRNAseq_cross_condition_Figs4-6/     hibernating vs non-hibernating hamster POA nuclei
03_Xenium_STx_Figs7-9/                   Xenium ROI extraction, clustering, RCTD, IEG analyses
04_cluster_correspondence_Figs5c_6c_7h/  correspondence of the three cluster systems (Figs. 5c, 6c, 7h)
   ├── label_transfer_Fig5c_6c/          snRNA-seq <-> snRNA-seq (Seurat label transfer)
   ├── xenium_rctd_Fig7h/                Xenium <-> cross-species atlas (RCTD)
   └── source_data/                      tables behind Figs. 5c, 6c and 7h
05_deposited_Seurat_objects/             construction of the GEO-deposited .rds objects
```

## Figure-to-script map

| Figure panel | Script | Folder |
|---|---|---|
| Fig. 1c–j | `CrossSpecies_ham_clean.R` | 01 |
| Fig. 2a–c | `CrossSpecies_ham_clean.R` | 01 |
| Fig. 2d | `Fig2d.R` | 01 |
| Fig. 3a–c | `CrossSpecies_ham_clean.R` | 01 |
| Fig. 4b–i | `hamster_sumwin_clean.R` | 02 |
| Fig. 5a, b, d–i | `hamster_sumwin_clean.R` | 02 |
| **Fig. 5c** (excitatory cluster correspondence) | `cluster_correspondence_labeltransfer.R` → `replot_correspondence_from_csv.R` (CLASS=Ex) | 04/label_transfer |
| Fig. 6a, b, d–f | `hamster_sumwin_clean.R` | 02 |
| **Fig. 6c** (inhibitory cluster correspondence) | `cluster_correspondence_labeltransfer.R` → `replot_correspondence_from_csv.R` (CLASS=Inh) | 04/label_transfer |
| Fig. 7d–g | `STx_hamster_clean.R` (after `20260710_ROIextraction.Rmd`) | 03 |
| **Fig. 7h** (Xenium cluster correspondence) | `rerun_rctd_stx.R` → `rctd_counts_singlet.R` → `cluster_correspondence_STx.R` | 04/xenium_rctd |
| Fig. 8a | `STx_hamster_clean.R` (pseudobulk IEG test, Welch's t-test + Benjamini–Hochberg) | 03 |
| Fig. 8b | `Fig8b.R` | 03 |
| Fig. 9a, b | `STx_hamster_clean.R` | 03 |
| Fig. 9c, d | `Fig9cd.R` | 03 |

Panels not listed (Fig. 1a, b; Fig. 4a; Fig. 7a–c) are schematics or microscopy images with no associated code.

> Note on panel letters inside the scripts: the `# FigNx` comments in `hamster_sumwin_clean.R` were written before the correspondence panels (Fig. 5c, Fig. 6c) were inserted. Within that script, `# Fig5c`–`# Fig5g,h` correspond to the current Fig. 5d–i and `# Fig6c`–`# Fig6d,e` to the current Fig. 6d–f.

## Order of execution

1. `01_…/CrossSpecies_ham_clean.R` builds the cross-species objects `MPN_merged.neuron_vGluT2_QC.rds` (excitatory) and `MPN_merged.neuron_vGat_QC.rds` (inhibitory) with the subtype labels `celltype_display` (e1–e44, i1–i48). These are the **reference** for everything in folders 03 and 04.
2. `02_…/hamster_sumwin_clean.R` builds the cross-condition objects `sumwinHam_merged.neuron_vGluT2.rds` and `sumwinHam_merged.neuron_vGat_afterQC.rds` (cluster ids in `seurat_clusters`).
3. `03_…` processes the Xenium runs (ROI extraction → clustering → RCTD against the folder-01 objects → IEG analyses) and builds `merged.obj_EX_RCTD.rds` / `merged.obj_IN_RCTD.rds`.
4. `04_…` relates the three cluster systems to one another (Figs. 5c, 6c, 7h). It consumes the objects from steps 1–3 and writes the tables in `04_…/source_data/`.
5. `05_…` assembles the deposited `.rds` objects, adding the correspondence results from step 4 to the metadata.

## Software

R (≥ 4.3; the correspondence and deposit scripts were run with R 4.5.1), Seurat v5, SeuratObject, spacexr (RCTD), dplyr, tidyr, ggplot2, patchwork, EnhancedVolcano, clusterProfiler, sf (ROI polygons).
Scripts in folders 04 and 05 are command-line `Rscript` programs configured through environment variables (documented in each script header and folder README); scripts in folders 01–03 are interactive analysis sessions in RStudio.

## Citation

Sato H. *et al.* A single-nucleus and spatial transcriptomic atlas of the Syrian hamster preoptic area across hibernation states. *Scientific Data* (in revision).
