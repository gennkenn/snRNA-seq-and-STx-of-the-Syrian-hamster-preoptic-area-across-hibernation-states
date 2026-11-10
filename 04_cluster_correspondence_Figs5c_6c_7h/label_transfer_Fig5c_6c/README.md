# Label transfer: cross-condition clusters ↔ cross-species subtypes (Fig. 5c, Fig. 6c)

Command-line scripts (`Rscript`), configured by environment variables. One run per neuronal class (`CLASS=Ex` or `CLASS=Inh`).

| Script | Input | Output |
|---|---|---|
| `cluster_correspondence_labeltransfer.R` | REF: `MPN_merged.neuron_vGluT2_QC.rds` / `MPN_merged.neuron_vGat_QC.rds` (folder 01; labels `celltype_display`, mouse + hamster nuclei). QRY: `sumwinHam_merged.neuron_vGluT2.rds` / `sumwinHam_merged.neuron_vGat_afterQC.rds` (folder 02; labels `seurat_clusters`) | per-cell predictions (`percell_predictions_query_<CLASS>.csv`), forward/reverse matrices, `correspondence_table_<CLASS>.csv`, `reference_label_representation_<CLASS>.csv`, a first draft of the figure |
| `replot_correspondence_from_csv.R` | the CSVs above | publication panel `fig_correspondence_rowNumberColSorted_<CLASS>.eps/.png` (= Fig. 5c / 6c) and `sourcedata_fig_correspondence_<CLASS>.csv` |

Method summary (details in the script header):

1. Both objects are restricted to their shared genes and re-normalised from raw counts with one common pipeline (LogNormalize; the query's SCT assay is not used).
2. `FindTransferAnchors` / `TransferData` (PCA projection) in both directions: atlas → cross-condition (forward) and cross-condition → atlas (reverse).
3. The default reference is **all** atlas nuclei (mouse + hamster; `REF_MODE=all`). `REF_MODE=hamster` restricts the reference to hamster nuclei and was run as a sensitivity check (best match unchanged for 45/51 excitatory and 47/53 inhibitory clusters).
4. Best matches, reciprocal best hits and relation classes are derived from the forward and reverse matrices; the figure places query clusters on the x axis in native order and atlas subtypes on the y axis ordered so that best matches fall on the diagonal, with the per-cluster confidence (median prediction score) as a bar above the heatmap.

Example:

```sh
REF_RDS=MPN_merged.neuron_vGluT2_QC.rds QRY_RDS=sumwinHam_merged.neuron_vGluT2.rds \
CLASS=Ex REF_MODE=all OUTDIR=results_Ex Rscript cluster_correspondence_labeltransfer.R
OUTDIR=results_Ex CLASS=Ex Rscript replot_correspondence_from_csv.R
```
