# 05 — Construction of the deposited Seurat objects

Command-line scripts (`Rscript`) that turn the analysis objects of folders 01–03 into the six `.rds` files deposited with the GEO records, and attach the cluster-correspondence results (folder 04) to their metadata. Run in this order:

| # | Script | What it does |
|---|---|---|
| 1 | `make_deposit_objects.R` | Reads each analysis object, keeps the RNA counts, embeddings and curated metadata, renames metadata columns to a common dictionary, and writes the object under its deposit name. Writes `deposit_metadata_dictionary.csv` (one row per metadata column per object: name, type, meaning). |
| 2 | `attach_rctd_EX.R` | Attaches the RCTD re-run results (`rctd_results_EX.csv.gz` from folder 04) to the excitatory Xenium deposit object as metadata columns `spot_class`, `first_type`, `second_type`, `min_score`, `singlet_score`, `conv_all`, `conv_doublet` (the inhibitory object received the same columns during step 1). |
| 3 | `attach_labeltransfer.R` | Attaches the label-transfer results (folder 04, `percell_predictions_query_<CLASS>.csv`) to the two cross-condition deposit objects as `crossspecies_subtype_transferred` and `crossspecies_transfer_score`. |
| 4 | `sync_predicted_celltype.R` | Sets `predicted.celltype` in the Xenium deposit objects equal to the RCTD `first_type`, so that the deposited label is exactly the one used for Fig. 7h. |

| Deposit file | GEO | Source object (folder) |
|---|---|---|
| `GSE345405_MPN_crossspecies_vGluT2_deposit.rds` | GSE345405 | `MPN_merged.neuron_vGluT2_QC.rds` (01) |
| `GSE345405_MPN_crossspecies_vGat_deposit.rds` | GSE345405 | `MPN_merged.neuron_vGat_QC.rds` (01) |
| `GSE346081_sumwin_vGluT2_deposit.rds` | GSE346081 | `sumwinHam_merged.neuron_vGluT2.rds` (02) |
| `GSE346081_sumwin_vGat_deposit.rds` | GSE346081 | `sumwinHam_merged.neuron_vGat_afterQC.rds` (02) |
| `GSE345359_Xenium_EX_deposit.rds` | GSE345359 | `merged.obj_EX_RCTD.rds` (03) |
| `GSE345359_Xenium_IN_deposit.rds` | GSE345359 | `merged.obj_IN_RCTD.rds` (03) |

`deposit_metadata_dictionary.csv` in this folder documents every metadata column of the deposited objects.
