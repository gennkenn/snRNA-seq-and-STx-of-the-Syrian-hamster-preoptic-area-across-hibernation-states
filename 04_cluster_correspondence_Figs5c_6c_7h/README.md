# 04 — Correspondence between the three cluster systems (Figs. 5c, 6c, 7h)

The three datasets were clustered independently and carry three label systems:

| Dataset | Labels | Metadata column |
|---|---|---|
| Cross-species atlas (folder 01) | e1–e44, i1–i48 | `celltype_display` |
| Cross-condition snRNA-seq (folder 02) | excitatory 0–50, inhibitory 0–52 | `seurat_clusters` |
| Xenium STx (folder 03) | excitatory 0–43, inhibitory 0–56 | `seurat_clusters` |

The cross-species subtypes (e/i) serve as the common vocabulary. The two sub-folders map the other two systems onto it:

* `label_transfer_Fig5c_6c/` — cross-condition cluster ↔ cross-species subtype, by bidirectional Seurat label transfer (snRNA-seq to snRNA-seq, ~13–14 k shared genes). Produces Fig. 5c (excitatory) and Fig. 6c (inhibitory).
* `xenium_rctd_Fig7h/` — Xenium cluster ↔ cross-species subtype, from the RCTD assignment of each Xenium cell (100-gene panel; only cells RCTD classified as `singlet` are used). Produces Fig. 7h.
* `source_data/` — every value plotted in Figs. 5c, 6c and 7h, plus per-cluster summary tables and one master table.

## Vocabulary used in all tables

* **forward** (row-normalised): for a query cluster, the fraction of its cells assigned to each atlas subtype → identity of the cluster.
* **reverse** (column-normalised): for an atlas subtype, the fraction of its assigned cells that fall in each query cluster → where the subtype went.
* **reciprocal best**: the pair is the best match in both directions. **one-directional best**: best in only one direction.
* **relation**: `1:1`, `1:many` (one query cluster collects several subtypes), `many:1` (several query clusters share one subtype), `unresolved` (no dominant subtype).
* **frac_best**: fraction of the query cluster's cells carrying the best-matching label. **median score**: median `prediction.score.max` of the cluster (label transfer only). **frac_singlet**: fraction of the Xenium cluster that RCTD called a confident singlet (RCTD only).

## Headline numbers (current tables)

| | Excitatory | Inhibitory |
|---|---|---|
| snRNA-seq ↔ snRNA-seq: reciprocal best / query clusters | 31 / 51 | 30 / 53 |
| relation classes 1:1 / 1:many / many:1 / unresolved | 20 / 15 / 12 / 4 | 24 / 13 / 11 / 5 |
| Xenium ↔ atlas: reciprocal best / Xenium clusters | 27 / 44 | 30 / 57 |

## source_data/

| File | Used in |
|---|---|
| `sourcedata_fig_correspondence_Ex.csv`, `..._Inh.csv` | Fig. 5c, Fig. 6c (long-form: every tile value, best-match flags, confidence bar) |
| `sourcedata_fig_correspondence_STx_EX.csv`, `..._STx_IN.csv` | Fig. 7h |
| `correspondence_table_Ex.csv`, `..._Inh.csv` | per cross-condition cluster: best atlas subtype, frac_best, median score, reciprocal flag, relation class |
| `correspondence_table_STx_EX.csv`, `..._STx_IN.csv` | per Xenium cluster: best atlas subtype, frac_best, frac_singlet, reciprocal flag, relation class |
| `correspondence_master_table.csv` | all four tables stacked (205 rows), one row per query cluster |
