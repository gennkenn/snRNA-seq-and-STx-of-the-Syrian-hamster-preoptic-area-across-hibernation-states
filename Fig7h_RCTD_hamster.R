MPN_merged.neuron_vGluT2_afterQC<-readRDS("~/singlecellRNAseq/hamster_mouse_nuclei/integration_CCA_LogNormalize/MPN_merged.neuron_vGluT2_afterQC.rds")
MPN_merged.neuron_vGat_afterQC<-readRDS("~/singlecellRNAseq/hamster_mouse_nuclei/integration_CCA_LogNormalize/MPN_merged.neuron_vGat_afterQC.rds")


# RCTD
merged.obj_EX@active.assay='Xenium'
merged.obj_EX <- JoinLayers(merged.obj_EX)
query.counts <- GetAssayData(merged.obj_EX, assay = "Xenium", slot = "counts")

# ctrl1
coords <- GetTissueCoordinates(merged.obj_EX[["fov"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))


# remove CR cells because there aren't enough of them for annotation
Idents(MPN_merged.neuron_vGluT2_afterQC) <- "celltype_display"
MPN_merged.neuron_vGluT2_afterQC@active.assay='RNA'
counts <- GetAssayData(MPN_merged.neuron_vGluT2_afterQC, assay = "RNA", slot = "counts")
cluster <- as.factor(MPN_merged.neuron_vGluT2_afterQC$celltype_display)
names(cluster) <- colnames(MPN_merged.neuron_vGluT2_afterQC)
nUMI <- MPN_merged.neuron_vGluT2_afterQC$nCount_RNA
names(nUMI) <- colnames(MPN_merged.neuron_vGluT2_afterQC)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD.ctrl1_Ex <- create.RCTD(query, reference, max_cores = 8)
RCTD.ctrl1_Ex <- run.RCTD(RCTD.ctrl1_Ex, doublet_mode = "doublet")

annotations.df <- RCTD.ctrl1_Ex@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
merged.obj_EX$predicted.celltype <- annotations

query.counts <- GetAssayData(merged.obj_EX, assay = "Xenium", slot = "counts")


# ctrl2
coords <- GetTissueCoordinates(merged.obj_EX[["fov.Ctrl2"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# run RCTD with many cores
RCTD.ctrl2_Ex <- create.RCTD(query, reference, max_cores = 8)
RCTD.ctrl2_Ex <- run.RCTD(RCTD.ctrl2_Ex, doublet_mode = "doublet")

annotations.df <- RCTD.ctrl2_Ex@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
merged.obj_EX$predicted.celltype <- annotations

# winter1
query.counts <- GetAssayData(merged.obj_EX, assay = "Xenium", slot = "counts")
coords <- GetTissueCoordinates(merged.obj_EX[["fov.Winter1"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))
# run RCTD with many cores
RCTD.winter1_Ex <- create.RCTD(query, reference, max_cores = 8)
RCTD.winter1_Ex <- run.RCTD(RCTD.winter1_Ex, doublet_mode = "doublet")
annotations.df <- RCTD.winter1_Ex@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
merged.obj_EX$predicted.celltype <- annotations

# winter2
query.counts <- GetAssayData(merged.obj_EX, assay = "Xenium", slot = "counts")
coords <- GetTissueCoordinates(merged.obj_EX[["fov.Winter2"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))
# run RCTD with many cores
RCTD.winter2_Ex <- create.RCTD(query, reference, max_cores = 8)
RCTD.winter2_Ex <- run.RCTD(RCTD.winter2_Ex, doublet_mode = "doublet")

annotations.df <- RCTD.winter2_Ex@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
merged.obj_EX$predicted.celltype <- annotations


# inhibitory neurons
merged.obj_IN@active.assay='Xenium'
merged.obj_IN <- JoinLayers(merged.obj_IN)
query.counts <- GetAssayData(merged.obj_IN, assay = "Xenium", slot = "counts")
# crop <- Crop(merged.obj_ViscMN_afterQC[["fov"]])
# merged.obj_ViscMN_afterQC[["crop"]] <- crop
coords <- GetTissueCoordinates(merged.obj_IN[["fov"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))


Idents(MPN_merged.neuron_vGat_afterQC) <- "celltype_display"
MPN_merged.neuron_vGat_afterQC@active.assay='RNA'
counts <- GetAssayData(MPN_merged.neuron_vGat_afterQC, assay = "RNA", slot = "counts")
cluster <- as.factor(MPN_merged.neuron_vGat_afterQC$celltype_display)
names(cluster) <- colnames(MPN_merged.neuron_vGat_afterQC)
nUMI <- MPN_merged.neuron_vGat_afterQC$nCount_RNA
names(nUMI) <- colnames(MPN_merged.neuron_vGat_afterQC)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)


# run RCTD with many cores
RCTD.ctrl1_In <- create.RCTD(query, reference, max_cores = 8)
RCTD.ctrl1_In <- run.RCTD(RCTD.ctrl1_In, doublet_mode = "doublet")

annotations.df <- RCTD.ctrl1_In@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
merged.obj_IN$predicted.celltype <- annotations


# ctrl2
query.counts <- GetAssayData(merged.obj_IN, assay = "Xenium", slot = "counts")
# crop <- Crop(merged.obj_ViscMN_afterQC[["fov"]])
# merged.obj_ViscMN_afterQC[["crop"]] <- crop
coords <- GetTissueCoordinates(merged.obj_IN[["fov.Ctrl2"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# run RCTD with many cores
RCTD.ctrl2_In <- create.RCTD(query, reference, max_cores = 8)
RCTD.ctrl2_In <- run.RCTD(RCTD.ctrl2_In, doublet_mode = "doublet")

annotations.df <- RCTD.ctrl2_In@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
merged.obj_IN$predicted.celltype <- annotations

# winter1
query.counts <- GetAssayData(merged.obj_IN, assay = "Xenium", slot = "counts")
# crop <- Crop(merged.obj_ViscMN_afterQC[["fov"]])
# merged.obj_ViscMN_afterQC[["crop"]] <- crop
coords <- GetTissueCoordinates(merged.obj_IN[["fov.Winter1"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

Idents(MPN_merged.neuron_vGat_afterQC) <- "celltype_display"
# MPN_merged.neuron_vGat_afterQC_filtered <- subset(MPN_merged.neuron_vGat_afterQC, idents = c("Int-i84"), invert = TRUE)
# "Int-i0","Int-i6","Int-i7","Int-i8","Int-i10","Int-i12","Int-i14","Int-i15","Int-i19","Int-i20","Int-i21","Int-i25","Int-i28","Int-i34","Int-i43","Int-i55","Int-i58","Int-i60","Int-i65","Int-i73","Int-i80"
MPN_merged.neuron_vGat_afterQC@active.assay='RNA'
counts <- GetAssayData(MPN_merged.neuron_vGat_afterQC, assay = "RNA", slot = "counts")
cluster <- as.factor(MPN_merged.neuron_vGat_afterQC$celltype_display)
names(cluster) <- colnames(MPN_merged.neuron_vGat_afterQC)
nUMI <- MPN_merged.neuron_vGat_afterQC$nCount_RNA
names(nUMI) <- colnames(MPN_merged.neuron_vGat_afterQC)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD.winter1_In <- create.RCTD(query, reference, max_cores = 8)
RCTD.winter1_In <- run.RCTD(RCTD.winter1_In, doublet_mode = "doublet")

annotations.df <- RCTD.winter1_In@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
merged.obj_IN$predicted.celltype <- annotations


# winter2
query.counts <- GetAssayData(merged.obj_IN, assay = "Xenium", slot = "counts")
# crop <- Crop(merged.obj_ViscMN_afterQC[["fov"]])
# merged.obj_ViscMN_afterQC[["crop"]] <- crop
coords <- GetTissueCoordinates(merged.obj_IN[["fov.Winter2"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))



# run RCTD with many cores
RCTD.winter2_In <- create.RCTD(query, reference, max_cores = 8)
RCTD.winter2_In <- run.RCTD(RCTD.winter2_In, doublet_mode = "doublet")

annotations.df <- RCTD.winter2_In@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
merged.obj_IN$predicted.celltype <- annotations

annotations.df <- RCTD.ctrl1_Ex@results$results_df
annotations1 <- annotations.df$first_type
names(annotations1) <- rownames(annotations.df)

annotations.df <- RCTD.ctrl2_Ex@results$results_df
annotations2 <- annotations.df$first_type
names(annotations2) <- rownames(annotations.df)

annotations.df <- RCTD.winter1_Ex@results$results_df
annotations3 <- annotations.df$first_type
names(annotations3) <- rownames(annotations.df)

annotations.df <- RCTD.winter2_Ex@results$results_df
annotations4 <- annotations.df$first_type
names(annotations4) <- rownames(annotations.df)

annotations.df <- RCTD.ctrl1_In@results$results_df
annotations5 <- annotations.df$first_type
names(annotations5) <- rownames(annotations.df)

annotations.df <- RCTD.ctrl2_In@results$results_df
annotations6 <- annotations.df$first_type
names(annotations6) <- rownames(annotations.df)

annotations.df <- RCTD.winter1_In@results$results_df
annotations7 <- annotations.df$first_type
names(annotations7) <- rownames(annotations.df)

annotations.df <- RCTD.winter2_In@results$results_df
annotations8 <- annotations.df$first_type
names(annotations8) <- rownames(annotations.df)

annotations_Ex<-c(annotations1,annotations2,annotations3,annotations4)
annotations_In<-c(annotations5,annotations6,annotations7,annotations8)
merged.obj_EX$predicted.celltype <- annotations_Ex
merged.obj_IN$predicted.celltype <- annotations_In



library(tidyr)
meta_df <- merged.obj_EX@meta.data %>%
  dplyr::select(predicted.celltype, seurat_clusters) %>%
  drop_na()  

links <- dplyr::count(meta_df, predicted.celltype, seurat_clusters, name = "value")
nodes <- data.frame(name = c(unique(links$predicted.celltype),
                             unique(links$seurat_clusters)))
links <- links %>%
  mutate(source = match(predicted.celltype, nodes$name) - 1,
         target = match(seurat_clusters, nodes$name) - 1)

meta_df <- merged.obj_EX@meta.data[, c("predicted.celltype", "seurat_clusters")]
meta_df$predicted.celltype <- as.character(meta_df$predicted.celltype)
meta_df$seurat_clusters <- as.character(meta_df$seurat_clusters)

meta_df <- na.omit(meta_df)
plot_df <- meta_df %>%
  group_by(predicted.celltype, seurat_clusters) %>%
  summarise(value = n(), .groups = "drop")


library(ggplot2)
library(ggalluvial)
plot_df_long <- plot_df %>%
  tidyr::pivot_longer(cols = c(predicted.celltype, seurat_clusters),
                      names_to = "axis",
                      values_to = "category")
ggplot(
  links,
  aes(
    axis1 = predicted.celltype,
    axis2 = seurat_clusters,
    y = value
  )
) +
  geom_alluvium(aes(fill = predicted.celltype), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Predicted", "Subtype"), expand = c(.1, .1)) +
  theme_minimal()

ggplot(
  links,
  aes(
    axis1 = predicted.celltype,
    axis2 = seurat_clusters,
    y = value
  )
) +
  geom_alluvium(
    aes(fill = predicted.celltype),
    width = 1/12,
    alpha = 1        #  ←これで不透明
  ) +
  geom_stratum(
    width = 1/12,
    fill = "grey80",
    color = "black",
    alpha = 1        #  ← stratums も不透明に
  ) +
  geom_label(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    alpha = 1        #  ← ラベル背景も不透明
  ) +
  scale_x_discrete(limits = c("Predicted", "Subtype"), expand = c(.1, .1)) +
  theme_minimal()


# inhibitory neuron
# Seuratオブジェクトからメタデータ抽出
meta_df <- merged.obj_IN_RCTD@meta.data %>%
  dplyr::select(predicted.celltype, seurat_clusters) %>%
  drop_na()  # NAを除去
# meta_df <- merged.obj_ViscMN_afterQC@meta.data[, c("predicted.celltype", "SPN_subtype")]
# meta_df <- na.omit(meta_df)


# 各組み合わせごとの頻度をカウント
links <- dplyr::count(meta_df, predicted.celltype, seurat_clusters, name = "value")


# ノードの一覧を作成
nodes <- data.frame(name = c(unique(links$predicted.celltype),
                             unique(links$seurat_clusters)))

# source と target のIDに変換
links <- links %>%
  mutate(source = match(predicted.celltype, nodes$name) - 1,
         target = match(seurat_clusters, nodes$name) - 1)

# 必要な列だけ抽出
meta_df <- merged.obj_IN_RCTD@meta.data[, c("predicted.celltype", "seurat_clusters")]

# factor → character（安全のため）
meta_df$predicted.celltype <- as.character(meta_df$predicted.celltype)
meta_df$seurat_clusters <- as.character(meta_df$seurat_clusters)

# NA除去
meta_df <- na.omit(meta_df)

# 組み合わせごとに頻度集計
plot_df <- meta_df %>%
  group_by(predicted.celltype, seurat_clusters) %>%
  summarise(value = n(), .groups = "drop")


library(ggplot2)
library(ggalluvial)
plot_df_long <- plot_df %>%
  tidyr::pivot_longer(cols = c(predicted.celltype, seurat_clusters),
                      names_to = "axis",
                      values_to = "category")
ggplot(
  links,
  aes(
    axis1 = predicted.celltype,
    axis2 = seurat_clusters,
    y = value
  )
) +
  geom_alluvium(aes(fill = predicted.celltype), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Predicted", "Subtype"), expand = c(.1, .1)) +
  theme_minimal()

ggplot(
  links,
  aes(
    axis1 = predicted.celltype,
    axis2 = seurat_clusters,
    y = value
  )
) +
  geom_alluvium(
    aes(fill = predicted.celltype),
    width = 1/12,
    alpha = 1        #  ←これで不透明
  ) +
  geom_stratum(
    width = 1/12,
    fill = "grey80",
    color = "black",
    alpha = 1        #  ← stratums も不透明に
  ) +
  geom_label(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    alpha = 1        #  ← ラベル背景も不透明
  ) +
  scale_x_discrete(limits = c("Predicted", "Subtype"), expand = c(.1, .1)) +
  theme_minimal()

