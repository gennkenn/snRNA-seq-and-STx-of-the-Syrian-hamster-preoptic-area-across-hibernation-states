#ライブラリの読み込みs
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(sets)
library(tibble)
library(dplyr)
library(ggplot2)
library(scales)
library(EnhancedVolcano)
library(clusterProfiler)


# Load the dataset
mouse_F1.data <- Read10X(data.dir = "~/2024_Kaplan/CR_P65F_1/filtered_feature_bc_matrix")
mouse_F2.data <- Read10X(data.dir = "~/2024_Kaplan/CR_P65F_2/filtered_feature_bc_matrix")
mouse_M1.data <- Read10X(data.dir = "~/2024_Kaplan/CR_P65M_1/filtered_feature_bc_matrix")
mouse_M2.data <- Read10X(data.dir = "~/2024_Kaplan/CR_P65M_2/filtered_feature_bc_matrix")

mouse_F1<- CreateSeuratObject(counts = mouse_F1.data, project = "mouseF1", min.cells = 3, min.features = 250)
mouse_F1
mouse_F2<- CreateSeuratObject(counts = mouse_F2.data, project = "mouseF2", min.cells = 3, min.features = 250)
mouse_F2
mouse_M1<- CreateSeuratObject(counts = mouse_M1.data, project = "mouseM1", min.cells = 3, min.features = 250)
mouse_M1
mouse_M2<- CreateSeuratObject(counts = mouse_M2.data, project = "mouseM2", min.cells = 3, min.features = 250)
mouse_M2



#QC
mouse_F1[["percent.mt"]] <- PercentageFeatureSet(mouse_F1, pattern = "^mt-")
mouse_F2[["percent.mt"]] <- PercentageFeatureSet(mouse_F2, pattern = "^mt-")
mouse_M1[["percent.mt"]] <- PercentageFeatureSet(mouse_M1, pattern = "^mt-")
mouse_M2[["percent.mt"]] <- PercentageFeatureSet(mouse_M2, pattern = "^mt-")
mouse_F1[["percent.rp"]] <- PercentageFeatureSet(mouse_F1, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
mouse_F2[["percent.rp"]] <- PercentageFeatureSet(mouse_F2, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
mouse_M1[["percent.rp"]] <- PercentageFeatureSet(mouse_M1, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
mouse_M2[["percent.rp"]] <- PercentageFeatureSet(mouse_M2, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")

# Visualize QC metrics as a violin plot
VlnPlot(mouse_F1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
VlnPlot(mouse_F2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
VlnPlot(mouse_M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
VlnPlot(mouse_M2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)

plot1 <- FeatureScatter(mouse_F1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_F1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(mouse_F2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_F2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(mouse_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(mouse_M2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


mouse_F1 <- subset(mouse_F1, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)
mouse_F2 <- subset(mouse_F2, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)
mouse_M1 <- subset(mouse_M1, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)
mouse_M2 <- subset(mouse_M2, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)



# Visualize QC metrics as a violin plot
VlnPlot(mouse_F1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(mouse_F2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(mouse_M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(mouse_M2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)

plot1 <- FeatureScatter(mouse_F1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_F1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(mouse_F2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_F2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(mouse_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(mouse_M2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


mouse_F1$Species <- 'mouse'
mouse_F1$Sex <- 'F'
mouse_F2$Species <- 'mouse'
mouse_F2$Sex <- 'F'
mouse_M1$Species <- 'mouse'
mouse_M1$Sex <- 'M'
mouse_M2$Species <- 'mouse'
mouse_M2$Sex <- 'M'



#Rのバージョンの確認（Seuratとかはver4で作られている）
sessionInfo()

# Load the dataset (cellrangerのraw dataを読み込み、UMIの行列に変換genes,matrix,barcodesの３点セットを集めておく)
hamster2.data <- Read10X(data.dir = "~/CR_Ham_MPN0311/filtered_feature_bc_matrix_ham")
hamster1.data <- Read10X(data.dir = "~/CR_Ham_MPN0118/outs/filtered_feature_bc_matrix")

hamster_M1<- CreateSeuratObject(counts = hamster1.data, project = "hamsterM1", min.cells = 3, min.features = 250)
hamster_M1

hamster_M2<- CreateSeuratObject(counts = hamster2.data, project = "hamsterM2", min.cells = 3, min.features = 250)
hamster_M2

hamster_mouse_orthologs <- read.csv('~/hamster_mouse_orthologs_v8.csv')


# batchの情報をオブジェクトに入れる
hamster_M1$Species <- 'hamster'
hamster_M2$Species <- 'hamster'
hamster_M1$Sex <- 'M'
hamster_M2$Sex <- 'M'


hamster_M1[["percent.mt"]] <- PercentageFeatureSet(hamster_M1, pattern = "ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|CYTB|ATP8|ATP6")
hamster_M2[["percent.mt"]] <- PercentageFeatureSet(hamster_M2, pattern = "ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|CYTB|ATP8|ATP6")

hamster_M1[["percent.rp"]] <- PercentageFeatureSet(hamster_M1, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
hamster_M2[["percent.rp"]] <- PercentageFeatureSet(hamster_M2, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")


VlnPlot(hamster_M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(hamster_M2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)


# RNAのある細胞を除去
# この時代は核ではなく細胞体seqだったので、ミトコンドリア遺伝子の多い低品質
plot1 <- FeatureScatter(hamster_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hamster_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(hamster_M2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hamster_M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#nFeatureが多すぎる細胞、5％以上ミトコンドリアをcut
hamster_M1 <- subset(hamster_M1, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)
hamster_M2 <- subset(hamster_M2, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)

VlnPlot(hamster_M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(hamster_M2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)


# RNAのある細胞を除去
# この時代は核ではなく細胞体seqだったので、ミトコンドリア遺伝子の多い低品質
plot1 <- FeatureScatter(hamster_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hamster_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(hamster_M2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hamster_M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2





### 2. Get matrix data and identify common genes
# get the mouse counts matrix, keep genes that have human orthologs
hamster_M1_matrix <- GetAssayData(hamster_M1, assay = "RNA", layer = "counts")
hamster_M2_matrix <- GetAssayData(hamster_M2, assay = "RNA", layer = "counts")

hamster_M1_matrix <- hamster_M1_matrix[rownames(hamster_M1_matrix) %in% 
                                         hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol,]
hamster_M2_matrix <- hamster_M2_matrix[rownames(hamster_M2_matrix) %in% 
                                         hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol,]

matched.genes1 <- match(rownames(hamster_M1_matrix), 
                        hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol)
converted_genes1 <- hamster_mouse_orthologs[matched.genes1,"Mus_musculus_gene_symbol"]
rownames(hamster_M1_matrix) <- converted_genes1

matched.genes2 <- match(rownames(hamster_M2_matrix), 
                        hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol)
converted_genes2 <- hamster_mouse_orthologs[matched.genes2,"Mus_musculus_gene_symbol"]
rownames(hamster_M2_matrix) <- converted_genes2


# get the mouse counts matrix
mouse_F1_matrix <- GetAssayData(mouse_F1, assay = "RNA", layer = "counts")
mouse_F2_matrix <- GetAssayData(mouse_F2, assay = "RNA", layer = "counts")
mouse_M1_matrix <- GetAssayData(mouse_M1, assay = "RNA", layer = "counts")
mouse_M2_matrix <- GetAssayData(mouse_M2, assay = "RNA", layer = "counts")

mouse_F1_matrix <- mouse_F1_matrix[rownames(mouse_F1_matrix) %in% 
                                     hamster_mouse_orthologs$Mus_musculus_gene_symbol,]


mouse_F2_matrix <- mouse_F2_matrix[rownames(mouse_F2_matrix) %in% 
                                     hamster_mouse_orthologs$Mus_musculus_gene_symbol,]

mouse_M1_matrix <- mouse_M1_matrix[rownames(mouse_M1_matrix) %in% 
                                     hamster_mouse_orthologs$Mus_musculus_gene_symbol,]


mouse_M2_matrix <- mouse_M2_matrix[rownames(mouse_M2_matrix) %in% 
                                     hamster_mouse_orthologs$Mus_musculus_gene_symbol,]




# # Identify common genes between both datasets
# common.genes1 <- intersect(rownames(mouse_F1_matrix), rownames(mouse_F2_matrix))
# common.genes2 <- intersect(rownames(mouse_M1_matrix), rownames(mouse_M2_matrix))
# common.genes3 <- intersect(rownames(hamster_M1_matrix), rownames(hamster_M2_matrix))
# common.genes <- intersect(common.genes1,common.genes2)
# common.genes <- intersect(common.genes,common.genes3)

# # only keep genes that are common in the mouse and human datasets
# hamster_M1_matrix <- hamster_M1_matrix[common.genes,]
# hamster_M2_matrix <- hamster_M2_matrix[common.genes,]
# mouse_F1_matrix <- mouse_F1_matrix[common.genes,]
# mouse_F2_matrix <- mouse_F2_matrix[common.genes,]
# mouse_M1_matrix <- mouse_M1_matrix[common.genes,]
# mouse_M2_matrix <- mouse_M2_matrix[common.genes,]
# 


### 4. Create merged dataset
# Create new mouse seurat obj
mouse_F1_renamed <- CreateSeuratObject(counts= mouse_F1_matrix, project = "mouseF1")
mouse_F1_renamed$Species <- 'mouse'
mouse_F1_renamed$Sex <- 'F'
mouse_F2_renamed <- CreateSeuratObject(counts= mouse_F2_matrix, project = "mouseF2")
mouse_F2_renamed$Species <- 'mouse'
mouse_F2_renamed$Sex <- 'F'
mouse_M1_renamed <- CreateSeuratObject(counts= mouse_M1_matrix, project = "mouseM1")
mouse_M1_renamed$Species <- 'mouse'
mouse_M1_renamed$Sex <- 'M'
mouse_M2_renamed <- CreateSeuratObject(counts= mouse_M2_matrix, project = "mouseM2")
mouse_M2_renamed$Species <- 'mouse'
mouse_M2_renamed$Sex <- 'M'
hamster_M1_renamed <- CreateSeuratObject(counts= hamster_M1_matrix, project = "hamsterM1")
hamster_M1_renamed$Species <- 'hamster'
hamster_M1_renamed$Sex <- 'M'
hamster_M2_renamed <- CreateSeuratObject(counts= hamster_M2_matrix, project = "hamsterM2")
hamster_M2_renamed$Species <- 'hamster'
hamster_M2_renamed$Sex <- 'M'



# renamedクラスターでミトコンドリア遺伝子名やリボソーム関連タンパクの割合を出す
hamster_M1_renamed[["percent.mt"]] <- PercentageFeatureSet(hamster_M1_renamed, pattern = "^mt-")
hamster_M2_renamed[["percent.mt"]] <- PercentageFeatureSet(hamster_M2_renamed, pattern = "^mt-")

hamster_M1_renamed[["percent.rp"]] <- PercentageFeatureSet(hamster_M1_renamed, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
hamster_M2_renamed[["percent.rp"]] <- PercentageFeatureSet(hamster_M2_renamed, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")

mouse_F1_renamed[["percent.mt"]] <- PercentageFeatureSet(mouse_F1_renamed, pattern = "^mt-")
mouse_F2_renamed[["percent.mt"]] <- PercentageFeatureSet(mouse_F2_renamed, pattern = "^mt-")
mouse_M1_renamed[["percent.mt"]] <- PercentageFeatureSet(mouse_M1_renamed, pattern = "^mt-")
mouse_M2_renamed[["percent.mt"]] <- PercentageFeatureSet(mouse_M2_renamed, pattern = "^mt-")

mouse_F1_renamed[["percent.rp"]] <- PercentageFeatureSet(mouse_F1_renamed, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
mouse_F2_renamed[["percent.rp"]] <- PercentageFeatureSet(mouse_F2_renamed, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
mouse_M1_renamed[["percent.rp"]] <- PercentageFeatureSet(mouse_M1_renamed, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
mouse_M2_renamed[["percent.rp"]] <- PercentageFeatureSet(mouse_M2_renamed, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")

VlnPlot(hamster_M2_renamed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)


# merge:
MPN_merged <- merge(mouse_F1_renamed, y = c(mouse_F2_renamed,mouse_M1_renamed,mouse_M2_renamed,hamster_M1_renamed,hamster_M2_renamed), add.cell.ids = c("mF1","mF2","mM1","mM2","hM1","hM2"))




#version5のSeurat Objectにする
MPN_merged <- UpdateSeuratObject(MPN_merged)

# SCTransform以下、いつもの処理
MPN_merged <- NormalizeData(MPN_merged)
MPN_merged <- FindVariableFeatures(MPN_merged)
MPN_merged <- ScaleData(MPN_merged)
MPN_merged <- RunPCA(MPN_merged,npcs=150)
ElbowPlot(MPN_merged, ndims = 150)
MPN_merged <- IntegrateLayers(object = MPN_merged, method = CCAIntegration, normalization.method = "LogNormalize",
                              verbose = FALSE)

MPN_merged <- FindNeighbors(MPN_merged, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)
# MPN_merged <- FindClusters(MPN_merged, resolution = 0.05)
MPN_merged <- FindClusters(MPN_merged, resolution = 0.4)
MPN_merged <- RunUMAP(MPN_merged, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)

# labeling of clusters
MPN_merged <- RenameIdents(MPN_merged,'0'='Oligodendrocyte','1'='Astrocyte','2'='Neuron','3'='Neuron','4'='Neuron','5'='Neuron','6'='Neuron','7'='Neuron','8'='Neuron','9'='Mixed','10'='Neuron','11'='Neuron','12'='Immature Oligodendrocyte','13'='Neuron','14'='Neuron','15'='Microglia','16'='Neuron','17'='Neuron','18'='Fibroblast','19'='Mixed','20'='Mixed','21'='Newly Formed Oligodendrocyte',"22"="Mixed","23"="Neuron","24"="Neuron")

MPN_merged[["Int_cellltype"]] <- Idents(MPN_merged)


Idents(MPN_merged) <- factor(
  Idents(MPN_merged),
  levels = c(
    "Neuron",
    "Oligodendrocyte",
    "Immature Oligodendrocyte",
    "Newly Formed Oligodendrocyte",
    "Microglia",
    "Astrocyte",
    "Fibroblast",
    "Mixed"
  )
)

MPN_merged$Species <- factor(
  MPN_merged$Species,
  levels = c("mouse", "hamster")
)


# Fig1c
DimPlot(MPN_merged, reduction = "umap", label=TRUE, split.by = "Species",label.size=3)
DimPlot(MPN_merged, group.by='Species', label=TRUE,label.size=3,pt.size = 0.00001,cols = c(
  "mouse" = "#00BFC4",
  "hamster" = "#F8766D"
))

# Fig1d
VlnPlot(
  MPN_merged,
  features = c("Camk2a", "Rbfox3","Mog", "Pdgfra", "Bcas1", "P2ry12", "Aqp4","Dcn","Col3a1"), 
  split.by  = "Species",
  cols = c(
    "mouse" = "#00BFC4",
    "hamster" = "#F8766D"
  ), stack = T,flip = T,group.by="Int_cellltype"
)



# Fig1e
library(tibble)
library(dplyr)
library(ggplot2)
library(scales)
df <- table(
  MPN_merged@meta.data$Int_cellltype,
  MPN_merged@meta.data$Species
)
df <- as.data.frame(df)
df_long <- df %>%
  rownames_to_column("Int_cellltype")
df_long <- df_long %>%
  group_by(Var2) %>%
  mutate(Fraction = Freq / sum(Freq))
levels_order <- c(
  "Neuron",
  "Oligodendrocyte",
  "Immature Oligodendrocyte",
  "Newly Formed Oligodendrocyte",
  "Microglia",
  "Astrocyte",
  "Fibroblast",
  "Mixed"
)

df_long$Var1 <- factor(
  df_long$Var1,
  levels = levels_order
)
ggplot(
  df_long,
  aes(
    x = Var2,
    y = Fraction,
    fill = Var1
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Species",
    y = "fraction"
  ) +
  theme_classic()

# Fig1i
total_marker <- FindMarkers(MPN_merged, assay = 'RNA', slot = 'data',fc.slot = "data", ident.1 = "hamster" ,ident.2 = "mouse", group.by = "Species", min.pct = 0.25, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)
genes_to_label <- c("Atg2a", "Atg5", "Atg7", "Atg10", "Atg13", "Atg14")
keyvals <- ifelse(
  total_marker$avg_log2FC > 1 & total_marker$p_val_adj < 0.05,
  "#F8766D",      # Up（赤）
  ifelse(
    total_marker$avg_log2FC < -1 & total_marker$p_val_adj < 0.05,
    "#00BFC4",    # Down（青）
    "grey80"      # その他
  )
)

names(keyvals) <- ifelse(
  total_marker$avg_log2FC > 1 & total_marker$p_val_adj < 0.05,
  "Up",
  ifelse(
    total_marker$avg_log2FC < -1 & total_marker$p_val_adj < 0.05,
    "Down",
    "NS"
  )
)

EnhancedVolcano(
  total_marker,
  lab = rownames(total_marker),
  x = "avg_log2FC",
  y = "p_val_adj",
  xlab = bquote(~Log[2]~ "Fold change"),
  ylab = bquote(~-Log[10]~ "adjusted p value"),
  selectLab = genes_to_label,
  pointSize = 1,
  labSize = 6.0,
  labCol = "black",
  labFace = "bold",
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 0.8,
  colConnectors = "black",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  colAlpha = 1,
  legendLabSize = 12,
  legendIconSize = 4.0,
  legendPosition = "top",
  colCustom = keyvals
)


#Fig1j
total_marker_up <- total_marker[total_marker$p_val_adj < 0.05 & total_marker$avg_log2FC> 1,]

DEG.gene_up <- rownames(total_marker_up)
DEG.gene_up

ego <- enrichGO(gene = DEG.gene_up,
                ont ="BP",
                OrgDb='org.Mm.eg.db', 
                keyType = "SYMBOL",
                readable      = TRUE)
dotplot(ego)

ego.simple <- simplify(ego)
dotplot(ego.simple)


total_marker_down <- total_marker[total_marker$p_val_adj < 0.05 & total_marker$avg_log2FC< -1,]

DEG.gene_down <- rownames(total_marker_down)
DEG.gene_down

ego <- enrichGO(gene = DEG.gene_down,
                ont ="BP",
                OrgDb='org.Mm.eg.db', 
                keyType = "SYMBOL",
                readable      = TRUE)
dotplot(ego)

ego.simple <- simplify(ego)
dotplot(ego.simple)


### Neuron
MPN_merged.neuron <- subset(MPN_merged, idents = c("Neuron"))
DimPlot(MPN_merged.neuron)
MPN_merged.neuron <- NormalizeData(MPN_merged.neuron)
MPN_merged.neuron <- FindVariableFeatures(MPN_merged.neuron)
MPN_merged.neuron <- ScaleData(MPN_merged.neuron)
MPN_merged.neuron <- RunPCA(MPN_merged.neuron,npcs = 150)
ElbowPlot(MPN_merged.neuron, ndims = 150)
MPN_merged.neuron <- IntegrateLayers(object = MPN_merged.neuron, method = CCAIntegration, normalization.method = "LogNormalize",
                                     verbose = FALSE)
MPN_merged.neuron <- FindNeighbors(MPN_merged.neuron, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)
MPN_merged.neuron <- FindClusters(MPN_merged.neuron, resolution = 0.4)
MPN_merged.neuron <- RunUMAP(MPN_merged.neuron, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)

# labeling of clusters
MPN_merged.neuron <- RenameIdents(MPN_merged.neuron,'0'='In','1'='In','2'='In','3'='In','4'='In','5'='Mixed','6'='In','7'='Ex','8'='In','9'='Ex','10'='Ex','11'='Ex','12'='In','13'='In','14'='In','15'='In','16'='Ex','17'='In','18'='Ex','19'='In','20'='In','21'='In','22'='In','23'='Ex','24'='Ex','25'='Ex','26'='Ex','27'='In')
MPN_merged.neuron[["NeuronType"]] <- Idents(MPN_merged.neuron)
Idents(MPN_merged.neuron) <- factor(
  Idents(MPN_merged.neuron),
  levels = c(
    "Ex",
    "In",
    "Mixed"
  )
)

MPN_merged.neuron$Species <- factor(
  MPN_merged.neuron$Species,
  levels = c("mouse", "hamster")
)


# Fig1f
DimPlot(MPN_merged.neuron, reduction = "umap", split.by = "Species", label=TRUE,cols = c(
  "Ex" = "#61D04F",
  "In" = "#2297E6",
  "Mixed" = "gray62" 
))
DimPlot(MPN_merged.neuron, group.by='Species',cols = c(
  "mouse" = "#00BFC4",
  "hamster" = "#F8766D"
), label=TRUE)

# Fig1g
VlnPlot(
  MPN_merged.neuron,
  features = c("Camk2a", "Rbfox3","Slc17a6", "Gad2", "Mog"), 
  split.by  = "Species",
  cols = c(
    "mouse" = "#00BFC4",
    "hamster" = "#F8766D"
  ), stack = T,flip = T
)


# Fig1h
library(tibble)
library(dplyr)
library(ggplot2)
library(scales)

df <- table(
  MPN_merged.neuron@meta.data$NeuronType,
  MPN_merged.neuron@meta.data$Species
)

df <- as.data.frame(df)

df_long <- df %>%
  rownames_to_column("NeuronType")
df_long <- df_long %>%
  group_by(Var2) %>%
  mutate(Fraction = Freq / sum(Freq))

levels_order <- c(
  "Ex",
  "In",
  "Mixed"
)

df_long$Var1 <- factor(
  df_long$Var1,
  levels = levels_order
)

cols <- c(
  "Ex" = "#61D04F",
  "In" = "#2297E6",
  "Mixed" = "gray62"
)

ggplot(
  df_long,
  aes(
    x = Var2,
    y = Fraction,
    fill = Var1
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Species",
    y = "fraction"
  ) +
  theme_classic()+
  scale_fill_manual(
    values = cols
  ) 



### 興奮性神経細胞のクラスター
MPN_merged.neuron_vGluT2 <- subset(MPN_merged.neuron, idents = c("Ex"))
DimPlot(MPN_merged.neuron_vGluT2)
MPN_merged.neuron_vGluT2 <- NormalizeData(MPN_merged.neuron_vGluT2)
MPN_merged.neuron_vGluT2 <- FindVariableFeatures(MPN_merged.neuron_vGluT2)
MPN_merged.neuron_vGluT2 <- ScaleData(MPN_merged.neuron_vGluT2)
MPN_merged.neuron_vGluT2 <- RunPCA(MPN_merged.neuron_vGluT2,npcs = 150)
ElbowPlot(MPN_merged.neuron_vGluT2, ndims = 150)
MPN_merged.neuron_vGluT2 <- IntegrateLayers(object = MPN_merged.neuron_vGluT2, method = CCAIntegration, normalization.method = "LogNormalize",
                                            verbose = FALSE)

MPN_merged.neuron_vGluT2 <- FindNeighbors(MPN_merged.neuron_vGluT2, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)
MPN_merged.neuron_vGluT2 <- FindClusters(MPN_merged.neuron_vGluT2, resolution =　5.5
)
MPN_merged.neuron_vGluT2 <- RunUMAP(MPN_merged.neuron_vGluT2, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)


# QC→exclude#47 due to Mog expression
VlnPlot(
  MPN_merged.neuron_vGluT2_QC,
  features = c("Camk2a", "Rbfox3","Mog", "Pdgfra", "Bcas1", "P2ry12", "Aqp4","Dcn","Col3a1","Gad1","Gad2"), 
  split.by  = "Species", 
  cols = c(
    "mouse" = "#00BFC4",
    "hamster" = "#F8766D"
  ), stack = T,flip = T
)
VlnPlot(MPN_merged.neuron_vGluT2, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp"), ncol = 2)



MPN_merged.neuron_vGluT2_QC <- subset(MPN_merged.neuron_vGluT2, idents = 47, invert= T)
# manual integration based on their marker genes and UMAP distribution
MPN_merged.neuron_vGluT2_QC <- RenameIdents(MPN_merged.neuron_vGluT2_QC,"29" = "7", "15"="9", "46"="43", "44"="34","3"="0")
MPN_merged.neuron_vGluT2_QC$Ex_celltype <- Idents(MPN_merged.neuron_vGluT2_QC)
cluster_size <- sort(
  table(MPN_merged.neuron_vGluT2_QC@active.ident),
  decreasing = TRUE
)

cluster_size
#  order the culsters accroding to the number of neurons
Idents(MPN_merged.neuron_vGluT2_QC) <- factor(
  Idents(MPN_merged.neuron_vGluT2_QC),
  levels = c(0,9,7,1,2,4,5,6,8,10,11,12,13,14,16,17,18,19,20,21,34,22,23,24,25,26,27,28,30,31,32,33,43,35,36,37,38,39,40,41,42,45,48,49)
)

MPN_merged.neuron_vGluT2_QC$Ex_celltype <-Idents(MPN_merged.neuron_vGluT2_QC)
current_levels <- levels(factor(MPN_merged.neuron_vGluT2_QC$Ex_celltype))
new_labels <- paste0("e", seq_along(current_levels))
MPN_merged.neuron_vGluT2_QC$celltype_display <- factor(
  MPN_merged.neuron_vGluT2_QC$Ex_celltype,
  levels = current_levels,
  labels = new_labels
)

# Fig2a
DimPlot(MPN_merged.neuron_vGluT2_QC, group.by='Species',cols = c(
  "mouse" = "#00BFC4",
  "hamster" = "#F8766D"
), label=TRUE)

DimPlot(MPN_merged.neuron_vGluT2_QC,group.by = "celltype_display", split.by ='Species', label=TRUE)

# Fig2b
markers.to.plot.Ex<-c("Slit3","Meis2","Fign","Otof","Pmfbp1","Col2a1","Ebf2","Tcf7l2","Opn5","Pou6f2","Bcl11a","Kcng2","Fbxl7"
                      ,"Kcnab1","Thrb","Prkd1","Grik3","Nfib","Spon1","Chn2","Synpo2","Sst"
                      ,"Adarb2","Esr1","Ttn","Prlr","Ndst4","Nms","Meis1","Six3","Trhde"
                      ,"Ptprt","Penk","Cacna2d3","Adamtsl1","Dach1","Nos1","Onecut2","Aifm3","Ghr"
                      ,"Cdh12","Reln","Ucn3","Trhr")
DotPlot(MPN_merged.neuron_vGluT2_QC, features = markers.to.plot.Ex,dot.scale = 2,group.by = "celltype_display")+RotatedAxis()+
  coord_flip()

Idents(MPN_merged.neuron_vGluT2_QC)<-"celltype_display"

# Fig2c
library(dplyr)
library(ggplot2)
library(scales)

meta <- MPN_merged.neuron_vGluT2_QC@meta.data %>%
  mutate(
    Species = factor(Species, levels = c("mouse", "hamster"))
  )

overall_hamster_fraction <- mean(meta$Species == "hamster", na.rm = TRUE)

plot_df <- meta %>%
  count(celltype_display, Species, name = "n") %>%
  group_by(celltype_display) %>%
  mutate(
    cluster_total = sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    total_cells = nrow(meta),
    species_fraction = n / total_cells
  )

# 各棒に入れる点線の位置
# = そのクラスター全体の高さ × 全体のmouse割合
line_df <- plot_df %>%
  distinct(celltype_display, cluster_total, total_cells) %>%
  mutate(
    cluster_fraction = cluster_total / total_cells,
    expected_hamster_boundary = cluster_fraction * overall_hamster_fraction
  )

cols <- c(
  "mouse"   = "#00BFC4",
  "hamster" = "#F8766D"
)

ggplot(
  plot_df,
  aes(
    x = celltype_display,
    y = species_fraction,
    fill = Species
  )
) +
  scale_fill_manual(values = cols) +  
  geom_col(width = 0.8) +
  geom_segment(
    data = line_df,
    aes(
      x = as.numeric(celltype_display) - 0.4,
      xend = as.numeric(celltype_display) + 0.4,
      y = expected_hamster_boundary,
      yend = expected_hamster_boundary
    ),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Cell type",
    y = "Fraction of all nuclei"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

library(dplyr)

meta <- MPN_merged.neuron_vGluT2_QC@meta.data %>%
  mutate(
    celltype_display = as.character(celltype_display),
    Species = factor(Species, levels = c("mouse", "hamster"))
  )

celltypes <- sort(unique(meta$celltype_display))

fisher_res <- lapply(celltypes, function(ct) {
  
  tab <- table(
    in_celltype = meta$celltype_display == ct,
    Species = meta$Species
  )
  
  ft <- fisher.test(tab)
  
  data.frame(
    celltype_display = ct,
    Mouse_in_celltype = tab["TRUE", "mouse"],
    Hamster_in_celltype = tab["TRUE", "hamster"],
    Mouse_other = tab["FALSE", "mouse"],
    Hamster_other = tab["FALSE", "hamster"],
    odds_ratio_Mouse_vs_Hamster = unname(ft$estimate),
    p_value = ft$p.value,
    enriched_species = ifelse(
      ft$estimate > 1,
      "mouse",
      "hamster"
    )
  )
}) %>%
  bind_rows() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )
fisher_res <- fisher_res %>%
  arrange(as.numeric(sub("^e", "", celltype_display)))

fisher_res

library(dplyr)
library(ggplot2)
library(stringr)
library(stringr)

cluster_order <- fisher_res %>%
  distinct(celltype_display) %>%
  mutate(
    cluster_num = as.numeric(str_extract(celltype_display, "\\d+$"))
  ) %>%
  arrange(cluster_num) %>%
  pull(celltype_display)
heat_p <- fisher_res %>%
  mutate(
    p_adj_plot = pmax(p_adj_BH, .Machine$double.xmin),
    neg_log10_padj = -log10(p_adj_plot),
    
    signed_score = ifelse(
      enriched_species == "mouse",
      neg_log10_padj,
      -neg_log10_padj
    ),
    
    celltype_display = factor(
      celltype_display,
      levels = cluster_order
    )
  )

ggplot(
  heat_p,
  aes(
    x = celltype_display,
    y = "Species enrichment",
    fill = signed_score
  )
) +
  geom_tile(color = "white") +
  
  geom_text(
    aes(label = sig),
    size = 5
  ) +
  
  scale_fill_gradient2(
    low = "#00BFC4",
    mid = "white",
    high = "#F8766D",
    midpoint = 0,
    name = expression(
      paste(
        "Hamster  <  ",
        -log[10],
        "(FDR)",
        "  >  Mouse"
      )
    )
  ) +
  
  labs(
    x = "Cell type",
    y = NULL
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.ticks.y = element_blank()
  )




















# Fig3a
VlnPlot(MPN_merged.neuron_vGluT2_QC, features = c('Adcyap1','Bdnf',"Lepr","Esr1",'Ptger3','Opn5','Brs3',"Samd3"),stack = T,flip = T,split.by="Species"
        ,idents = c("e2","e4","e9","e10","e12","e13","e14","e15","e18","e19","e24","e25","e26","e27","e28","e29","e32","e34","e35","e36"),cols = c(
          "mouse" = "#00BFC4",
          "hamster" = "#F8766D"
        ))
# Fig3a
library(Seurat)
library(dplyr)

##############################
## 対象クラスター
##############################

target_clusters <- c(
  "e2","e4","e9","e10","e12","e13","e14",
  "e15","e18","e19","e24","e25","e26","e27","e28","e29","e32",
  "e34","e35","e36"
)

##############################
## 対象遺伝子
##############################

target_genes <- c(
  "Adcyap1",
  "Bdnf",
  "Lepr",
  "Esr1",
  "Ptger3",
  "Opn5",
  "Brs3",
  "Samd3"
)

##############################
## SpeciesをIdentに設定
##############################

Idents(MPN_merged.neuron_vGluT2_QC) <- "Species"

##############################
## 結果保存
##############################

result_list <- list()

##############################
## クラスターごとに解析
##############################

for(cl in target_clusters){
  
  message("Cluster ", cl)
  
  obj_sub <- subset(
    MPN_merged.neuron_vGluT2_QC,
    subset = celltype_display == cl
  )
  
  ## mouse・hamster両方存在するか確認
  if(!all(c("mouse","hamster") %in% unique(obj_sub$Species))){
    message("  Skip: one species is missing")
    next
  }
  
  res <- FindMarkers(
    object = obj_sub,
    ident.1 = "mouse",
    ident.2 = "hamster",
    test.use = "wilcox",
    features = target_genes,
    logfc.threshold = 0,
    min.pct = 0
  )
  
  res$Gene <- rownames(res)
  res$Cluster <- cl
  
  result_list[[cl]] <- res
}

##############################
## 結果結合
##############################

result_df <- bind_rows(result_list)

##############################
## 見やすく整理
##############################

result_df <- result_df %>%
  select(
    Cluster,
    Gene,
    avg_log2FC,
    p_val,
    pct.1,
    pct.2,
    everything()
  )

##############################
## 表示
##############################

print(result_df)

# Fig3b
vGluT2_marker <- FindMarkers(MPN_merged.neuron_vGluT2_QC, assay = 'RNA', slot = 'data',fc.slot = "data", ident.1 = "hamster" ,ident.2 = "mouse", group.by = "Species", min.pct = 0.25, logfc.threshold = 0,subset.ident = "e9",test.use = "wilcox",pseudocount.use = 1)
# ポイントサイズを調整
# point_sizes <- rep(1, nrow(total_marker))
# names(point_sizes) <- rownames(total_marker)
# point_sizes[c("Atg2a","Atg5","Atg7","Atg10","Atg13","Atg14")] <- 6

genes_to_label <-  c("Ptger3","Lepr")

keyvals <- ifelse(
  vGluT2_marker$avg_log2FC > 1 & vGluT2_marker$p_val_adj < 0.05,
  "#F8766D",      # Up（赤）
  ifelse(
    vGluT2_marker$avg_log2FC < -1 & vGluT2_marker$p_val_adj < 0.05,
    "#00BFC4",    # Down（青）
    "grey80"      # その他
  )
)

names(keyvals) <- ifelse(
  vGluT2_marker$avg_log2FC > 1 & vGluT2_marker$p_val_adj < 0.05,
  "Up",
  ifelse(
    vGluT2_marker$avg_log2FC < -1 & vGluT2_marker$p_val_adj < 0.05,
    "Down",
    "NS"
  )
)

EnhancedVolcano(
  vGluT2_marker,
  lab = rownames(vGluT2_marker),
  x = "avg_log2FC",
  y = "p_val_adj",
  xlab = bquote(~Log[2]~ "Fold change"),
  ylab = bquote(~-Log[10]~ "adjusted p value"),
  selectLab = genes_to_label,
  pointSize = 1,
  labSize = 6.0,
  labCol = "black",
  labFace = "bold",
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 0.8,
  colConnectors = "black",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  colAlpha = 1,
  legendLabSize = 12,
  legendIconSize = 4.0,
  legendPosition = "top",
  colCustom = keyvals
)


# GO analysis
vGluT2_marker_up <- vGluT2_marker[vGluT2_marker$p_val_adj < 0.05 & vGluT2_marker$avg_log2FC> 1,]

DEG.gene_up <- rownames(vGluT2_marker_up)
DEG.gene_up

ego <- enrichGO(gene = DEG.gene_up,
                ont ="BP",
                OrgDb='org.Mm.eg.db', 
                keyType = "SYMBOL",
                readable      = TRUE)
dotplot(ego)

res <- summary(ego)

ego.simple <- simplify(ego)
dotplot(ego.simple)



#Fig3cright
vGluT2_marker_down <- vGluT2_marker[vGluT2_marker$p_val_adj < 0.05 & vGluT2_marker$avg_log2FC< -1,]

DEG.gene_down <- rownames(vGluT2_marker_down)
DEG.gene_down

ego_down <- enrichGO(gene = DEG.gene_down,
                     ont ="BP",
                     OrgDb='org.Mm.eg.db', 
                     keyType = "SYMBOL",
                     readable      = TRUE)
dotplot(ego_down)

res_down <- summary(ego_down)
ego_down.simple <- simplify(ego_down)
dotplot(ego_down.simple)



# inhibitory neuron
MPN_merged.neuron_vGat <- subset(MPN_merged.neuron, idents = c("In"))
DimPlot(MPN_merged.neuron_vGat)

#いつもの処理
# MPN_merged.neuron_vGat[["RNA"]] <- split(MPN_merged.neuron_vGat[["RNA"]], f = MPN_merged.neuron_vGat$orig.ident)
MPN_merged.neuron_vGat <- NormalizeData(MPN_merged.neuron_vGat)
MPN_merged.neuron_vGat <- FindVariableFeatures(MPN_merged.neuron_vGat)
MPN_merged.neuron_vGat <- ScaleData(MPN_merged.neuron_vGat)
MPN_merged.neuron_vGat <- RunPCA(MPN_merged.neuron_vGat,npcs = 150)
ElbowPlot(MPN_merged.neuron_vGat, ndims = 150)
MPN_merged.neuron_vGat <- IntegrateLayers(object = MPN_merged.neuron_vGat, method = CCAIntegration, normalization.method = "LogNormalize",
                                          verbose = FALSE)

MPN_merged.neuron_vGat <- FindNeighbors(MPN_merged.neuron_vGat, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)
MPN_merged.neuron_vGat <- FindClusters(MPN_merged.neuron_vGat, resolution = 4)#元々8.5
MPN_merged.neuron_vGat <- RunUMAP(MPN_merged.neuron_vGat, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)

# Exclude doublet clusters
MPN_merged.neuron_vGat_QC <- subset(MPN_merged.neuron_vGat, idents = c(0,5,12,20,24,42,50,61,54), invert= T)
MPN_merged.neuron_vGat_QC <- RenameIdents(MPN_merged.neuron_vGat_QC,"14"="6","38"="6","33"="29","53"="35","60"="26","49"="26","53"="35","60"="26","49"="26") 

MPN_merged.neuron_vGat_QC$In_celltype <-Idents(MPN_merged.neuron_vGat_QC)

current_levels <- levels(factor(MPN_merged.neuron_vGat_QC$In_celltype))
new_labels <- paste0("i", seq_along(current_levels))
MPN_merged.neuron_vGat_QC$celltype_display <- factor(
  MPN_merged.neuron_vGat_QC$In_celltype,
  levels = current_levels,
  labels = new_labels
)


# Fig2a
DimPlot(MPN_merged.neuron_vGat_QC, group.by='Species',cols = c(
  "mouse" = "#00BFC4",
  "hamster" = "#F8766D"
), label=TRUE)
DimPlot(MPN_merged.neuron_vGat_QC,group.by = "celltype_display", split.by ='Species', label=TRUE)


Idents(MPN_merged.neuron_vGat_QC)<-"celltype_display"
# markers.to.plot.In =  c('Esr1','Myo3a','Fign','Pde3a','Vipr','Trpc7','Onecut2','Gpc3','Marchf3','Cdon',
#                      'Scn5a','Lhx8','Prdm16','Calcr','Moxd1','Tafa1','Pou6f2','Nfib','Atpap1l','Map3k15',
#                      'Pard3b','Meis2','Prox1','Nmu','Tmem132c','Nfia','Sntb1','Six3','Drd2','Ddc',
#                      'Ttn','Kcnmb4','Satb2','Nr4a2','Sema3a','Lhx6','Slc18a2','Ebf1','Thsd4','Rec114',
#                      'Nxph2','Gpc4'
#                      )
markers.to.plot.Ex<-c("Slit3","Meis2","Fign","Otof","Tcf7l2","Col2a1","Opn5","Glra2","Pou6f2","Fbxl7"
                      ,"Kcnab1","Kcng2","Bcl11a","Grik3","Thrb","Chn2","Nfib","Spon1","Vav3","Zbtb16"
                      ,"Adarb2","Esr1","Prlr","Nms","Ttn","Ndst4","Cacna2d3","Six3","Meis1","Trhde"
                      ,"Ptprt","Pmfbp1","Penk","Dach1","Adamtsl1","Prkd1","Nos1","Onecut2","Aifm3","Ghr"
                      ,"Cdh12","Reln","Ebf2","Ucn3","Trhr")
markers.to.plot.In =  c('Esr1',"Sv2b",'Zeb2','Rorb','Cacna2d1',"Trpc4","Satb1",'Lhx8','Tac1',"Trdn"
                        ,'Penk',"Meis2",'Pde11a',"Adgrv1",'Map3k15',"Tshz1","Ghr","Pou6f2","Tcf4",'Marchf3'
                        ,"Ebf1","Drd3","Gli3",'Tmem132c','Pmfbp1',"Lama3",'Cdh20','Nfib',"Onecut2",'Prox1',
                        "Satb2","Lypd6","Ndst4",'Six3',"Trhde","Fign",'Kcnh8','Ptprm','Ptger3',"Otof",
                        'Sema3a','Arhgap6',"Sst",'Reln','Vipr2',"Il1rapl2",'Tafa1',"Crh"
                        
                        
)


DotPlot(MPN_merged.neuron_vGat_QC, features = markers.to.plot.In,dot.scale = 2,group.by = "celltype_display")+RotatedAxis()+
  coord_flip()



# Fig2c

library(dplyr)
library(ggplot2)
library(scales)

meta <- MPN_merged.neuron_vGat_QC@meta.data %>%
  mutate(
    Species = factor(Species, levels = c("mouse", "hamster"))
  )

# 全体におけるmouseの割合
overall_hamster_fraction <- mean(meta$Species == "hamster", na.rm = TRUE)

plot_df <- meta %>%
  count(celltype_display, Species, name = "n") %>%
  group_by(celltype_display) %>%
  mutate(
    cluster_total = sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    total_cells = nrow(meta),
    species_fraction = n / total_cells
  )

# 各棒に入れる点線の位置
# = そのクラスター全体の高さ × 全体のmouse割合
line_df <- plot_df %>%
  distinct(celltype_display, cluster_total, total_cells) %>%
  mutate(
    cluster_fraction = cluster_total / total_cells,
    expected_hamster_boundary = cluster_fraction * overall_hamster_fraction
  )

ggplot(
  plot_df,
  aes(
    x = celltype_display,
    y = species_fraction,
    fill = Species
  )
) +
  geom_col(width = 0.8) +
  geom_segment(
    data = line_df,
    aes(
      x = as.numeric(celltype_display) - 0.4,
      xend = as.numeric(celltype_display) + 0.4,
      y = expected_hamster_boundary,
      yend = expected_hamster_boundary
    ),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Cell type",
    y = "Fraction of all nuclei"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


meta <- MPN_merged.neuron_vGat_QC@meta.data %>%
  mutate(
    celltype_display = as.character(celltype_display),
    Species = factor(Species, levels = c("mouse", "hamster"))
  )

celltypes <- sort(unique(meta$celltype_display))

fisher_res <- lapply(celltypes, function(ct) {
  
  tab <- table(
    in_celltype = meta$celltype_display == ct,
    Species = meta$Species
  )
  
  ft <- fisher.test(tab)
  
  data.frame(
    celltype_display = ct,
    Mouse_in_celltype = tab["TRUE", "mouse"],
    Hamster_in_celltype = tab["TRUE", "hamster"],
    Mouse_other = tab["FALSE", "mouse"],
    Hamster_other = tab["FALSE", "hamster"],
    odds_ratio_Mouse_vs_Hamster = unname(ft$estimate),
    p_value = ft$p.value,
    enriched_species = ifelse(
      ft$estimate > 1,
      "mouse",
      "hamster"
    )
  )
}) %>%
  bind_rows() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )
fisher_res <- fisher_res %>%
  arrange(as.numeric(sub("^e", "", celltype_display)))

fisher_res

library(dplyr)
library(ggplot2)
library(stringr)
library(stringr)

cluster_order <- fisher_res %>%
  distinct(celltype_display) %>%
  mutate(
    cluster_num = as.numeric(str_extract(celltype_display, "\\d+$"))
  ) %>%
  arrange(cluster_num) %>%
  pull(celltype_display)
heat_p <- fisher_res %>%
  mutate(
    p_adj_plot = pmax(p_adj_BH, .Machine$double.xmin),
    neg_log10_padj = -log10(p_adj_plot),
    
    signed_score = ifelse(
      enriched_species == "mouse",
      neg_log10_padj,
      -neg_log10_padj
    ),
    
    celltype_display = factor(
      celltype_display,
      levels = cluster_order
    )
  )

ggplot(
  heat_p,
  aes(
    x = celltype_display,
    y = "Species enrichment",
    fill = signed_score
  )
) +
  geom_tile(color = "white") +
  
  geom_text(
    aes(label = sig),
    size = 5
  ) +
  
  scale_fill_gradient2(
    low = "#00BFC4",
    mid = "white",
    high = "#F8766D",
    midpoint = 0,
    name = expression(
      paste(
        "Hamster  <  ",
        -log[10],
        "(FDR)",
        "  >  Mouse"
      )
    )
  ) +
  
  labs(
    x = "Cell type",
    y = NULL
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.ticks.y = element_blank()
  )s