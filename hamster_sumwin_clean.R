library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(sets)
library(org.Mm.eg.db)
library(rvcheck)
library(clusterProfiler)
library(KEGGREST)
library(DOSE)
library(ReactomePA)
library(stringr)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(tibble)
library(ggplot2)
library(patchwork)
library(enrichplot)
library(GOSemSim)  # if using semantic similarity
library(pathview)


hamster_mouse_orthologs <- read.csv('~/one-to-one_ortholog_genes.csv')

# Load the dataset
summer_M1.data <- Read10X(data.dir = "~/CR_summerHam_MPN0127/outs/filtered_feature_bc_matrix")
winter_M1.data <- Read10X(data.dir = "~/CR_winterHam_MPN1015/outs/filtered_feature_bc_matrix")
winter_M2.data <- Read10X(data.dir = "~/CR_winterHam_MPN1017/outs/filtered_feature_bc_matrix")
winter_F1.data <- Read10X(data.dir = "~/CR_winterHam_MPN0210/outs/filtered_feature_bc_matrix")
summer_F1.data <- Read10X(data.dir = "~/CR_summerHam_MPN0529/outs/filtered_feature_bc_matrix")

summer_M1<- CreateSeuratObject(counts = summer_M1.data, project = "summerM1", min.cells = 3, min.features = 250)
summer_M1
winter_M1<- CreateSeuratObject(counts = winter_M1.data, project = "winterM1", min.cells = 3, min.features = 250)
winter_M1
winter_M2<- CreateSeuratObject(counts = winter_M2.data, project = "winterM2", min.cells = 3, min.features = 250)
winter_M2
winter_F1<- CreateSeuratObject(counts = winter_F1.data, project = "winterF1", min.cells = 3, min.features = 250)
winter_F1
summer_F1<- CreateSeuratObject(counts = summer_F1.data, project = "summerF1", min.cells = 3, min.features = 250)
summer_F1


#QC
summer_M1[["percent.mt"]] <- PercentageFeatureSet(summer_M1,pattern = "ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|CYTB|ATP8|ATP6")
winter_M1[["percent.mt"]] <- PercentageFeatureSet(winter_M1,pattern = "ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|CYTB|ATP8|ATP6")
winter_M2[["percent.mt"]] <- PercentageFeatureSet(winter_M2,pattern = "ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|CYTB|ATP8|ATP6")
winter_F1[["percent.mt"]] <- PercentageFeatureSet(winter_F1,pattern = "ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|CYTB|ATP8|ATP6")
summer_F1[["percent.mt"]] <- PercentageFeatureSet(summer_F1,pattern = "ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|CYTB|ATP8|ATP6")

summer_M1[["percent.rp"]] <- PercentageFeatureSet(summer_M1, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
winter_M1[["percent.rp"]] <- PercentageFeatureSet(winter_M1, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
winter_M2[["percent.rp"]] <- PercentageFeatureSet(winter_M2, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
winter_F1[["percent.rp"]] <- PercentageFeatureSet(winter_F1, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
summer_F1[["percent.rp"]] <- PercentageFeatureSet(summer_F1, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")


# Visualize QC metrics as a violin plot
VlnPlot(summer_M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
VlnPlot(winter_M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
VlnPlot(winter_M2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
VlnPlot(winter_F1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
VlnPlot(summer_F1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)


plot1 <- FeatureScatter(summer_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(summer_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(winter_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(winter_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(winter_M2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(winter_M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(winter_F1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(winter_F1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(summer_F1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(summer_F1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



#細胞の削減（UMI数が400-10000のものに絞っている）
summer_M1 <- subset(summer_M1, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)
winter_M1 <- subset(winter_M1, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)
winter_M2 <- subset(winter_M2, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)
winter_F1 <- subset(winter_F1, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)
summer_F1 <- subset(summer_F1, subset = nCount_RNA > 400 & nCount_RNA < 100000  & percent.mt < 10 & percent.rp < 10)



# Visualize QC metrics as a violin plot
VlnPlot(summer_M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(winter_M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(winter_M2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(winter_F1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
VlnPlot(summer_F1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)


# species, sexスロットの追加

summer_M1$Mode <- 'summer'
summer_M1$Sex <- 'M'
winter_M1$Mode <- 'winter'
winter_M1$Sex <- 'M'
winter_M2$Mode <- 'winter'
winter_M2$Sex <- 'M'
winter_F1$Mode <- 'winter'
winter_F1$Sex <- 'F'
summer_F1$Mode <- 'summer'
summer_F1$Sex <- 'F'


# merge
sumwinHam_merged <- merge(summer_M1, y = c(winter_M1,winter_M2,winter_F1,summer_F1), add.cell.ids = c("sM1","wM2","wM1","wF1","sF1"))

# SCTransform以下、いつもの処理
sumwinHam_merged <- SCTransform(sumwinHam_merged, vars.to.regress = c("percent.mt","percent.rp"))
sumwinHam_merged <- RunPCA(sumwinHam_merged,npcs=150)
ElbowPlot(sumwinHam_merged, ndims = 150)
sumwinHam_merged <- IntegrateLayers(object = sumwinHam_merged, method = CCAIntegration, normalization.method = "LogNormalize",
                                    verbose = FALSE)
sumwinHam_merged <- FindNeighbors(sumwinHam_merged, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)
sumwinHam_merged <- FindClusters(sumwinHam_merged, resolution = 0.3)
sumwinHam_merged <- RunUMAP(sumwinHam_merged, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)

sumwinHam_merged <- RenameIdents(sumwinHam_merged,'0'='Oligodendrocyte','1'='Astrocyte','2'='Neuron','3'='Neuron','4'='Neuron','5'='MG','6'='Neuron','7'='NFO','8'='Neuron','9'='Neuron','10'='Neuron','11'='Neuron','12'='Neuron','13'='Neuron','14'='Neuron','15'='Neuron','16'='Neuron','17'='Neuron','18'='Mixed','19'='Neuron','20'='Neuron','21'='Endothelial','22'='IO','23'='MG','24'='Mixed','25'='Neuron','26'='Mixed','27'='Neuron','28'='Mixed','29'='Fibroblast')
sumwinHam_merged[["Cell_class"]] <- Idents(sumwinHam_merged)

# Fig4b
DimPlot(sumwinHam_merged, group.by='Cell_class',cols = c(
  "Neuron" = "#F8766D",
  "Oligodendrocyte" = "#D89000",
  "IO"="#A3A500",
  "NFO"="#39B600",
  "Astrocyte"="#00B0F6",
  "MG"="#00BFC4",
  "Endothelial"="#E76BF3",
  "Fibroblast"="#9590FF",
  "Mixed"="gray50"
), label=TRUE,split.by="Mode")
DimPlot(sumwinHam_merged, reduction = "umap", group.by = "Mode",label.size=3,
        cols = c(
          "winter" = "#2166AC",
          "summer" = "#B2182B"
        ))

# Fig4c
VlnPlot(
  sumwinHam_merged,
  features = c("Camk2a", "Rbfox3","Mog", "Pdgfra", "Bcas1", "Csf1r", "Aqp4","Cldn5","Col3a1"), 
  split.by  = "Mode",group.by = "Cell_class",
  cols = c(
    "winter" = "#2166AC",
    "summer" = "#B2182B"
  ), stack = T,flip = T
)

# Fig4d
df<-table(sumwinHam_merged@meta.data$Cell_class,sumwinHam_merged@meta.data$Mode)
df <- as.data.frame(df) 
df_long <- df %>%
  rownames_to_column("CellType") 
df_long <- df_long %>%
  group_by(Var2) %>%
  mutate(Fraction = Freq / sum(Freq))
df_long


ggplot(df_long, aes(x = Var2, y = Fraction, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Fraction of each cell class in each condition",
    x = "Condition",
    y = "Fraction"
  ) +
  theme_minimal()


# Find marker genes in each cluster
sumwinHam_merged <- PrepSCTFindMarkers(sumwinHam_merged)
sumwinHam_merged@active.assay='SCT'

sumwinHam_merged.markers <- FindMarkers(sumwinHam_merged, ident.1 = 22, logfc.threshold = 0.25, only.pos = TRUE)
sumwinHam_merged.markers_top40 <- sumwinHam_merged.markers[sumwinHam_merged.markers$p_val_adj < 0.01 & sumwinHam_merged.markers$pct.1> 0.1,]%>% top_n(n = 20, wt = avg_log2FC)
# DotPlot(merged.obj_ViscMN, features = rownames(merged.obj_ViscMN.markers_top40))

plots <- VlnPlot(sumwinHam_merged, features = rownames(sumwinHam_merged.markers_top40),layer="data", assay='SCT',
                 stack = T,flip = T, pt.size = 0.3, split.by = "Mode")
wrap_plots(plots = plots, ncol = 1)



## extract Neurons
sumwinHam_merged.neuron <- subset(sumwinHam_merged, idents = c("Neuron"))
DimPlot(sumwinHam_merged.neuron)

sumwinHam_merged.neuron[["RNA"]] <- split(sumwinHam_merged.neuron[["RNA"]], f = sumwinHam_merged.neuron$orig.ident)
sumwinHam_merged.neuron <- SCTransform(sumwinHam_merged.neuron, vars.to.regress = c("percent.mt","percent.rp"))
sumwinHam_merged.neuron <- RunPCA(sumwinHam_merged.neuron,npcs = 50)
ElbowPlot(sumwinHam_merged.neuron, ndims = 50)
sumwinHam_merged.neuron <- IntegrateLayers(object = sumwinHam_merged.neuron, method = CCAIntegration, normalization.method = "SCT",
                                           verbose = FALSE)
sumwinHam_merged.neuron <- FindNeighbors(sumwinHam_merged.neuron, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)
sumwinHam_merged.neuron <- FindClusters(sumwinHam_merged.neuron, resolution = 0.2)
sumwinHam_merged.neuron <- RunUMAP(sumwinHam_merged.neuron, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)

DimPlot(sumwinHam_merged.neuron, reduction = "umap", split.by = "Mode", label=TRUE)

sumwinHam_merged.neuron<- RenameIdents(sumwinHam_merged.neuron,
                                       '0'='In','1'='In','2'='In','3'='Ex','4'='In','5'='In','6'='In','7'='Ex',
                                       '8'='In','9'='In','10'='In','11'='Ex','12'='In','13'='Ex','14'='Ex','15'='In',
                                       '16'='Mixed','17'='Ex','18'='Mixed','19'='Mixed','20'='Ex','21'='Ex',
                                       '22'='In','23'='Ex','24'='Mixed'
)

sumwinHam_merged.neuron[["Neuron_Class"]] <- Idents(sumwinHam_merged.neuron)
sumwinHam_merged.neuron$Neuron_Class <- factor(
  sumwinHam_merged.neuron $Neuron_Class,
  levels = c("Ex","In","Mixed")
)

# Fig4e
DimPlot(sumwinHam_merged.neuron, reduction = "umap", split.by = "Mode", label=TRUE,cols = c(
  "Ex" = "#61D04F",
  "In" = "#2297E6",
  "Mixed" = "gray62" 
))
DimPlot(sumwinHam_merged.neuron, reduction = "umap", group.by = "Mode",label.size=3,
        cols = c(
          "winter" = "#2166AC",
          "summer" = "#B2182B"
        ))


# Fig4f
VlnPlot(
  sumwinHam_merged.neuron,
  features = c("Camk2a", "Rbfox3","Slc17a6", "Gad2", "Bcas1","Aqp4","Csf1r"), 
  split.by  = "Mode",group.by = "Neuron_Class",
  cols = c(
    "winter" = "#2166AC",
    "summer" = "#B2182B"
  ), stack = T,flip = T
)

# Fig4g
df<-table(sumwinHam_merged.neuron@meta.data$Neuron_Class,sumwinHam_merged.neuron@meta.data$Mode)
df <- as.data.frame(df) 
df_long <- df %>%
  rownames_to_column("Neuron_Class") 
df_long <- df_long %>%
  group_by(Var2) %>%
  mutate(Fraction = Freq / sum(Freq))
df_long

ggplot(df_long, aes(x = Var2, y = Fraction, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "各Conditionにおける細胞種の割合",
    x = "Condition",
    y = "割合"
  ) +
  theme_minimal()


# Fig4h
sumwinHam_merged.neuron<-PrepSCTFindMarkers(sumwinHam_merged.neuron)

total_marker <- FindMarkers(sumwinHam_merged.neuron, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "winter" ,ident.2 = "summer", group.by = "Mode", min.pct = 0.25, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)
# ポイントサイズを調整
# point_sizes <- rep(1, nrow(total_marker))
# names(point_sizes) <- rownames(total_marker)
# point_sizes[c("Atg2a","Atg5","Atg7","Atg10","Atg13","Atg14")] <- 6


keyvals <- ifelse(
  total_marker$avg_log2FC > 1 & total_marker$p_val_adj < 0.05,
  "#2166AC",      # Up（赤）
  ifelse(
    total_marker$avg_log2FC < -1 & total_marker$p_val_adj < 0.05,
    "#B2182B",    # Down（青）
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
library(EnhancedVolcano)
EnhancedVolcano(
  total_marker,
  lab = rownames(total_marker),
  x = "avg_log2FC",
  y = "p_val_adj",
  xlab = bquote(~Log[2]~ "Fold change"),
  ylab = bquote(~-Log[10]~ "adjusted p value"),
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


# Fig4i
library(clusterProfiler)
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



### 興奮性神経細胞のクラスター
sumwinHam_merged.neuron_vGluT2 <- subset(sumwinHam_merged.neuron, idents = c("Ex"))
DimPlot(sumwinHam_merged.neuron_vGluT2)
sumwinHam_merged.neuron_vGluT2 <- SCTransform(sumwinHam_merged.neuron_vGluT2, vars.to.regress = c("percent.mt","percent.rp"), return.only.var.genes = TRUE)
sumwinHam_merged.neuron_vGluT2 <- RunPCA(sumwinHam_merged.neuron_vGluT2,npcs = 150)
ElbowPlot(sumwinHam_merged.neuron_vGluT2, ndims = 150)
sumwinHam_merged.neuron_vGluT2 <- IntegrateLayers(object = sumwinHam_merged.neuron_vGluT2, method = CCAIntegration, normalization.method = "SCT",
                                                  verbose = FALSE)

sumwinHam_merged.neuron_vGluT2 <- FindNeighbors(sumwinHam_merged.neuron_vGluT2, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)
sumwinHam_merged.neuron_vGluT2 <- FindClusters(sumwinHam_merged.neuron_vGluT2, resolution = 1.8)
sumwinHam_merged.neuron_vGluT2 <- RunUMAP(sumwinHam_merged.neuron_vGluT2, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)

# Fig5a
DimPlot(sumwinHam_merged.neuron_vGluT2, group.by='Species',pt.size = 0.1)
DimPlot(sumwinHam_merged.neuron_vGluT2, reduction = "umap", group.by = "Mode",label.size=3,
        cols = c(
          "winter" = "#2166AC",
          "summer" = "#B2182B"
        ))


# Fig5b
df<-table(sumwinHam_merged.neuron_vGluT2@meta.data$seurat_clusters,sumwinHam_merged.neuron_vGluT2@meta.data$Mode)
df <- as.data.frame(df) 
df_long <- df %>%
  rownames_to_column("seurat_clusters") 
df_long <- df_long %>%
  group_by(Var2) %>%
  mutate(Fraction = Freq / sum(Freq))
df_long


# 割合の棒グラフを描画
ggplot(df_long, aes(x = Var2, y = Fraction, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "各Conditionにおける細胞種の割合",
    x = "Condition",
    y = "割合"
  ) +
  theme_minimal()


# Fig5c
sumwinHam_merged.neuron_vGluT2<-PrepSCTFindMarkers(sumwinHam_merged.neuron_vGluT2)
total_marker <- FindMarkers(sumwinHam_merged.neuron_vGluT2, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "winter" ,ident.2 = "summer", group.by = "Mode", min.pct = 0.25, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)
# ポイントサイズを調整
# point_sizes <- rep(1, nrow(total_marker))
# names(point_sizes) <- rownames(total_marker)
# point_sizes[c("Atg2a","Atg5","Atg7","Atg10","Atg13","Atg14")] <- 6
write.csv2(total_marker,"totalexneuron_DEGs.csv")
# upとdownの数
up_genes <- total_marker %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > log2(2)
  )

down_genes <- total_marker %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC < log2(0.5)
  )

n_up <- nrow(up_genes)
n_down <- nrow(down_genes)

n_up
n_down


keyvals <- ifelse(
  total_marker$avg_log2FC > 1 & total_marker$p_val_adj < 0.05,
  "#2166AC",      # Up（赤）
  ifelse(
    total_marker$avg_log2FC < -1 & total_marker$p_val_adj < 0.05,
    "#B2182B",    # Down（青）
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
  lab = rep("", nrow(total_marker)),
  x = "avg_log2FC",
  y = "p_val_adj",
  xlab = bquote(~Log[2]~ "Fold change"),
  ylab = bquote(~-Log[10]~ "adjusted p value"),
  pointSize = 1,
  pCutoff = 0.05,
  FCcutoff = 1.0,
  colAlpha = 1,
  legendLabSize = 12,
  legendIconSize = 4.0,
  legendPosition = "top",
  colCustom = keyvals
)

# Fig5d
clusters <- levels(Idents(sumwinHam_merged.neuron_vGluT2))
deg_list <- list()

for(cl in clusters){
  
  deg <- FindMarkers(
    sumwinHam_merged.neuron_vGluT2,
    subset.ident = cl,
    group.by = "Mode",
    ident.1 = "winter",
    ident.2 = "summer",
    recorrect_umi = FALSE,
    logfc.threshold = 0,
    min.pct = 0
  )
  
  deg <- deg[
    abs(deg$avg_log2FC) > 1 &
      deg$p_val_adj < 0.05,
  ]
  
  deg$gene <- rownames(deg)
  
  deg_list[[cl]] <- deg
}

deg_direction <- do.call(
  rbind,
  lapply(names(deg_list), function(cl){
    
    tmp <- deg_list[[cl]]
    
    if(nrow(tmp) == 0){
      return(NULL)
    }
    
    tmp$Cluster <- cl
    
    tmp$Direction <- ifelse(
      tmp$avg_log2FC > 0,
      "Winter",
      "Summer"
    )
    
    tmp[, c(
      "gene",
      "Cluster",
      "Direction"
    )]
  })
)

gene_freq <- deg_direction %>%
  distinct(
    gene,
    Cluster,
    Direction
  ) %>%
  count(
    gene,
    Direction,
    name = "n_cluster"
  )

n_cluster <- length(clusters)

shared_cutoff <- ceiling(
  0.50 * n_cluster
)

specific_cutoff <- floor(
  0.10 * n_cluster
)

gene_freq$Category <- case_when(
  gene_freq$n_cluster >= shared_cutoff ~
    "Shared",
  
  gene_freq$n_cluster <= specific_cutoff ~
    "Cell type-specific",
  
  TRUE ~
    "Intermediate"
)

deg_direction <- merge(
  deg_direction,
  gene_freq,
  by = c(
    "gene",
    "Direction"
  )
)

plot_df2 <- deg_direction %>%
  count(
    Cluster,
    Direction,
    Category
  )

plot_df2$N <- plot_df2$n

plot_df2$N[
  plot_df2$Direction == "Summer"
] <-
  -plot_df2$N[
    plot_df2$Direction == "Summer"
  ]

plot_df2$Cluster <- factor(
  plot_df2$Cluster,
  levels = as.character(
    sort(
      as.numeric(
        unique(plot_df2$Cluster)
      )
    )
  )
)

plot_df2$FillGroup <- paste(
  plot_df2$Direction,
  plot_df2$Category,
  sep = "_"
)

ggplot(
  plot_df2,
  aes(
    x = Cluster,
    y = N,
    fill = FillGroup
  )
) +
  geom_col(width = 0.8) +
  
  geom_hline(yintercept = 0) +
  scale_fill_manual(
    values = c(
      
      # Winter系（赤）
      "Winter_Shared" = "#2166AC",
      "Winter_Intermediate" = "#67A9CF",
      "Winter_Cell type-specific" = "#D1E5F0",
      
      # Summer系（青）
      "Summer_Shared" = "#B2182B",
      "Summer_Intermediate" = "#EF8A62",
      "Summer_Cell type-specific" = "#FDDBC7"
    )
  ) +
  theme_classic()


# Fig5e
rownames(total_marker) <- ifelse(row.names(total_marker) %in% hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol, 
                                 hamster_mouse_orthologs$Mus_musculus_gene_symbol[match(row.names(total_marker), hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol)], 
                                 row.names(total_marker))
hamster_mouse_orthologs <- read.csv('~/one-to-one_ortholog_genes.csv')
hamster_mouse_orthologs <- hamster_mouse_orthologs %>% filter(!row_number() %in% c(match('GeneX',hamster_mouse_orthologs$Mus_musculus_gene_symbol)))
total_marker_GOlist_down <- total_marker[total_marker$p_val_adj < 0.05 & total_marker$avg_log2FC< -1,]

DEG.gene_down <- rownames(total_marker_GOlist_down)
ego_down <- enrichGO(gene = DEG.gene_down,
                     ont ="BP",
                     OrgDb='org.Mm.eg.db', 
                     keyType = "SYMBOL",
                     readable      = TRUE)
dotplot(ego_down)+ggtitle('Total_summer')
ego_down.simple <- simplify(ego_down)
dotplot(ego_down.simple, showCategory=10)+ggtitle(32)



#Fig5f

library(Seurat)
library(dplyr)

##############################
## 対象クラスター
##############################

target_clusters <- c(
  0,1,2,3,4,6,7,9,12,14,
  15,17,18,19,20,21,23,24,
  27,28,30,31,33,34,37,43,
  45,47
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
## ModeをIdentに設定
##############################

Idents(sumwinHam_merged.neuron_vGluT2) <- "Mode"

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
    sumwinHam_merged.neuron_vGluT2,
    subset = seurat_clusters == cl
  )
  
  ## summerとwinterの両方が存在するか確認
  if(!all(c("summer","winter") %in% unique(obj_sub$Mode))){
    next
  }
  
  res <- FindMarkers(
    object = obj_sub,
    ident.1 = "winter",
    ident.2 = "summer",
    test.use = "wilcox",
    features = target_genes,
    logfc.threshold = 0,
    min.pct = 0,
    recorrect_umi = F
  )
  
  res$Gene <- rownames(res)
  res$Cluster <- cl
  
  result_list[[as.character(cl)]] <- res
  
}

##############################
## 結果結合
##############################

result_df <- bind_rows(result_list)

##############################
## 必要な列のみ
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


VlnPlot(sumwinHam_merged.neuron_vGluT2,idents = c(0,1,2,3,4,6,7,9,12,14,15,17,18,19,20,21,23,24,27,28,30,31,33,34,37,43,45,47), features = target_genes, 
        pt.size = 0.3, split.by = "Mode",stack = T,flip = T)


# Fig5g,h
total_marker <- FindMarkers(sumwinHam_merged.neuron_vGluT2, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "winter" ,ident.2 = "summer", group.by = "Mode",subset.ident = 4,recorrect_umi=F, min.pct = 0.25, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

# upとdownの数
up_genes <- total_marker %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > log2(2)
  )

down_genes <- total_marker %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC < log2(0.5)
  )

n_up <- nrow(up_genes)
n_down <- nrow(down_genes)

n_up
n_down


keyvals <- ifelse(
  total_marker$avg_log2FC > 1 & total_marker$p_val_adj < 0.05,
  "#2166AC",      # Up（赤）
  ifelse(
    total_marker$avg_log2FC < -1 & total_marker$p_val_adj < 0.05,
    "#B2182B",    # Down（青）
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
  selectLab = "Adcyap1", 
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


rownames(total_marker) <- ifelse(row.names(total_marker) %in% hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol, 
                                 hamster_mouse_orthologs$Mus_musculus_gene_symbol[match(row.names(total_marker), hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol)], 
                                 row.names(total_marker))
# エラーが出たら以下のコードで特定の遺伝子を取り除く
hamster_mouse_orthologs <- hamster_mouse_orthologs %>% filter(!row_number() %in% c(match('GeneX',hamster_mouse_orthologs$Mus_musculus_gene_symbol)))


total_marker_GOlist_down <- total_marker[total_marker$p_val_adj < 0.05 & total_marker$avg_log2FC< -1.0,]

DEG.gene_down <- rownames(total_marker_GOlist_down)
ego_down <- enrichGO(gene = DEG.gene_down,
                     ont ="BP",
                     OrgDb='org.Mm.eg.db', 
                     keyType = "SYMBOL",
                     readable      = TRUE)
dotplot(ego_down)+ggtitle('Total_summer')
ego_down.simple <- simplify(ego_down)
dotplot(ego_down.simple, showCategory=10)


# inhibitory neurons
sumwinHam_merged.neuron_vGat <- subset(sumwinHam_merged.neuron, idents = c("In"))
DimPlot(sumwinHam_merged.neuron_vGat)
sumwinHam_merged.neuron_vGat <- SCTransform(sumwinHam_merged.neuron_vGat, vars.to.regress = c("percent.mt","percent.rp"), return.only.var.genes = TRUE)
sumwinHam_merged.neuron_vGat <- RunPCA(sumwinHam_merged.neuron_vGat,npcs = 150)
ElbowPlot(sumwinHam_merged.neuron_vGat, ndims = 150)
sumwinHam_merged.neuron_vGat <- IntegrateLayers(object = sumwinHam_merged.neuron_vGat, method = CCAIntegration, normalization.method = "SCT",
                                                verbose = FALSE)

sumwinHam_merged.neuron_vGat <- FindNeighbors(sumwinHam_merged.neuron_vGat, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)
sumwinHam_merged.neuron_vGat <- FindClusters(sumwinHam_merged.neuron_vGat, resolution = 1.6)
sumwinHam_merged.neuron_vGat <- RunUMAP(sumwinHam_merged.neuron_vGat, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)

VlnPlot(sumwinHam_merged.neuron_vGat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2, group.by = "SCT_snn_res.1.6")
plots <- VlnPlot(sumwinHam_merged.neuron_vGat, features = c("Camk2a", "Rbfox3","Mog", "Bcas1","Pdgfra","Aqp4","Csf1r","Cldn5",'Col3a1'),slot="data",pt.size = 0.3, split.by = "Mode",stack = T,flip = T)
wrap_plots(plots = plots, ncol = 1)


sumwinHam_merged.neuron_vGat_afterQC <- subset(sumwinHam_merged.neuron_vGat, idents = c(2,1,8), invert = T)
DimPlot(sumwinHam_merged.neuron_vGat_afterQC)
sumwinHam_merged.neuron_vGat_afterQC <- SCTransform(sumwinHam_merged.neuron_vGat_afterQC, vars.to.regress = c("percent.mt","percent.rp"), return.only.var.genes = TRUE)
sumwinHam_merged.neuron_vGat_afterQC <- RunPCA(sumwinHam_merged.neuron_vGat_afterQC,npcs = 150)
ElbowPlot(sumwinHam_merged.neuron_vGat_afterQC, ndims = 150)
sumwinHam_merged.neuron_vGat_afterQC <- IntegrateLayers(object = sumwinHam_merged.neuron_vGat_afterQC, method = CCAIntegration, normalization.method = "SCT",
                                                        verbose = FALSE)
sumwinHam_merged.neuron_vGat_afterQC <- FindNeighbors(sumwinHam_merged.neuron_vGat_afterQC, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)
sumwinHam_merged.neuron_vGat_afterQC <- FindClusters(sumwinHam_merged.neuron_vGat_afterQC, resolution = 1.2)
sumwinHam_merged.neuron_vGat_afterQC <- RunUMAP(sumwinHam_merged.neuron_vGat_afterQC, dims = 1:150, reduction = "integrated.dr", verbose = FALSE)


# Fig6a
DimPlot(sumwinHam_merged.neuron_vGat_afterQC, reduction = "umap", split.by = "Mode",label = T)
DimPlot(sumwinHam_merged.neuron_vGat_afterQC, reduction = "umap", group.by = "Mode",label.size=3,
        cols = c(
          "winter" = "#2166AC",
          "summer" = "#B2182B"
        ))

# Fig6b
df<-table(sumwinHam_merged.neuron_vGat_afterQC@meta.data$seurat_clusters,sumwinHam_merged.neuron_vGat_afterQC@meta.data$Mode)
df <- as.data.frame(df) 
df_long <- df %>%
  rownames_to_column("seurat_clusters") 
df_long <- df_long %>%
  group_by(Var2) %>%
  mutate(Fraction = Freq / sum(Freq))
df_long

ggplot(df_long, aes(x = Var2, y = Fraction, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "fraction of each subtype in each condition",
    x = "Condition",
    y = "fraction"
  ) +
  theme_minimal()

# Fig6c
sumwinHam_merged.neuron_vGat_afterQC<-PrepSCTFindMarkers(sumwinHam_merged.neuron_vGat_afterQC)
total_marker <- FindMarkers(sumwinHam_merged.neuron_vGat_afterQC, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "winter" ,ident.2 = "summer", group.by = "Mode", min.pct = 0.25, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)
# ポイントサイズを調整
# point_sizes <- rep(1, nrow(total_marker))
# names(point_sizes) <- rownames(total_marker)
# point_sizes[c("Atg2a","Atg5","Atg7","Atg10","Atg13","Atg14")] <- 6

# upとdownの数
up_genes <- total_marker %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > log2(2)
  )

down_genes <- total_marker %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC < log2(0.5)
  )

n_up <- nrow(up_genes)
n_down <- nrow(down_genes)

n_up
n_down


keyvals <- ifelse(
  total_marker$avg_log2FC > 1 & total_marker$p_val_adj < 0.05,
  "#2166AC",      # Up（赤）
  ifelse(
    total_marker$avg_log2FC < -1 & total_marker$p_val_adj < 0.05,
    "#B2182B",    # Down（青）
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



# Fig6e
rownames(total_marker) <- ifelse(row.names(total_marker) %in% hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol, 
                                 hamster_mouse_orthologs$Mus_musculus_gene_symbol[match(row.names(total_marker), hamster_mouse_orthologs$Mesocricetus_auratus_gene_symbol)], 
                                 row.names(total_marker))
hamster_mouse_orthologs <- hamster_mouse_orthologs %>% filter(!row_number() %in% c(match('GeneX',hamster_mouse_orthologs$Mus_musculus_gene_symbol)))
total_marker_GOlist_down <- total_marker[total_marker$p_val_adj < 0.05 & total_marker$avg_log2FC< -1.0,]
DEG.gene_down <- rownames(total_marker_GOlist_down)
ego_down <- enrichGO(gene = DEG.gene_down,
                     ont ="BP",
                     OrgDb='org.Mm.eg.db', 
                     keyType = "SYMBOL",
                     readable      = TRUE)
dotplot(ego_down)+ggtitle('Total_summer')
ego_down.simple <- simplify(ego_down)
dotplot(ego_down.simple, showCategory=10)


#Fig6d,e
deg_list <- list()
clusters <- levels(Idents(sumwinHam_merged.neuron_vGat_afterQC))

for(cl in clusters){
  
  deg <- FindMarkers(
    sumwinHam_merged.neuron_vGat_afterQC,
    subset.ident = cl,
    group.by = "Mode",
    ident.1 = "winter",
    ident.2 = "summer",
    recorrect_umi = FALSE,
    logfc.threshold = 0,
    min.pct = 0
  )
  
  deg <- deg[
    abs(deg$avg_log2FC) > 1 &
      deg$p_val_adj < 0.05,
  ]
  
  deg$gene <- rownames(deg)
  
  deg_list[[cl]] <- deg
}
deg_direction <- do.call(
  rbind,
  lapply(names(deg_list), function(cl){
    
    tmp <- deg_list[[cl]]
    
    if(nrow(tmp) == 0){
      return(NULL)
    }
    
    tmp$Cluster <- cl
    
    tmp$Direction <- ifelse(
      tmp$avg_log2FC > 0,
      "Winter",
      "Summer"
    )
    
    tmp[, c(
      "gene",
      "Cluster",
      "Direction"
    )]
  })
)

gene_freq <- deg_direction %>%
  distinct(
    gene,
    Cluster,
    Direction
  ) %>%
  count(
    gene,
    Direction,
    name = "n_cluster"
  )

n_cluster <- length(clusters)

shared_cutoff <- ceiling(
  0.50 * n_cluster
)

specific_cutoff <- floor(
  0.10 * n_cluster
)

gene_freq$Category <- case_when(
  gene_freq$n_cluster >= shared_cutoff ~
    "Shared",

  gene_freq$n_cluster <= specific_cutoff ~
    "Cell type-specific",

  TRUE ~
    "Intermediate"
)
deg_direction <- merge(
  deg_direction,
  gene_freq,
  by = c(
    "gene",
    "Direction"
  )
)

plot_df2 <- deg_direction %>%
  count(
    Cluster,
    Direction,
    Category
  )

plot_df2$N <- plot_df2$n

plot_df2$N[
  plot_df2$Direction == "Summer"
] <-
  -plot_df2$N[
    plot_df2$Direction == "Summer"
  ]

plot_df2$Cluster <- factor(
  plot_df2$Cluster,
  levels = as.character(
    sort(
      as.numeric(
        unique(plot_df2$Cluster)
      )
    )
  )
)

plot_df2$FillGroup <- paste(
  plot_df2$Direction,
  plot_df2$Category,
  sep = "_"
)

ggplot(
  plot_df2,
  aes(
    x = Cluster,
    y = N,
    fill = FillGroup
  )
) +
  geom_col(width = 0.8) +
  
  geom_hline(yintercept = 0) +
  scale_fill_manual(
    values = c(
      
      # Winter系（赤）
      "Winter_Shared" = "#2166AC",
      "Winter_Intermediate" = "#67A9CF",
      "Winter_Cell type-specific" = "#D1E5F0",
      
      # Summer系（青）
      "Summer_Shared" = "#B2182B",
      "Summer_Intermediate" = "#EF8A62",
      "Summer_Cell type-specific" = "#FDDBC7"
    )
  ) +
  theme_classic()
