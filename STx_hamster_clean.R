## Load Libraries
library(Seurat)
library(future)
plan("multisession", workers = 8)
library(ggplot2)
library(dplyr)
library(spacexr)
library(patchwork)
library(Seurat)
library(data.table)
library(sf)
library(geojsonsf)
library(RColorBrewer)



prep_FindMarkers <- function(obj, num_slices = 4) {
  # Loop through total number slices
  for (i in 1:num_slices) {
    slot(object = obj@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
  obj <- obj %>% PrepSCTFindMarkers()
  # Return the modified obj
  return(obj)}


ReadXenium <- function (data.dir, outs = c("matrix", "microns"), type = "centroids", 
                        mols.qv.threshold = 20) 
{
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), 
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"), 
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, 
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid, 
                                     y = cell_info$y_centroid, cell = cell_info$cell_id, 
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir, 
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {
      
      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)
      
      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location, 
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}


## Load Xenium Data
data.ctrl1 <- ReadXenium(data.dir = "./Xenium_hamster/output-XETG00312__0043148__XP010_L01_control1__20250509__063229", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.winter1<- ReadXenium(data.dir = "./Xenium_hamster/output-XETG00312__0043151__XP010_R01_hibernator1__20250509__063230", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.ctrl2 <- ReadXenium("./Xenium_hamster/output-XETG00312__0042638__XP011_L01_ctrl2__20250522__073324", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.winter2 <- ReadXenium("./Xenium_hamster/output-XETG00312__0042641__XP011_R01_hibernator2__20250522__073324", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))


gene_lists<-read.csv("./Xenium_hamster_analysis/4XBGWM_mHamster_MPN_100g_gene_list.csv")

# matrixの遺伝子名変更
matched_indices <- match( rownames(data.ctrl1$matrix[["Gene Expression"]]),gene_lists$Gene)
ensembl_ids <- gene_lists$Ensembl.ID[matched_indices]
ensembl_ids
rownames(data.ctrl1$matrix[["Gene Expression"]]) <- ensembl_ids

matched_indices_microns <- match(data.ctrl1[["microns"]][["gene"]],gene_lists$Gene)

ensembl_ids_microns <- gene_lists$Ensembl.ID[matched_indices_microns]
ensembl_or_gene <- ifelse(is.na(ensembl_ids_microns),
                          data.ctrl1[["microns"]][["gene"]],                 
                          ensembl_ids_microns)   
ensembl_or_gene
data.ctrl1[["microns"]][["gene"]] <- ensembl_or_gene



matched_indices <- match(rownames(data.ctrl2$matrix[["Gene Expression"]]),gene_lists$Gene)
ensembl_ids <- gene_lists$Ensembl.ID[matched_indices]
ensembl_ids
rownames(data.ctrl2$matrix[["Gene Expression"]]) <- ensembl_ids
matched_indices_microns <- match(data.ctrl2[["microns"]][["gene"]],gene_lists$Gene)
ensembl_ids_microns <- gene_lists$Ensembl.ID[matched_indices_microns]
ensembl_or_gene <- ifelse(is.na(ensembl_ids_microns),
                          data.ctrl2[["microns"]][["gene"]],                  
                          ensembl_ids_microns)   

ensembl_or_gene
data.ctrl2[["microns"]][["gene"]] <- ensembl_or_gene




segmentations.data.ctrl1 <- list(
  centroids = CreateCentroids(data.ctrl1$centroids),
  segmentation = CreateSegmentation(data.ctrl1$segmentations))
coords.ctrl1 <- CreateFOV(
  coords = segmentations.data.ctrl1, 
  type = c("segmentation", "centroids"), 
  molecules = data.ctrl1$microns, 
  assay = "Xenium")
xenium.obj.ctrl1 <- CreateSeuratObject(
  counts = data.ctrl1$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "Ctrl1")
xenium.obj.ctrl1[["BlankCodeword"]] <- CreateAssayObject(counts = data.ctrl1$matrix[["Unassigned Codeword"]])
xenium.obj.ctrl1[["ControlCodeword"]] <- CreateAssayObject(counts = data.ctrl1$matrix[["Negative Control Codeword"]])
xenium.obj.ctrl1[["ControlProbe"]] <- CreateAssayObject(counts = data.ctrl1$matrix[["Negative Control Probe"]])
xenium.obj.ctrl1[["fov"]] <- coords.ctrl1
xenium.obj.ctrl1$Condition <- 'Ctrl'

xenium.obj.ctrl1 <- subset(xenium.obj.ctrl1, subset = nCount_Xenium > 0)

print(xenium.obj.ctrl1)




segmentations.data.ctrl2 <- list(
  centroids = CreateCentroids(data.ctrl2$centroids),
  segmentation = CreateSegmentation(data.ctrl2$segmentations))
coords.ctrl2 <- CreateFOV(
  coords = segmentations.data.ctrl2,
  type = c("segmentation", "centroids"),
  molecules = data.ctrl2$microns,
  assay = "Xenium")
xenium.obj.ctrl2 <- CreateSeuratObject(
  counts = data.ctrl2$matrix[["Gene Expression"]],
  assay = "Xenium", project = "Ctrl2")
xenium.obj.ctrl2[["BlankCodeword"]] <- CreateAssayObject(counts = data.ctrl2$matrix[["Unassigned Codeword"]])
xenium.obj.ctrl2[["ControlCodeword"]] <- CreateAssayObject(counts = data.ctrl2$matrix[["Negative Control Codeword"]])
xenium.obj.ctrl2[["ControlProbe"]] <- CreateAssayObject(counts = data.ctrl2$matrix[["Negative Control Probe"]])
xenium.obj.ctrl2[["fov"]] <- coords.ctrl2
xenium.obj.ctrl2$Condition <- 'Ctrl'
xenium.obj.ctrl2 <- subset(xenium.obj.ctrl2, subset = nCount_Xenium > 0)
print(xenium.obj.ctrl2)


matched_indices <- match( rownames(data.winter1$matrix[["Gene Expression"]]),gene_lists$Gene)
ensembl_ids <- gene_lists$Ensembl.ID[matched_indices]
ensembl_ids
rownames(data.winter1$matrix[["Gene Expression"]]) <- ensembl_ids
matched_indices_microns <- match(data.winter1[["microns"]][["gene"]],gene_lists$Gene)
ensembl_ids_microns <- gene_lists$Ensembl.ID[matched_indices_microns]
ensembl_or_gene <- ifelse(is.na(ensembl_ids_microns),
                          data.winter1[["microns"]][["gene"]],                  
                          ensembl_ids_microns)   

ensembl_or_gene
data.winter1[["microns"]][["gene"]] <- ensembl_or_gene


matched_indices <- match(rownames(data.winter2$matrix[["Gene Expression"]]),gene_lists$Gene)
ensembl_ids <- gene_lists$Ensembl.ID[matched_indices]
ensembl_ids
rownames(data.winter2$matrix[["Gene Expression"]]) <- ensembl_ids
matched_indices_microns <- match(data.winter2[["microns"]][["gene"]],gene_lists$Gene)
ensembl_ids_microns <- gene_lists$Ensembl.ID[matched_indices_microns]
ensembl_or_gene <- ifelse(is.na(ensembl_ids_microns),
                          data.winter2[["microns"]][["gene"]],                 
                          ensembl_ids_microns)   
ensembl_or_gene
data.winter2[["microns"]][["gene"]] <- ensembl_or_gene



segmentations.data.winter1 <- list(
  centroids = CreateCentroids(data.winter1$centroids),
  segmentation = CreateSegmentation(data.winter1$segmentations))
coords.winter1 <- CreateFOV(
  coords = segmentations.data.winter1,
  type = c("segmentation", "centroids"),
  molecules = data.winter1$microns,
  assay = "Xenium")
xenium.obj.winter1 <- CreateSeuratObject(
  counts = data.winter1$matrix[["Gene Expression"]],
  assay = "Xenium", project = "Winter1")
xenium.obj.winter1[["BlankCodeword"]] <- CreateAssayObject(counts = data.winter1$matrix[["Unassigned Codeword"]])
xenium.obj.winter1[["ControlCodeword"]] <- CreateAssayObject(counts = data.winter1$matrix[["Negative Control Codeword"]])
xenium.obj.winter1[["ControlProbe"]] <- CreateAssayObject(counts = data.winter1$matrix[["Negative Control Probe"]])
xenium.obj.winter1[["fov"]] <- coords.winter1
xenium.obj.winter1$Condition <- 'Winter'
xenium.obj.winter1 <- subset(xenium.obj.winter1, subset = nCount_Xenium > 0)

print(xenium.obj.winter1)

segmentations.data.winter2 <- list(
  centroids = CreateCentroids(data.winter2$centroids),
  segmentation = CreateSegmentation(data.winter2$segmentations))
coords.winter2 <- CreateFOV(
  coords = segmentations.data.winter2,
  type = c("segmentation", "centroids"),
  molecules = data.winter2$microns,
  assay = "Xenium")
xenium.obj.winter2 <- CreateSeuratObject(
  counts = data.winter2$matrix[["Gene Expression"]],
  assay = "Xenium", project = "Winter2")
xenium.obj.winter2[["BlankCodeword"]] <- CreateAssayObject(counts = data.winter2$matrix[["Unassigned Codeword"]])
xenium.obj.winter2[["ControlCodeword"]] <- CreateAssayObject(counts = data.winter2$matrix[["Negative Control Codeword"]])
xenium.obj.winter2[["ControlProbe"]] <- CreateAssayObject(counts = data.winter2$matrix[["Negative Control Probe"]])
xenium.obj.winter2[["fov"]] <- coords.winter2
xenium.obj.winter2$Condition <- 'Winter'
xenium.obj.winter2 <- subset(xenium.obj.winter2, subset = nCount_Xenium > 0)

print(xenium.obj.winter2)


## Integration#
merged.obj <- merge(x = xenium.obj.ctrl1, y = c(xenium.obj.winter1,xenium.obj.ctrl2,xenium.obj.winter2), add.cell.ids = c("Ctrl1","Winter1","Ctrl2","Winter2"))

## SCTransform, PCA, UMAP, Clustering
options(future.globals.maxSize = 64 * 1024 ^ 3) 
merged.obj <- SCTransform(merged.obj, assay = "Xenium")
merged.obj <- RunPCA(merged.obj, npcs = 50, features = rownames(merged.obj))
merged.obj <- RunUMAP(merged.obj, dims = 1:50)
merged.obj <- FindNeighbors(merged.obj, dims = 1:50)
merged.obj <- FindClusters(merged.obj, resolution = 0.3)



# くらすたー情報の保存（Ependymalはまーかー遺伝子ではなく、局在で判断した）
merged.obj <- RenameIdents(merged.obj,'0'='Astrocyte','1'='MO','2'='Ex','3'='Microglia','4'='IO','5'='In','6'='In','7'='In','8'='Fibroblast','9'='Endothelial','10'='In','11'='In','12'='Endothelial','13'='In','14'='Ambiguous','15'='Ambiguous','16'='Ependymal','17'='Ambiguous','18'='Ambiguous','19'='Ex','20'='In','21'='In','22'='Ambiguous','23'='MO','24'='Ex','25'='Endothelial','26'='Ex','27'='MO','28'='In','29'='Ex','30'='Ambiguous','31'='Ambiguous','32'='Ambiguous','33'='Ambiguous','34'='Ambiguous','35'='Ambiguous','36'='Ambiguous')
merged.obj[["Cell_class"]] <- Idents(merged.obj)


merged.obj$Cell_class <- factor(
  merged.obj$Cell_class,
  levels = c("Ex", "In", "MO", "IO", "Microglia", "Astrocyte", "Fibroblast", "Endothelial", "Ependymal", "Ambiguous")
)
merged.obj@active.ident <- merged.obj$Cell_class

# Fig7d
DimPlot(merged.obj, reduction = "umap", label = F,split.by = "Condition")
DimPlot(merged.obj, reduction = "umap", group.by = "Condition",label.size=3,
        cols = c(
          "Winter" = "#67A9CF",
          "Ctrl" = "#EF8A62"
        ))

# Fig7e
VlnPlot(merged.obj, features = c("Camk2a", "Rbfox3", "Slc17a6",  "Gad1","Slc32a1","Plp1","Opalin","Pdgfra", "Bcas1","Csf1r","Aqp4","Slc39a12","Col3a1","Cldn5","Adgrv1"), 
        pt.size = 0, split.by = "Condition",flip = T, stack = T)


# Fig7g
df<-table(merged.obj@meta.data$Cell_class,merged.obj@meta.data$Condition)
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
    x = "Condition",
    y = "fraction"
  ) +
  theme_minimal()









# inhibitory neurons
merged.obj_IN <- subset(merged.obj, idents = "In")
merged.obj_IN <- SCTransform(merged.obj_IN, assay = "Xenium")
merged.obj_IN <- RunPCA(merged.obj_IN, npcs = 50, features = rownames(merged.obj_IN))
merged.obj_IN <- RunUMAP(merged.obj_IN, dims = 1:50)
merged.obj_IN <- FindNeighbors(merged.obj_IN, dims = 1:50)
merged.obj_IN <- FindClusters(merged.obj_IN, resolution = 0.9)

# Fig7f
DimPlot(merged.obj_IN, reduction = "umap", label = F,split.by="Condition")
DimPlot(merged.obj_IN, reduction = "umap", group.by = "Condition",label.size=3,
        cols = c(
          "Winter" = "#67A9CF",
          "Ctrl" = "#EF8A62"
        ))

merged.obj_IN <- AddModuleScore(object = merged.obj_IN,
                                     features = list(c("Arc","Fos","Fosb","Nr4a1","Egr1","Homer1","Junb","Dusp1")),
                                     name = "IEG_score",ctrl = 4)

# Fig7g
df<-table(merged.obj_IN@meta.data$seurat_clusters,merged.obj_IN@meta.data$Condition)
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
    title = "各Conditionにおける細胞種の割合",
    x = "Condition",
    y = "割合"
  ) +
  theme_minimal()



# Fig9a
FeaturePlot(merged.obj_IN,features = c("Fos"), split.by = "Condition")


# Fig9b
merged.obj_IN<-prep_FindMarkers(merged.obj_IN)
vGat_marker <- FindMarkers(merged.obj_IN, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "Winter" ,ident.2 = "Ctrl", group.by = "Condition", min.pct = 0, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

EnhancedVolcano(vGat_marker,
                lab = rownames(vGat_marker),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'Fold change'),
                ylab = bquote(~Log[10]~ 'p val'),
                pCutoff = 0.1*10^(-2),
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 5/5,
                
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')


# excitatoy neuron
merged.obj_EX <- subset(merged.obj, idents = c("Ex"))
merged.obj_EX <- SCTransform(merged.obj_EX, assay = "Xenium")
merged.obj_EX <- RunPCA(merged.obj_EX, npcs = 50, features = rownames(merged.obj_EX))
merged.obj_EX <- RunUMAP(merged.obj_EX, dims = 1:50)
merged.obj_EX <- FindNeighbors(merged.obj_EX, dims = 1:50)
merged.obj_EX <- FindClusters(merged.obj_EX, resolution = 0.6)
DefaultAssay(merged.obj_EX) <- "SCT"


# Fig7f
DimPlot(merged.obj_EX, reduction = "umap", split.by ="Condition",shuffle = TRUE , label=FALSE)
DimPlot(merged.obj_EX, reduction = "umap", group.by = "Condition",label.size=3,
        cols = c(
          "Winter" = "#67A9CF",
          "Ctrl" = "#EF8A62"
        ))

merged.obj_EX <- AddModuleScore(object = merged.obj_EX,
                                features = list(c("Arc","Fos","Fosb","Nr4a1","Egr1","Homer1","Junb","Dusp1")),
                                name = "IEG_score",ctrl = 4)

# Fig9a
FeaturePlot(merged.obj_EX,"Fos",split.by="Condition", col=c("aquamarine", "tomato2"),  pt.size = 0.2, min.cutoff = 0)

# Fig9b
merged.obj_EX <-prep_FindMarkers(merged.obj_EX)
vGluT2_marker <- FindMarkers(merged.obj_EX, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "Winter" ,ident.2 = "Ctrl", group.by = "Condition", min.pct = 0, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

EnhancedVolcano(vGluT2_marker,
                lab = rownames(vGluT2_marker),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'Fold change'),
                ylab = bquote(~Log[10]~ 'p val'),
                pCutoff = 0.1*10^(-2),
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 5/5,
                
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')



# Fig7g
df<-table(merged.obj_EX@meta.data$seurat_clusters,merged.obj_EX@meta.data$Condition)
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
    title = "各Conditionにおける細胞種の割合",
    x = "Condition",
    y = "割合"
  ) +
  theme_minimal()


# Fig8a
#To add the “animal” column to the metadata, see 20260710_ROIextraction.Rmd.
plan('sequential')
# Peudobulkにして発現変動遺伝子の検定
pseudo_merged.obj_EX <- AggregateExpression(merged.obj_EX, assays = "SCT", return.seurat = T, group.by = c("Condition", "animal", "seurat_clusters"))
pseudo_merged.obj_EX$celltype.condition <- paste(pseudo_merged.obj_EX$seurat_clusters, pseudo_merged.obj_EX$Condition, sep = "_")
Idents(pseudo_merged.obj_EX) <- "celltype.condition"

library(Seurat)
library(dplyr)
library(tidyr)

##############################
# 対象クラスター
##############################

target_clusters <- c(
  3,4,5,6,20,21,23,25,
  27,28,29,32,33,34,37,41
)

##############################
# 対象遺伝子
##############################

target_genes <- c(
  "Adcyap1",
  "Bdnf",
  "Lepr",
  "Ptger3",
  "Opn5",
  "Brs3"
)

##############################
# 発現データ取得
##############################

expr <- FetchData(
  pseudo_merged.obj_EX,
  vars = c(
    target_genes,
    "Condition",
    "animal",
    "seurat_clusters"
  )
)

##############################
# 結果保存用
##############################

results <- list()

##############################
# ループ
##############################

for(cl in target_clusters){
  
  dat_cluster <- expr %>%
    filter(seurat_clusters == cl)
  
  for(gene in target_genes){
    
    dat <- dat_cluster %>%
      select(
        Condition,
        animal,
        Expression = all_of(gene)
      )
    
    ## Condition確認
    if(length(unique(dat$Condition)) != 2){
      next
    }
    
    ## Wilcoxon
    wt <- wilcox.test(
      Expression ~ Condition,
      data = dat,
      exact = FALSE
    )
    
    ## 出力
    results[[length(results)+1]] <- data.frame(
      
      Cluster = cl,
      Gene = gene,
      
      Ctrl_mean =
        mean(dat$Expression[
          dat$Condition=="Ctrl"
        ]),
      
      Winter_mean =
        mean(dat$Expression[
          dat$Condition=="Winter"
        ]),
      
      Ctrl_SD =
        sd(dat$Expression[
          dat$Condition=="Ctrl"
        ]),
      
      Winter_SD =
        sd(dat$Expression[
          dat$Condition=="Winter"
        ]),
      
      Ctrl1 =
        dat$Expression[
          dat$animal=="Ctrl1"
        ],
      
      Ctrl2 =
        dat$Expression[
          dat$animal=="Ctrl2"
        ],
      
      Ctrl3 =
        dat$Expression[
          dat$animal=="Ctrl3"
        ],
      
      Ctrl4 =
        dat$Expression[
          dat$animal=="Ctrl4"
        ],
      
      Winter1 =
        dat$Expression[
          dat$animal=="Winter1"
        ],
      
      Winter2 =
        dat$Expression[
          dat$animal=="Winter2"
        ],
      
      Winter3 =
        dat$Expression[
          dat$animal=="Winter3"
        ],
      
      Winter4 =
        dat$Expression[
          dat$animal=="Winter4"
        ],
      
      p.value = wt$p.value
      
    )
    
  }
  
}

##############################
# 結果結合
##############################

result_df <- bind_rows(results)


##############################
# Bonferroni補正
##############################

result_df <- result_df %>%
  mutate(
    p.value.bonferroni = p.adjust(
      p.value,
      method = "bonferroni"
    )
  )

##############################
# 表示
##############################

result_df

##############################
# 保存
##############################

write.csv(
  result_df,
  "PseudoBulk_Wilcoxon_Results.csv",
  row.names = FALSE
)




# 描画用コード
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

##############################
## 設定
##############################

target_clusters <- c(
  3,4,5,6,20,21,23,25,
  27,28,29,32,33,34,37,41
)

target_genes <- c(
  "Adcyap1",
  "Bdnf",
  "Lepr",
  "Ptger3",
  "Opn5",
  "Brs3"
)

##############################
## データ取得
##############################

plot_df <- FetchData(
  pseudo_merged.obj_EX,
  vars = c(
    target_genes,
    "Condition",
    "animal",
    "seurat_clusters"
  )
)

##############################
## long形式へ
##############################

plot_df <- plot_df %>%
  filter(seurat_clusters %in% target_clusters) %>%
  pivot_longer(
    cols = all_of(target_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

##############################
## factor順
##############################

plot_df$seurat_clusters <- factor(
  plot_df$seurat_clusters,
  levels = target_clusters
)

plot_df$Gene <- factor(
  plot_df$Gene,
  levels = target_genes
)

plot_df$Condition <- factor(
  plot_df$Condition,
  levels = c("Ctrl","Winter")
)

##############################
## 平均・SEM
##############################

summary_df <- plot_df %>%
  group_by(
    Gene,
    seurat_clusters,
    Condition
  ) %>%
  summarise(
    Mean = mean(Expression),
    SEM = sd(Expression) / sqrt(n()),
    .groups = "drop"
  )

##############################
## Plot
##############################

p <- ggplot(
  plot_df,
  aes(
    x = Condition,
    y = Expression,
    colour = Condition
  )
) +
  
  ## animalごとの点
  geom_quasirandom(
    width = 0.18,
    size = 2.6,
    alpha = 1
  ) +
  
  ## SEM
  geom_errorbar(
    data = summary_df,
    aes(
      x = Condition,
      ymin = Mean - SEM,
      ymax = Mean + SEM
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.5
  ) +
  
  ## Mean
  geom_crossbar(
    data = summary_df,
    aes(
      x = Condition,
      y = Mean,
      ymin = Mean,
      ymax = Mean
    ),
    inherit.aes = FALSE,
    width = 0.45,
    linewidth = 0.6,
    colour = "black"
  ) +
  
  facet_grid(
    Gene ~ seurat_clusters,
    scales = "free_y"
  ) +
  
  ## 必ず0スタート
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  scale_colour_manual(
    values = c(
      Ctrl = "#EF8A62",
      Winter = "#67A9CF"
    )
  ) +
  
  labs(
    x = NULL,
    y = "SCT expression"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    
    legend.position = "none",
    
    strip.background = element_blank(),
    
    strip.text.x = element_text(
      face = "bold",
      size = 10
    ),
    
    strip.text.y = element_text(
      face = "italic",
      size = 10
    ),
    
    panel.spacing.x = unit(
      0.15,
      "lines"
    ),
    
    panel.spacing.y = unit(
      0.35,
      "lines"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    
    axis.text.y = element_text(
      size = 8
    ),
    
    axis.title.y = element_text(
      size = 11
    ),
    
    axis.title.x = element_blank()
    
  )

p