deg_list <- list()

clusters <- levels(Idents(MPN_merged.neuron_vGluT2_QC))

for (cl in clusters) {
  
  # クラスター内の細胞を取得
  obj_sub <- subset(MPN_merged.neuron_vGluT2_QC, idents = cl)
  
  # Speciesごとの細胞数
  sp_count <- table(obj_sub$Species)
  
  hamster_n <- ifelse("hamster" %in% names(sp_count), sp_count["hamster"], 0)
  mouse_n   <- ifelse("mouse" %in% names(sp_count), sp_count["mouse"], 0)
  
  # どちらかが3細胞未満ならスキップ
  if (hamster_n < 3 | mouse_n < 3) {
    message(
      "Skipping cluster ", cl,
      " (hamster=", hamster_n,
      ", mouse=", mouse_n, ")"
    )
    next
  }
  
  deg <- FindMarkers(
    MPN_merged.neuron_vGluT2_QC,
    subset.ident = cl,
    group.by = "Species",
    ident.1 = "hamster",
    ident.2 = "mouse",
    logfc.threshold = 0,
    min.pct = 0
  )
  
  deg <- subset(
    deg,
    abs(avg_log2FC) > 1 &
      p_val_adj < 0.05
  )
  
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
      "hamster",
      "mouse"
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
  plot_df2$Direction == "mouse"
] <-
  -plot_df2$N[
    plot_df2$Direction == "mouse"
  ]

library(dplyr)
library(stringr)

plot_df2 <- plot_df2 %>%
  mutate(
    Cluster = factor(
      Cluster,
      levels = paste0(
        "e",
        sort(as.numeric(str_extract(unique(Cluster), "\\d+")))
      )
    )
  ) %>%
  arrange(Cluster)


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
      "hamster_Shared" = "#F8766D",
      "hamster_Intermediate" = "#F9A39D",
      "hamster_Cell type-specific" = "#FCD4D1",
      
      "mouse_Shared" = "#00BFC4",
      "mouse_Intermediate" = "#66D6D9",
      "mouse_Cell type-specific" = "#BFEDEF"
    ) )+
  theme_classic()




clusters <- levels(Idents(MPN_merged.neuron_vGat_QC))
deg_list <- list()

clusters <- levels(Idents(MPN_merged.neuron_vGat_QC))

for (cl in clusters) {
  
  # クラスター内の細胞を取得
  obj_sub <- subset(MPN_merged.neuron_vGat_QC, idents = cl)
  
  # Speciesごとの細胞数
  sp_count <- table(obj_sub$Species)
  
  hamster_n <- ifelse("hamster" %in% names(sp_count), sp_count["hamster"], 0)
  mouse_n   <- ifelse("mouse" %in% names(sp_count), sp_count["mouse"], 0)
  
  # どちらかが3細胞未満ならスキップ
  if (hamster_n < 3 | mouse_n < 3) {
    message(
      "Skipping cluster ", cl,
      " (hamster=", hamster_n,
      ", mouse=", mouse_n, ")"
    )
    next
  }
  
  deg <- FindMarkers(
    MPN_merged.neuron_vGat_QC,
    subset.ident = cl,
    group.by = "Species",
    ident.1 = "hamster",
    ident.2 = "mouse",
    logfc.threshold = 0,
    min.pct = 0
  )
  
  deg <- subset(
    deg,
    abs(avg_log2FC) > 1 &
      p_val_adj < 0.05
  )
  
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
      "hamster",
      "mouse"
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
  plot_df2$Direction == "mouse"
] <-
  -plot_df2$N[
    plot_df2$Direction == "mouse"
  ]

library(dplyr)
library(stringr)

plot_df2 <- plot_df2 %>%
  mutate(
    Cluster = factor(
      Cluster,
      levels = paste0(
        "i",
        sort(as.numeric(str_extract(unique(Cluster), "\\d+")))
      )
    )
  ) %>%
  arrange(Cluster)

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
      "hamster_Shared" = "#F8766D",
      "hamster_Intermediate" = "#F9A39D",
      "hamster_Cell type-specific" = "#FCD4D1",
      
      "mouse_Shared" = "#00BFC4",
      "mouse_Intermediate" = "#66D6D9",
      "mouse_Cell type-specific" = "#BFEDEF"
    ) )+
  theme_classic()