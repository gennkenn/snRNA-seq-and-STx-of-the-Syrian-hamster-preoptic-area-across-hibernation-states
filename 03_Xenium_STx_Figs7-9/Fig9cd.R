library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

##############################
# 1. Fos発現データの取得
##############################

fos_df <- FetchData(
  merged.obj_EX,
  vars = c(
    "Fos",
    "Condition",
    "animal",
    "seurat_clusters"
  )
)

##############################
# 2. Fos-positive閾値の計算
##############################

fos_threshold <- mean(
  fos_df$Fos,
  na.rm = TRUE
) +
  2 * sd(
    fos_df$Fos,
    na.rm = TRUE
  )

fos_threshold

fos_df <- fos_df %>%
  mutate(
    Fos_positive = Fos >= fos_threshold,
    
    Condition = factor(
      Condition,
      levels = c("Ctrl", "Winter")
    )
  )

##############################
# 3. cluster × animalごとの
#    Fos-positive fraction
##############################

animal_fos_df <- fos_df %>%
  filter(
    Condition %in% c("Ctrl", "Winter"),
    animal %in% c(
      "Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4",
      "Winter1", "Winter2", "Winter3", "Winter4"
    )
  ) %>%
  group_by(
    seurat_clusters,
    Condition,
    animal
  ) %>%
  summarise(
    n_cells = n(),
    n_Fos_positive = sum(
      Fos_positive,
      na.rm = TRUE
    ),
    Fos_positive_fraction =
      n_Fos_positive / n_cells,
    Fos_positive_percent =
      100 * Fos_positive_fraction,
    .groups = "drop"
  ) %>%
  mutate(
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num,
    Condition,
    animal
  )

animal_fos_df

##############################
# 4. clusterごとの
#    two-sided Welch's t-test
##############################

welch_res <- animal_fos_df %>%
  group_by(
    seurat_clusters
  ) %>%
  nest() %>%
  mutate(
    n_Ctrl = map_int(
      data,
      ~ sum(.x$Condition == "Ctrl")
    ),
    
    n_Winter = map_int(
      data,
      ~ sum(.x$Condition == "Winter")
    ),
    
    test = map(
      data,
      ~ tryCatch(
        t.test(
          Fos_positive_fraction ~ Condition,
          data = .x,
          alternative = "two.sided",
          var.equal = FALSE
        ),
        error = function(e) NULL
      )
    ),
    
    mean_Ctrl = map_dbl(
      data,
      ~ mean(
        .x$Fos_positive_fraction[
          .x$Condition == "Ctrl"
        ],
        na.rm = TRUE
      )
    ),
    
    mean_Winter = map_dbl(
      data,
      ~ mean(
        .x$Fos_positive_fraction[
          .x$Condition == "Winter"
        ],
        na.rm = TRUE
      )
    ),
    
    mean_difference = mean_Winter - mean_Ctrl,
    
    t_statistic = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        unname(.x$statistic)
      }
    ),
    
    df = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        unname(.x$parameter)
      }
    ),
    
    p_value = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        .x$p.value
      }
    )
  ) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    
    significance = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE             ~ "ns"
    ),
    
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num
  ) %>%
  select(
    seurat_clusters,
    n_Ctrl,
    n_Winter,
    mean_Ctrl,
    mean_Winter,
    mean_difference,
    t_statistic,
    df,
    p_value,
    p_adj_BH,
    significance
  )

welch_res

library(dplyr)
library(ggplot2)

##############################
# 平均値とSEMを計算
##############################

bar_df <- animal_fos_df %>%
  group_by(
    seurat_clusters,
    Condition
  ) %>%
  summarise(
    n_animals = n(),
    
    mean_Fos_positive_percent =
      mean(
        Fos_positive_percent,
        na.rm = TRUE
      ),
    
    sd_Fos_positive_percent =
      sd(
        Fos_positive_percent,
        na.rm = TRUE
      ),
    
    sem_Fos_positive_percent =
      sd_Fos_positive_percent / sqrt(n_animals),
    
    .groups = "drop"
  ) %>%
  mutate(
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num,
    Condition
  )

##############################
# clusterの順番を設定
##############################

cluster_order <- bar_df %>%
  distinct(
    seurat_clusters,
    cluster_num
  ) %>%
  arrange(
    cluster_num
  ) %>%
  pull(
    seurat_clusters
  ) %>%
  as.character()

bar_df <- bar_df %>%
  mutate(
    seurat_clusters = factor(
      seurat_clusters,
      levels = cluster_order
    ),
    
    Condition = factor(
      Condition,
      levels = c(
        "Ctrl",
        "Winter"
      )
    )
  )

##############################
# 棒グラフ＋SEM
##############################

p_fos_bar <- ggplot(
  bar_df,
  aes(
    x = seurat_clusters,
    y = mean_Fos_positive_percent,
    fill = Condition
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(
      ymin =
        mean_Fos_positive_percent -
        sem_Fos_positive_percent,
      
      ymax =
        mean_Fos_positive_percent +
        sem_Fos_positive_percent
    ),
    position = position_dodge(
      width = 0.8
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      Ctrl = "#EF8A62",
      Winter = "#67A9CF"
    )
  ) +
  
  labs(
    x = "Seurat cluster",
    y = "Fos-positive cells (%)",
    fill = "Condition"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

p_fos_bar


# IEG score

library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

##############################
# 1. IEG_score1取得
##############################

ieg_df <- FetchData(
  merged.obj_EX,
  vars = c(
    "IEG_score1",
    "Condition",
    "animal",
    "seurat_clusters"
  )
)

##############################
# 2. IEG-positive閾値
##############################

ieg_threshold <- mean(
  ieg_df$IEG_score1,
  na.rm = TRUE
) +
  2 * sd(
    ieg_df$IEG_score1,
    na.rm = TRUE
  )

ieg_threshold

ieg_df <- ieg_df %>%
  mutate(
    IEG_positive = IEG_score1 >= ieg_threshold,
    
    Condition = factor(
      Condition,
      levels = c("Ctrl", "Winter")
    )
  )

##############################
# 3. cluster × animalごとの
#    IEG-positive fraction
##############################

animal_ieg_df <- ieg_df %>%
  filter(
    Condition %in% c("Ctrl", "Winter"),
    animal %in% c(
      "Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4",
      "Winter1", "Winter2", "Winter3", "Winter4"
    )
  ) %>%
  group_by(
    seurat_clusters,
    Condition,
    animal
  ) %>%
  summarise(
    n_cells = n(),
    
    n_IEG_positive = sum(
      IEG_positive,
      na.rm = TRUE
    ),
    
    IEG_positive_fraction =
      n_IEG_positive / n_cells,
    
    IEG_positive_percent =
      100 * IEG_positive_fraction,
    
    .groups = "drop"
  ) %>%
  mutate(
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num,
    Condition,
    animal
  )

animal_ieg_df

##############################
# 4. clusterごとの
#    Welch's t-test
##############################

welch_res <- animal_ieg_df %>%
  group_by(
    seurat_clusters
  ) %>%
  nest() %>%
  mutate(
    
    n_Ctrl = map_int(
      data,
      ~ sum(.x$Condition == "Ctrl")
    ),
    
    n_Winter = map_int(
      data,
      ~ sum(.x$Condition == "Winter")
    ),
    
    test = map(
      data,
      ~ tryCatch(
        t.test(
          IEG_positive_fraction ~ Condition,
          data = .x,
          alternative = "two.sided",
          var.equal = FALSE
        ),
        error = function(e) NULL
      )
    ),
    
    mean_Ctrl = map_dbl(
      data,
      ~ mean(
        .x$IEG_positive_fraction[
          .x$Condition == "Ctrl"
        ],
        na.rm = TRUE
      )
    ),
    
    mean_Winter = map_dbl(
      data,
      ~ mean(
        .x$IEG_positive_fraction[
          .x$Condition == "Winter"
        ],
        na.rm = TRUE
      )
    ),
    
    mean_difference =
      mean_Winter - mean_Ctrl,
    
    t_statistic = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        unname(.x$statistic)
      }
    ),
    
    df = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        unname(.x$parameter)
      }
    ),
    
    p_value = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        .x$p.value
      }
    )
  ) %>%
  ungroup() %>%
  mutate(
    
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    
    significance = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(cluster_num) %>%
  select(
    seurat_clusters,
    n_Ctrl,
    n_Winter,
    mean_Ctrl,
    mean_Winter,
    mean_difference,
    t_statistic,
    df,
    p_value,
    p_adj_BH,
    significance
  )

welch_res

##############################
# 5. 棒グラフ用データ
##############################

bar_df <- animal_ieg_df %>%
  group_by(
    seurat_clusters,
    Condition
  ) %>%
  summarise(
    
    n_animals = n(),
    
    mean_IEG_positive_percent =
      mean(
        IEG_positive_percent,
        na.rm = TRUE
      ),
    
    sd_IEG_positive_percent =
      sd(
        IEG_positive_percent,
        na.rm = TRUE
      ),
    
    sem_IEG_positive_percent =
      sd_IEG_positive_percent /
      sqrt(n_animals),
    
    .groups = "drop"
  ) %>%
  mutate(
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num,
    Condition
  )

##############################
# 6. cluster順設定
##############################

cluster_order <- bar_df %>%
  distinct(
    seurat_clusters,
    cluster_num
  ) %>%
  arrange(cluster_num) %>%
  pull(seurat_clusters) %>%
  as.character()

bar_df <- bar_df %>%
  mutate(
    seurat_clusters = factor(
      seurat_clusters,
      levels = cluster_order
    ),
    
    Condition = factor(
      Condition,
      levels = c(
        "Ctrl",
        "Winter"
      )
    )
  )

##############################
# 7. 棒グラフ
##############################

p_IEG_bar <- ggplot(
  bar_df,
  aes(
    x = seurat_clusters,
    y = mean_IEG_positive_percent,
    fill = Condition
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(
      ymin =
        mean_IEG_positive_percent -
        sem_IEG_positive_percent,
      
      ymax =
        mean_IEG_positive_percent +
        sem_IEG_positive_percent
    ),
    position = position_dodge(
      width = 0.8
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      Ctrl = "#EF8A62",
      Winter = "#67A9CF"
    )
  ) +
  
  labs(
    x = "Seurat cluster",
    y = "IEG-positive cells (%)",
    fill = "Condition"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

p_fos_bar <- ggplot(
  bar_df,
  aes(
    x = seurat_clusters,
    y = mean_IEG_positive_percent,
    fill = Condition
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(
      ymin =
        mean_IEG_positive_percent -
        sem_IEG_positive_percent,
      
      ymax =
        mean_IEG_positive_percent +
        sem_IEG_positive_percent
    ),
    position = position_dodge(
      width = 0.8
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      Ctrl = "#EF8A62",
      Winter = "#67A9CF"
    )
  ) +
  
  labs(
    x = "Seurat cluster",
    y = "IEG-positive cells (%)",
    fill = "Condition"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

p_fos_bar


# inhibitory neuron
library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

##############################
# 1. Fos発現データの取得
##############################

fos_df <- FetchData(
  merged.obj_IN,
  vars = c(
    "Fos",
    "Condition",
    "animal",
    "seurat_clusters"
  )
)

##############################
# 2. Fos-positive閾値の計算
##############################

fos_threshold <- mean(
  fos_df$Fos,
  na.rm = TRUE
) +
  2 * sd(
    fos_df$Fos,
    na.rm = TRUE
  )

fos_threshold

fos_df <- fos_df %>%
  mutate(
    Fos_positive = Fos >= fos_threshold,
    
    Condition = factor(
      Condition,
      levels = c("Ctrl", "Winter")
    )
  )

##############################
# 3. cluster × animalごとの
#    Fos-positive fraction
##############################

animal_fos_df <- fos_df %>%
  filter(
    Condition %in% c("Ctrl", "Winter"),
    animal %in% c(
      "Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4",
      "Winter1", "Winter2", "Winter3", "Winter4"
    )
  ) %>%
  group_by(
    seurat_clusters,
    Condition,
    animal
  ) %>%
  summarise(
    n_cells = n(),
    n_Fos_positive = sum(
      Fos_positive,
      na.rm = TRUE
    ),
    Fos_positive_fraction =
      n_Fos_positive / n_cells,
    Fos_positive_percent =
      100 * Fos_positive_fraction,
    .groups = "drop"
  ) %>%
  mutate(
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num,
    Condition,
    animal
  )

animal_fos_df

##############################
# 4. clusterごとの
#    two-sided Welch's t-test
##############################

welch_res <- animal_fos_df %>%
  group_by(
    seurat_clusters
  ) %>%
  nest() %>%
  mutate(
    n_Ctrl = map_int(
      data,
      ~ sum(.x$Condition == "Ctrl")
    ),
    
    n_Winter = map_int(
      data,
      ~ sum(.x$Condition == "Winter")
    ),
    
    test = map(
      data,
      ~ tryCatch(
        t.test(
          Fos_positive_fraction ~ Condition,
          data = .x,
          alternative = "two.sided",
          var.equal = FALSE
        ),
        error = function(e) NULL
      )
    ),
    
    mean_Ctrl = map_dbl(
      data,
      ~ mean(
        .x$Fos_positive_fraction[
          .x$Condition == "Ctrl"
        ],
        na.rm = TRUE
      )
    ),
    
    mean_Winter = map_dbl(
      data,
      ~ mean(
        .x$Fos_positive_fraction[
          .x$Condition == "Winter"
        ],
        na.rm = TRUE
      )
    ),
    
    mean_difference = mean_Winter - mean_Ctrl,
    
    t_statistic = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        unname(.x$statistic)
      }
    ),
    
    df = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        unname(.x$parameter)
      }
    ),
    
    p_value = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        .x$p.value
      }
    )
  ) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    
    significance = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE             ~ "ns"
    ),
    
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num
  ) %>%
  select(
    seurat_clusters,
    n_Ctrl,
    n_Winter,
    mean_Ctrl,
    mean_Winter,
    mean_difference,
    t_statistic,
    df,
    p_value,
    p_adj_BH,
    significance
  )

welch_res

##############################
# 5. cluster順の設定
##############################

cluster_order <- animal_fos_df %>%
  distinct(
    seurat_clusters,
    cluster_num
  ) %>%
  arrange(
    cluster_num
  ) %>%
  pull(
    seurat_clusters
  ) %>%
  as.character()

animal_fos_df <- animal_fos_df %>%
  mutate(
    seurat_clusters = factor(
      seurat_clusters,
      levels = cluster_order
    )
  )


library(dplyr)
library(ggplot2)

##############################
# 平均値とSEMを計算
##############################

bar_df <- animal_fos_df %>%
  group_by(
    seurat_clusters,
    Condition
  ) %>%
  summarise(
    n_animals = n(),
    
    mean_Fos_positive_percent =
      mean(
        Fos_positive_percent,
        na.rm = TRUE
      ),
    
    sd_Fos_positive_percent =
      sd(
        Fos_positive_percent,
        na.rm = TRUE
      ),
    
    sem_Fos_positive_percent =
      sd_Fos_positive_percent / sqrt(n_animals),
    
    .groups = "drop"
  ) %>%
  mutate(
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num,
    Condition
  )

##############################
# clusterの順番を設定
##############################

cluster_order <- bar_df %>%
  distinct(
    seurat_clusters,
    cluster_num
  ) %>%
  arrange(
    cluster_num
  ) %>%
  pull(
    seurat_clusters
  ) %>%
  as.character()

bar_df <- bar_df %>%
  mutate(
    seurat_clusters = factor(
      seurat_clusters,
      levels = cluster_order
    ),
    
    Condition = factor(
      Condition,
      levels = c(
        "Ctrl",
        "Winter"
      )
    )
  )

##############################
# 棒グラフ＋SEM
##############################

p_fos_bar <- ggplot(
  bar_df,
  aes(
    x = seurat_clusters,
    y = mean_Fos_positive_percent,
    fill = Condition
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(
      ymin =
        mean_Fos_positive_percent -
        sem_Fos_positive_percent,
      
      ymax =
        mean_Fos_positive_percent +
        sem_Fos_positive_percent
    ),
    position = position_dodge(
      width = 0.8
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      Ctrl = "#EF8A62",
      Winter = "#67A9CF"
    )
  ) +
  
  labs(
    x = "Seurat cluster",
    y = "Fos-positive cells (%)",
    fill = "Condition"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

p_fos_bar




##############################
# 1. IEG_score1取得
##############################

ieg_df <- FetchData(
  merged.obj_IN,
  vars = c(
    "IEG_score1",
    "Condition",
    "animal",
    "seurat_clusters"
  )
)

##############################
# 2. IEG-positive閾値
##############################

ieg_threshold <- mean(
  ieg_df$IEG_score1,
  na.rm = TRUE
) +
  2 * sd(
    ieg_df$IEG_score1,
    na.rm = TRUE
  )

ieg_threshold

ieg_df <- ieg_df %>%
  mutate(
    IEG_positive = IEG_score1 >= ieg_threshold,
    
    Condition = factor(
      Condition,
      levels = c("Ctrl", "Winter")
    )
  )

##############################
# 3. cluster × animalごとの
#    IEG-positive fraction
##############################

animal_ieg_df <- ieg_df %>%
  filter(
    Condition %in% c("Ctrl", "Winter"),
    animal %in% c(
      "Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4",
      "Winter1", "Winter2", "Winter3", "Winter4"
    )
  ) %>%
  group_by(
    seurat_clusters,
    Condition,
    animal
  ) %>%
  summarise(
    n_cells = n(),
    
    n_IEG_positive = sum(
      IEG_positive,
      na.rm = TRUE
    ),
    
    IEG_positive_fraction =
      n_IEG_positive / n_cells,
    
    IEG_positive_percent =
      100 * IEG_positive_fraction,
    
    .groups = "drop"
  ) %>%
  mutate(
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num,
    Condition,
    animal
  )

animal_ieg_df

##############################
# 4. clusterごとの
#    Welch's t-test
##############################

welch_res <- animal_ieg_df %>%
  group_by(
    seurat_clusters
  ) %>%
  nest() %>%
  mutate(
    
    n_Ctrl = map_int(
      data,
      ~ sum(.x$Condition == "Ctrl")
    ),
    
    n_Winter = map_int(
      data,
      ~ sum(.x$Condition == "Winter")
    ),
    
    test = map(
      data,
      ~ tryCatch(
        t.test(
          IEG_positive_fraction ~ Condition,
          data = .x,
          alternative = "two.sided",
          var.equal = FALSE
        ),
        error = function(e) NULL
      )
    ),
    
    mean_Ctrl = map_dbl(
      data,
      ~ mean(
        .x$IEG_positive_fraction[
          .x$Condition == "Ctrl"
        ],
        na.rm = TRUE
      )
    ),
    
    mean_Winter = map_dbl(
      data,
      ~ mean(
        .x$IEG_positive_fraction[
          .x$Condition == "Winter"
        ],
        na.rm = TRUE
      )
    ),
    
    mean_difference =
      mean_Winter - mean_Ctrl,
    
    t_statistic = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        unname(.x$statistic)
      }
    ),
    
    df = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        unname(.x$parameter)
      }
    ),
    
    p_value = map_dbl(
      test,
      ~ if (is.null(.x)) {
        NA_real_
      } else {
        .x$p.value
      }
    )
  ) %>%
  ungroup() %>%
  mutate(
    
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    
    significance = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(cluster_num) %>%
  select(
    seurat_clusters,
    n_Ctrl,
    n_Winter,
    mean_Ctrl,
    mean_Winter,
    mean_difference,
    t_statistic,
    df,
    p_value,
    p_adj_BH,
    significance
  )

welch_res

##############################
# 5. 棒グラフ用データ
##############################

bar_df <- animal_ieg_df %>%
  group_by(
    seurat_clusters,
    Condition
  ) %>%
  summarise(
    
    n_animals = n(),
    
    mean_IEG_positive_percent =
      mean(
        IEG_positive_percent,
        na.rm = TRUE
      ),
    
    sd_IEG_positive_percent =
      sd(
        IEG_positive_percent,
        na.rm = TRUE
      ),
    
    sem_IEG_positive_percent =
      sd_IEG_positive_percent /
      sqrt(n_animals),
    
    .groups = "drop"
  ) %>%
  mutate(
    cluster_num =
      as.numeric(
        as.character(seurat_clusters)
      )
  ) %>%
  arrange(
    cluster_num,
    Condition
  )

##############################
# 6. cluster順設定
##############################

cluster_order <- bar_df %>%
  distinct(
    seurat_clusters,
    cluster_num
  ) %>%
  arrange(cluster_num) %>%
  pull(seurat_clusters) %>%
  as.character()

bar_df <- bar_df %>%
  mutate(
    seurat_clusters = factor(
      seurat_clusters,
      levels = cluster_order
    ),
    
    Condition = factor(
      Condition,
      levels = c(
        "Ctrl",
        "Winter"
      )
    )
  )

##############################
# 7. 棒グラフ
##############################
p_fos_bar <- ggplot(
  bar_df,
  aes(
    x = seurat_clusters,
    y = mean_IEG_positive_percent,
    fill = Condition
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(
      ymin =
        mean_IEG_positive_percent -
        sem_IEG_positive_percent,
      
      ymax =
        mean_IEG_positive_percent +
        sem_IEG_positive_percent
    ),
    position = position_dodge(
      width = 0.8
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      Ctrl = "#EF8A62",
      Winter = "#67A9CF"
    )
  ) +
  
  labs(
    x = "Seurat cluster",
    y = "IEG-positive cells (%)",
    fill = "Condition"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

p_fos_bar