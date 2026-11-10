# Ex clusterをmerged.objに入れる
ex_clusters <- as.character(Idents(merged.obj_EX))
names(ex_clusters) <- colnames(merged.obj_EX)

merged.obj$Ex_cluster <- NA
merged.obj$Ex_cluster[colnames(merged.obj_EX)] <- ex_clusters

# ① 座標取得
coords <- merged.obj@images[["fov.Ctrl2"]]@boundaries[["centroids"]]@coords

# ② FOVに対応するcell名取得
fov_cells <- Cells(merged.obj[["fov.Ctrl2"]])

# ③ data.frame作成 & cell名付与（順番対応）
df <- as.data.frame(coords)
df$cell <- fov_cells

# ④ cluster情報付与
df$Ex_cluster <- merged.obj$Ex_cluster[df$cell]
df$Ex_cluster <- as.character(df$Ex_cluster)

# ⑤ 指定クラスター
clusters_use <- c("25","27","29")

# ⑥ ハイライトラベル作成
df$highlight <- ifelse(df$Ex_cluster %in% clusters_use,
                       df$Ex_cluster,
                       "other")

# ⑦ プロット

# Ctrl2,0um 1750^1500
ggplot(df, aes(x = y, y = x)) +
  geom_point(color = "grey85", size = 0.3) +
  geom_point(data = subset(df, highlight != "other"),
             aes(color = highlight),
             size = 0.5) +
  theme_classic() +
  coord_fixed(
    xlim = c(12250, 14000),
    ylim = c(9300, 10800),
    ratio = 1
  ) 

# Ctrl2,150um
ggplot(df, aes(x = y, y = x)) +
  geom_point(color = "grey85", size = 0.3) +
  geom_point(data = subset(df, highlight != "other"),
             aes(color = highlight),
             size = 0.5) +
  theme_classic() +
  coord_fixed(
    xlim = c(8000, 9750),
    ylim = c(400, 1900),
    ratio = 1
  ) 
# Ctrl2,300um
ggplot(df, aes(x = y, y = x)) +
  geom_point(color = "grey85", size = 0.3) +
  geom_point(data = subset(df, highlight != "other"),
             aes(color = highlight),
             size = 0.5) +
  theme_classic() +
  coord_fixed(
    xlim = c(7875, 9625),
    ylim = c(2600, 4100),
    ratio = 1
  ) 
# Ctrl2,450um 1750^1500
ggplot(df, aes(x = y, y = x)) +
  geom_point(color = "grey85", size = 0.3) +
  geom_point(data = subset(df, highlight != "other"),
             aes(color = highlight),
             size = 0.5) +
  theme_classic() +
  coord_fixed(
    xlim = c(8125, 9875),
    ylim = c(4950, 6450),
    ratio = 1
  ) 
# Ctrl2,600um 1750^1500
ggplot(df, aes(x = y, y = x)) +
  geom_point(color = "grey85", size = 0.3) +
  geom_point(data = subset(df, highlight != "other"),
             aes(color = highlight),
             size = 0.5) +
  theme_classic() +
  coord_fixed(
    xlim = c(8250, 10000),
    ylim = c(7250, 8750),
    ratio = 1
  ) 
# Ctrl2,750um 1750^1500
ggplot(df, aes(x = y, y = x)) +
  geom_point(color = "grey85", size = 0.3) +
  geom_point(data = subset(df, highlight != "other"),
             aes(color = highlight),
             size = 0.5) +
  theme_classic() +
  coord_fixed(
    xlim = c(8525, 10275),
    ylim = c(9300, 10800),
    ratio = 1
  ) 
