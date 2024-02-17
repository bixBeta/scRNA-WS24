sobj <- readRDS("git-repo/RDS/01_sobj.merged.RDS")

VlnPlot(sobj, features = c("nCount_RNA", "nFeature_RNA" , "log10GenesPerUMI", "percent.mt"),
        pt.size = 0, group.by = "orig.ident", ncol = 4)

metadata <- sobj@meta.data

# lets save all of our plots to a variable, so we can plot them together in one plot (the last line in this code chunk) ----

n_UMI = metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("nCount_RNA or nUMI's (logScale)")

n_Feature = metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("nFeature_RNA (logScale)")



log10GenePerUMI = metadata %>% 
  ggplot(aes(color=orig.ident, x=log10GenesPerUMI, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  ylab("Cell density") +
  xlab("Novelty Score (logScale)")


percent.mt =  metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("Mitochodrial Percentage (logScale)")




(n_UMI + n_Feature) / (log10GenePerUMI + percent.mt)

sobj.filtered <- subset(sobj, subset = nFeature_RNA > 1000 & nFeature_RNA < 9000 & percent.mt < 20 & log10GenesPerUMI > 0.80)
# run standard anlaysis workflow ---- 
# Using LogNormalize Method to Normalize our data ----

sobj.filtered <- NormalizeData(sobj.filtered, normalization.method = "LogNormalize", verbose = T)

# Finding Top 3000 High Variable features using the Variance Stabilized Transformation ----

sobj.filtered <- FindVariableFeatures(sobj.filtered, nfeatures = 3000, selection.method = "vst", verbose = T)

# Scaling (Calculating Z-scores), these are helpful to use for plotting e.g. Heatmaps ----
# By default, ScaleData only scales variable features, but we can do this for all genes 
# This helps with visualization of genes that are not part of HVG's 

all.genes <- rownames(sobj.filtered)
sobj.filtered <- ScaleData(sobj.filtered, features = all.genes)

# Running Principal Components Analysis ---- 
sobj.filtered <- RunPCA(sobj.filtered, npcs = 50, verbose = T)

ElbowPlot(sobj.filtered, ndims = 50, reduction = "pca")

# Calculating Shared Nearest Neighbors (SNN's) ----
sobj.filtered <- FindNeighbors(sobj.filtered, dims = 1:50, reduction = "pca")

# Finding Clusters using the Louvain Clustering Algorithm (default) followed by UMAP ----
sobj.filtered <- FindClusters(sobj.filtered, resolution = 0.4, cluster.name = "unintegrated_clusters", verbose = T)
sobj.filtered <- RunUMAP(sobj.filtered, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")

# Plotting UMAP
DimPlot(sobj.filtered, group.by = "seurat_clusters", label = T)


