# sobj <- readRDS("git-repo/RDS/01_sobj.merged.RDS")

dirs = list.files(path = "10X_Data/", full.names = T, pattern = "filtered_feature_bc_matrix.h5$", recursive = T)

# run dirs to confirm the path of h5 filtered matrices
dirs

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(dplyr)

# first we initialize an empty R list object and load in all h5 matrices per sample as a list using a for loop 

# initialize an empty list 
h5.list = list()

# looping over total number of samples (in this case 3) and storing the h5 matrices in our h5.list object
for (i in 1:length(dirs)) {
  h5.list[[i]] <- Read10X_h5(filename = dirs[i] )
  names(h5.list)[[i]] <- strsplit(dirname(dirs[i]), split = "/")[[1]][3] # might need tweeking for when the 10X outs are outside of Rstudio project dir
}


# create seurat objects RNA assay ----
sobj.list = lapply(h5.list, FUN = function(x){
  
  x = CreateSeuratObject(counts = x, min.cells = 3, min.features = 200, assay = "RNA")
  
})

# by default SampleName which is represented as the orig.ident metadata variable in a seurat object will be named to 'SeuratProject', 
# we will use the following loop to overwrite that and rename with the sampleID

for (i in 1:length(sobj.list)) {
  sobj.list[[i]]$orig.ident <- names(sobj.list)[i]  
}

sobj.list[[1]]$orig.ident %>% head()
sobj.list[[2]]$orig.ident %>% head()
sobj.list[[3]]$orig.ident %>% head()

head(sobj.list[[1]]@meta.data)

# add % Mito ----
sobj.list = lapply(sobj.list, function(x){
  AddMetaData(object = x, metadata = PercentageFeatureSet(x, "^MT-"), col.name = "percent.mt")
})

# add log10GenesPerUMI ----
sobj.list = lapply(sobj.list, function(x){
  AddMetaData(object = x, metadata = log10(x$nFeature_RNA) / log10(x$nCount_RNA), col.name = "log10GenesPerUMI")
})

sobj <- merge(x = sobj.list[[1]], y = sobj.list[2:length(sobj.list)], merge.data=TRUE)

saveRDS(sobj, "git-repo/RDS/01_sobj.merged.RDS")


# sobj <- readRDS("01_sobj.merged.RDS")

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

saveRDS(sobj.filtered, "git-repo/RDS/02_sobj_filtered_res0.4.RDS")
