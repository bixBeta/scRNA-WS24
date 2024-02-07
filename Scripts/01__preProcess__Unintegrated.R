# setup ----
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)

dirs = list.files(path = "10X_Data", full.names = T, pattern = ".h5$",recursive = T)

# import cellranger h5's ----
h5.list = list()

for (i in 1:length(dirs)) {
  h5.list[[i]] <- Read10X_h5(filename = dirs[i] )
  names(h5.list)[[i]] <- strsplit(dirname(dirs[i]), split = "/")[[1]][2]
}


# create seurat objects RNA assay ----
sobj.list = lapply(h5.list, FUN = function(x){
  
  x = CreateSeuratObject(counts = x, min.cells = 3, min.features = 200, assay = "RNA")
  
})


for (i in 1:length(sobj.list)) {
  
  sobj.list[[i]]$orig.ident <- names(sobj.list)[i]  
}

# add % Mito ----
sobj.list = lapply(sobj.list, function(x){
  AddMetaData(object = x, metadata = PercentageFeatureSet(x, "^MT-"), col.name = "percent.mt")
})

# add log10GenesPerUMI ----
sobj.list = lapply(sobj.list, function(x){
  AddMetaData(object = x, metadata = log10(x$nFeature_RNA) / log10(x$nCount_RNA), col.name = "log10GenesPerUMI")
})

# filter seurat object nFeature > 500 < 5000 etc. ----
sobj.list.filtered = lapply(sobj.list, function(x){
  x = subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 30  & log10GenesPerUMI > 0.80)
  
})

sobj.merged.filtered <- merge(x = sobj.list.filtered[[1]], y = sobj.list.filtered[2:length(sobj.list.filtered)], merge.data=TRUE)


saveRDS(sobj.merged.filtered, "RDS/01_sobj.merged.filtered.RDS")



# seurat pipeline on merged sobj ----

sobj <- readRDS("RDS/01_sobj.merged.filtered.RDS")
sobj

# run standard anlaysis workflow
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", verbose = T)
sobj <- FindVariableFeatures(sobj, nfeatures = 3000)

all.genes <- rownames(sobj)

sobj <- ScaleData(sobj, features = all.genes)
sobj <- RunPCA(sobj, npcs = 50, verbose = T)

sobj <- FindNeighbors(sobj, dims = 1:50, reduction = "pca")
sobj <- FindClusters(sobj, resolution = 0.8, cluster.name = "unintegrated_clusters", verbose = T)
sobj <- RunUMAP(sobj, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
sobj.meta = sobj@meta.data

png(filename = "figures/UnIntegrated_UMAP_res0.8.png", width = 2920, height = 1080,res = 150)
DimPlot(sobj, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"), pt.size = .1)
dev.off()

saveRDS(sobj, "RDS/02_UnIntegrated_res0.8.RDS")
