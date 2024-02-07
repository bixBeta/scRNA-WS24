# load in the unintegrated + seurat object ----

sobj <- readRDS("RDS/02_UnIntegrated_res0.8.RDS")
sobj.meta <- sobj@meta.data
sobj


all.genes = rownames(sobj)
var.feats = VariableFeatures(sobj)
DefaultAssay(sobj) <- "RNA"


sobj <- IntegrateLayers(object = sobj, method = CCAIntegration, orig.reduction = "pca", 
                        new.reduction = "pca_cca",
                        verbose = T)

sobj <- IntegrateLayers(object = sobj, method = HarmonyIntegration, orig.reduction = "pca", 
                        new.reduction = "harmony_pca",
                        verbose = T)


# re-join layers after integration ----
sobj[["RNA"]] <- JoinLayers(sobj[["RNA"]])

sobj.meta <- sobj@meta.data

# Run CCA ----
sobj <- FindNeighbors(sobj, reduction = "pca_cca", dims = 1:50)
sobj <- FindClusters(sobj, resolution = 0.8, cluster.name = "cca_clusters")
sobj <- RunUMAP(sobj, dims = 1:50, reduction = "pca_cca", reduction.name = "integrated.cca")


# Run HARMONY ----
sobj <- FindNeighbors(sobj, reduction = "harmony_pca", dims = 1:50)
sobj <- FindClusters(sobj, resolution = 0.8, cluster.name = "harmony_clusters")
sobj <- RunUMAP(sobj, dims = 1:50, reduction = "harmony_pca", reduction.name = "harmony.integrated")

sobj.meta <- sobj@meta.data


# UMAP Plots ----
u1 = DimPlot(sobj, reduction = "umap.unintegrated", group.by = c("orig.ident"), pt.size = .1)
u2 = DimPlot(sobj, reduction = "umap.unintegrated", group.by = c("orig.ident", "unintegrated_clusters"), pt.size = .1, split.by = "orig.ident")

i1 = DimPlot(sobj, reduction = "integrated.cca", group.by = c("orig.ident"), pt.size = .1)
i2 = DimPlot(sobj, reduction = "integrated.cca", group.by = c("orig.ident", "cca_clusters"), pt.size = .1, split.by = "orig.ident")

h1 = DimPlot(sobj, reduction = "harmony.integrated", group.by = c("orig.ident"), pt.size = .1)
h2 = DimPlot(sobj, reduction = "harmony.integrated", group.by = c("orig.ident", "harmony_clusters"), pt.size = .1, split.by = "orig.ident")


png(filename = "figures/UnIntegrated+CCA+HARMONY__UMAPS__res0.8.png", width = 2600, height = 3000, res = 150)
(u1 + u2) / (i1 +i2) / (h1 + h2)
dev.off()


saveRDS(sobj, "RDS/05__HARMONY__CCA__UNIntegrated__res0.8.RDS")


# FindALLMarkers ----
table(Idents(sobj) == sobj$harmony_clusters)

all.markers.res = FindAllMarkers(object = sobj, assay = "RNA", min.pct = 0.1, only.pos = T)

saveRDS(all.markers.res, "RDS/all.markers.results.RDS")

