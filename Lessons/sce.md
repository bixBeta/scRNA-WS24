Another most common data class that is widely used to store single cell assays is known as the `SingleCellExperiment` Class. <br>
This usually comes in handy when dealing with scRNA-seq data outside of the seurat framework. <br> For Example: when running doublet finder packages or integration using methods that are not currently supported within the seurat framework.


To convert a seurat object that is created using seurat v5 or above, please use the following guidleine to properly create the sce data class object within R. 


First we create a downgraded assay class within our seurat object. <br>

> assay3 structure was used in older versions of seurat and does not support layers.

To achieve this, let us first create a new assay named `RNA3` in our seurat object. 

```
library(Seurat)
library(SeuratWrappers)

# create a v3 assay and change the default assay to RNA3
sobj[["RNA3"]] <- as(sobj[["RNA"]], Class = "Assay")
DefaultAssay(sobj.sample) <- "RNA3"
```

In the above chunk when we specify Class = "Assay", the conversion of our assay object happens in a manner that is backwards compatible with older versions of seurat. <br>

Now we use the `as.SingleCellExperiment` function to convert our seurat object into a single cell experiment class object. 

```
# create sce from sobj.sample
sce = as.SingleCellExperiment(sobj, assay = "RNA3")
saveRDS(sce, "SingleCellExperiment_Format.RDS")
```
