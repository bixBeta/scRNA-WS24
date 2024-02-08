> [!Tip]
Useful Links:
[Seurat Version 5 Command Cheat Sheet](https://satijalab.org/seurat/articles/essential_commands)


### Cellranger Output

This week we start with the cellranger outputs that we created in the previous exercise. <br>
Conveniently, the output folders from cellranger are also available at the following path: `/path/to/cellranger-outs/`

It is important to note that **cellranger** creates a new directory for each sample and there are many files/subfolders that are created within that directory.
The count matrices produced by cellranger are located in the `./sampleName/outs/` directory. <br> For Seurat analysis we will be using the file named `filtered_feature_bc_matrix.h5` for each sample.

<hr>
<br>


# 1. Create Seurat Object

In an Rstudio session create a new project.

Create a new R script using the File --> New File --> Rscript or Shift/Cmd|Ctrl/N <br> 
To get started, we will need the path where the 10x outputs are present. Copy and paste the following line to specify the directory of cellranger outs:


```
# this assumes your 10x data is in Rstudio_Project/10X_Data/sampleName/outs/ directory 
# edit the path argument as needed 

dirs = list.files(path = "10X_Data", full.names = T, pattern = ".h5$", recursive = T)

# run dirs to confirm the path of h5 filtered matrices
dirs

```

Use the following code chunk to load all the R libraries that are needed to run the analysis.

```
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(dplyr)

```


Once all libraries are loaded, we will create a seurat object from all samples using `Seurat::Read10X_h5()` and `SeuratObject::CreateSeuratObject()` functions. 

```
# first we initialize an empty R list object and load in all h5 matrices per sample as a list using a for loop 

# initialize an empty list 
h5.list = list()

# looping over total number of samples (in this case 3) and storing the h5 matrices in our h5.list object
for (i in 1:length(dirs)) {
  h5.list[[i]] <- Read10X_h5(filename = dirs[i] )
  names(h5.list)[[i]] <- strsplit(dirname(dirs[i]), split = "/")[[1]][2]
}


```


> [!Tip] 
In R console type `?SeuratObject::CreateSeuratObject()` <br>
In the Help panel in RStudio, you will find a list of all available parameters that CreateSeuratObject uses. <br>
Scroll down to the Arguments section, please go through the description of each argument to understand what each argument/parameter does. <br>
You may use the `?` at any point in the R console followed by the name of the function to get help. e.g. ?RunPCA etc. 
![after](../images/m1.png)



Using the h5.list, we will now create seurat objects for each sample. 

```
# create seurat objects RNA assay ----
sobj.list = lapply(h5.list, FUN = function(x){
  
  x = CreateSeuratObject(counts = x, min.cells = 3, min.features = 200, assay = "RNA")
  
})

```

If we run the following code, you may notice that each cellbarcode (regardless of the sample it belongs to) is associated with a generic `SeuratProject` identifier.
```
sobj.list[[1]]$orig.ident %>% head()
sobj.list[[2]]$orig.ident %>% head()
sobj.list[[3]]$orig.ident %>% head()
```

![before](../images/i1.png)



We do not want this as if we were to merge these seurat objects together, we will loose the cellbarcode --> sample linkage and therefore will not be able to identify barcodes that belong to a given sample/condition. <br>
Therefore we want each cellbarcode to be associated with the correct sample it came from.  To correct this we will use the following helper loop to fix the issue. 

```
# by default SampleName which is represented as the orig.ident metadata variable in a seurat object will be named to 'SeuratProject', 
# we will use the following loop to overwrite that and rename with the sampleID

for (i in 1:length(sobj.list)) {
  sobj.list[[i]]$orig.ident <- names(sobj.list)[i]  
}


```

Run the same code chunk again and see how the results have changed.

```
sobj.list[[1]]$orig.ident %>% head()
sobj.list[[2]]$orig.ident %>% head()
sobj.list[[3]]$orig.ident %>% head()
```

![after](../images/i2.png)

<br>


Let us add some other useful metadata to each seurat object and merge all 3 seurat objects into 1 object that we will perform all of the downstream analysis on. These metrics will help us filter out noisey cells. 


```
# add % Mito ----
sobj.list = lapply(sobj.list, function(x){
  AddMetaData(object = x, metadata = PercentageFeatureSet(x, "^MT-"), col.name = "percent.mt")
})

# add log10GenesPerUMI ----
sobj.list = lapply(sobj.list, function(x){
  AddMetaData(object = x, metadata = log10(x$nFeature_RNA) / log10(x$nCount_RNA), col.name = "log10GenesPerUMI")
})
```

## 1.1 Merge Seurat Object

```
sobj <- merge(x = sobj.list[[1]], y = sobj.list[2:length(sobj.list)], merge.data=TRUE)

```

Save the merged seurat object for future access. 

```
saveRDS(sobj, "01_sobj.merged.RDS")
```


# 2. Initial QC

>[!Note]
Any saved RDS objects can be re-imported in R using the `readRDS` function.
e.g. to reload previously saved seurat object one may use the following command <br> `sobj <- readRDS("01_sobj.merged.RDS")`


We will first plot some violin plots and assess how the data looks. 

To plot a violin plot using seurat use the following code chunk:

```

VlnPlot(sobj, features = c("nCount_RNA", "nFeature_RNA" , "log10GenesPerUMI", "percent.mt"), pt.size = 0.1, group.by = "orig.ident", ncol = 4)

```

If you do not want to display points/dots on the violin plot change `pt.size = 0`




<hr>

<details>
  <summary> For Advanced Users </summary>

</details>
