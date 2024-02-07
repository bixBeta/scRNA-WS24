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


### 1. Create Seurat Object

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
Here you will find a list of all available parameters that the CreateSeuratObject uses. In the Arguments section, please go through the description to understand what each argument/parameter means. 
You may use the `?` at any point in the R console followed by the name of the function to get help. e.g. ?RunPCA etc. 


<hr>

<details>
  <summary> For Advanced Users </summary>

</details>
