> [!Tip]
Useful Links:
[Seurat Version 5 Command Cheat Sheet](https://satijalab.org/seurat/articles/essential_commands)


### Cellranger Output

This week we start with the cellranger outputs that we created in the previous exercise. <br>
Conveniently, the output folders from cellranger are also available at the following path: `/path/to/cellranger-outs/`

It is important to note that **cellranger** creates a new directory for each sample and there are many files/subfolders that are created within that directory.
The count matrices produced by cellranger are located in the `./sampleName/outs/` directory. <br> For Seurat analysis we will be using the file named `filtered_feature_bc_matrix.h5` for each sample.


### Create Seurat Object

In an Rstudio session create a new project.

Create a new R script using the File --> New File --> Rscript or Shift/Cmd|Ctrl/N <br> 
To get started, we will need the path where the 10x outputs are present. Copy and paste the following line to specify the directory of cellranger outs:

```dirs = list.files(path = "10X_Data", full.names = T, pattern = ".h5$", recursive = T)```


Use the following code chunk to load all the R libraries that are needed. 

```
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(dplyr)

```










<details>
  <summary> For Advanced Users </summary>

</details>
