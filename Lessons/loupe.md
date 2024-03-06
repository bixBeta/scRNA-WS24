
Loupe Browser 7 allows a user to import seurat objects into its browser framework. <br>
You can use majority of features of the Loupe Browser in the Loupe browser GUI. It can be a great resource to visualize marker genes/ custom gene lists in an interactive way. 

<hr>


### Getting Started: 

To get started, please install the loupeR package using the following code chunk:

```
remotes::install_github("10XGenomics/loupeR")
loupeR::setup()
```

Once installed, please load the loupeR library and follow the following steps. 

```
# import the library
library("loupeR")
```

First we will extract the assay of interest and save it to a variable

```
# Gene Expression RNA assay, if working with SCT, please change the "RNA" to "SCT"
assay <- seurat_obj[["RNA"]]
```

Next we isolate the count data from the assay object that we just created. 

```
# get counts matrix from either the old or newer formats of assay
counts <- counts_matrix_from_assay(assay)
```

Lastly, we will use the `create_loupe` function to create a loupe file named `mySeuratObject.cloupe`

```
# convert the count matrix, clusters, and projections into a Loupe file
create_loupe(
    counts,
    clusters = select_clusters(seurat_obj),
    projections = select_projections(seurat_obj),
    output_name = "mySeuratObject"
)
```

Once finished, open up the loupe browser and upload the `mySeuratObject.cloupe` file. 

Happy Exploring!!!
