> [!Tip]
Useful Links:
[Seurat Version 5 Command Cheat Sheet](https://satijalab.org/seurat/articles/essential_commands)


### Cellranger Output

This week we start with the cellranger outputs that we created in the previous exercise. <br>
Conveniently, the output folders from cellranger are also available at the following path: `/path/to/cellranger-outs/`

It is important to note that **cellranger** creates a new directory for each sample and there are many files/subfolders that are created within that directory.
The count matrices produced by cellranger are located in the `./sampleName/outs/` directory. <br> For Seurat analysis we will be using the file named `filtered_feature_bc_matrix.h5` for each sample.














<details>
  <summary> For Advanced Users </summary>

</details>
