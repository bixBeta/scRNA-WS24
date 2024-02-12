### Week 1. Single cell technology, from samples to raw data to count matrix

#### 1. Connect to the server

Find your assigned server name on [this website](https://biohpc.cornell.edu/ww/machines.aspx?i=165). Instructions to connect to the server can be found on [this page](https://biohpc.cornell.edu/lab/doc/remote_access.pdf).

#### 2. Set up your working directory and copy data files for this exercise.

```
mkdir /workdir/$USER

cd /workdir/$USER

cp -r /workdir/sc_workshop_2024/cellranger .
```

<details>
<summary>Expand for Linux details</summary>

##### What do those commands mean?
+ `mkdir` make directory
+ `$USER` replace with your netid (i.e. BioHPC username)
+ `cd` change directory
+ `cp -r` copy directory and all of its contents to `.` (dot = here, your current location)

##### A few more examples:

+ Where am I?

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`pwd`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;= path of the present working directory


+ What files/directories are here?

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ls -l`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;= list (long format)

&nbsp;&nbsp;&nbsp;&nbsp;<i>you should now see the 'cellranger' directory copied to your working directory</i>

##### We will use a lot of 'paths' describing the location of files
+ /workdir/$USER: is the path to the working directory you set up for this workshop on your assigned server
+ For many commands, you can use the full path (starting from the root /) or a relative path, but either way it needs to be correct

##### Pay close attention to syntax and characters!
+ forward-slashes denote nested directory and file names
+ back-slashes have a very different meaning in Linux, they define 'escape' sequences to override literal character meanings. A common use of backslashes in this workshop will be to 'escape' the end-of-line character when code breaks across several lines. If the code is on a single line, the backslashes are not needed.
</details>



#### 3. Run cellranger to process single cell RNAseq data

As it could take several hours to run cellranger on the original data files, special fastq files were prepared with 1/2000 downsampling.  Within each sub-directory, you will find a README file with a description of the sample and a link to the original data.

Under the data directory you just copied, there are 6 sub-directories containing fastq files:

1. IgG1d:	 **Singleplex** single cell gene expression library from GEO series GSE201999.
2. IgG4:           **Singleplex** single cell gene expression library from GEO series GSE201999.
3. UT:              **Singleplex** single cell gene expression library from GEO series GSE201999.
4. cellplex:      **Multiplex** single cell gene expression library.
5. cellplex_fb: **Multiplex** single cell gene expression library with **antibody capture** 
6. refdata-gex-GRCh38-2020-A: **human reference genome required for cellranger**

<details>
	<summary>Expand for file details</summary>
    
##### GEO series GSE201999
The first 3 directories contain data from an experiment described in this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10580117/) and will comprise the main dataset used in this workshop.

##### Multiplex datasets
Two additional datasets are provided as examples to run 'cellranger multi' on datasets with hashed (multiplexed) samples and feature barocoding (CITE-seq, antibody capture).

##### Human reference genome
This directory was downloaded from 10x Genomics as a pre-indexed reference genome (+transcriptome) in the format required by cellranger. 10x Genomics provides a limited set of pre-indexed genomes; additional indexed genomes can be created with 'cellranger mkref' for gene expression (GEX) analysis of other species.

##### What files are located in these directories?
+ Try listing the contents of each directory
+ Add a new parameter to list all of the contents of all subdirectories in one command:
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ls -lR`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;R = recursive list
</details>

##### 3.1 Run "cellranger count" on a singleplex data set.


```
cd /workdir/$USER/cellranger

export PATH=/programs/cellranger-7.2.0:$PATH

cellranger count --id=run_IgG1d --sample=IgG1d --transcriptome=refdata-gex-GRCh38-2020-A --fastqs=IgG1d  --localcores=8 --localmem=24

```


<details>
 <summary>Expand for Linux details</summary>
    
##### What do these commands mean?
+ cd: change directory
    + You were at the top level of /workdir/$USER/, and moved into the subdirectory named 'cellranger'
    + You need to be in the 'cellranger' subdirectory for the code to work.
    + If you are not in the right place, figure out where you are ('pwd') and move ('cd') to /workdir/$USER/cellranger

+ export PATH: this makes sure your code can find the cellranger program easily

+ cellranger count: starts the 'cellranger count' program, with the following parameters:
    + `--id=` &nbsp;&nbsp;&nbsp;&nbsp;name of the cellranger run <b>required</b>
    + `--sample=` &nbsp;&nbsp;&nbsp;&nbsp;name of the sample to analyze (must match the beginning of the fastq filename) <b>required</b>
    + `--transcriptome=` &nbsp;&nbsp;&nbsp;&nbsp;name of the directory containing the reference genome (index formatted for cellranger) <b>required</b>
    + `--fastqs=` &nbsp;&nbsp;&nbsp;&nbsp;location of the directory containing the fastq files (which happens to be the same as the sample name) <b>required</b>
    + `--localcores=` &nbsp;&nbsp;&nbsp;&nbsp;number of CPUs (cores) that your run of cellranger can use <b>recommended for shared servers</b>
    + `--localmem=` &nbsp;&nbsp;&nbsp;&nbsp;amount of memory that your run of cellranger can use <b>recommended for shared servers</b>

</details>


Notes:

- This step takes about 5-10 minutes. While you are waiting, you can do Step 4: Run Loupe browser;

- The output files are in the directory run_IgG1d/out. You can browse the output directory. Files of interest include : 

  - web_summary.html: a summery web page

  - possorted_genome_bam.bam: position-sorted read alignment;

  - filtered_feature_bc_matrix.h5: HDF5 formatted single gene expression count matrix;

  - filtered_feature_bc_matrix: directory containing MEX formatted matrix files. 

    Either filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix can be used for downstream data analysis with Seurat or Scanpy.

- The parameters "--localcores=8 --localmem=24" restrict the memory cpu core usage by the job. 

- Cellranger results from the original data of the same library are provided in the directory "/workdir/sc_workshop_2024/GSE201999_output". In Step 4, you will examine the full dataset (without downsampling). 



##### 3.2 Run "cellranger multi" on a multiplex data set

3.2.1 Inspect and modify the cellplex/config.ori.csv file

Replace all "xxxxx" with your BioHPC userID, and write the modified content to a new file config.csv.   You can use the LINUX sed command to do this.

```
cd /workdir/$USER/cellranger

sed "s/xxxxx/$USER/" cellplex/config.ori.csv > cellplex/config.csv

cat cellplex/config.csv
```

Or use a text editor such as `nano`.

3.2.2 run "cellranger multi"

```

cellranger multi --id=run_plex --csv=cellplex/config.csv --localcores=8 --localmem=40
```

<details>
<summary>Expand for Linux details</summary>

##### What do these commands mean?
+ cellranger multi: starts the 'cellranger multi' program, with the following parameters:
    + `--id=` &nbsp;&nbsp;&nbsp;&nbsp;name of the cellranger run <b>required</b>
    + `--csv=` &nbsp;&nbsp;&nbsp;&nbsp;name of the config file for the run <b>required</b>
    + `--localcores=` &nbsp;&nbsp;&nbsp;&nbsp;number of CPUs (cores) that your run of cellranger can use <b>recommended for shared servers</b>
    + `--localmem=` &nbsp;&nbsp;&nbsp;&nbsp;amount of memory that your run of cellranger can use <b>recommended for shared servers</b>

##### How is running 'cellranger multi' different than 'cellranger count'?
+ what information is in the config file?

</details>

Notes:

- This step takes about 10-15 minutes. 
- The output directory is located in run_plex/outs. You would find a "per_sample_outs" directory, with results from each de-multiplexed sample.



##### 3.3 Run "cellranger count" on a multiplex data set with antibody capture

3.3.1 Inspect and modify the cellplex_fb/config.ori.csv file

Replace all "xxxxx" with your BioHPC userID, and write the modified content to a new file config.csv.   You can use the LINUX sed command to do this.

```
cd /workdir/$USER/cellranger

sed "s/xxxxx/$USER/" cellplex_fb/config.ori.csv > cellplex_fb/config.csv

cat cellplex_fb/config.csv
```

3.3.2 run "cellranger multi"

```

cellranger multi --id=run_plex_fb --csv=cellplex_fb/config.csv --localcores=8 --localmem=40
```

Notes:

- This step takes about 10-15 minutes. 
- The output directory is located in run_plex_fb/outs. Similar to the previous run, you would find a "per_sample_outs" directory, with results from each demultiplexed sample. 
- The protein expression levels are  included in the gene expression count matrix, as additional features. 

#### 4. Review 'web_summary' QC and run Loupe browser

##### 4.1 Download the following directory to your laptop or computer

User Filezilla (or other sftp client) to download cellranger results from the original data.

**Server**: your assigned server

**Directory**: /workdir/sc_workshop_2024/GSE201999_output

**Port:** 22

<b>Transfer the full directory with Filezilla to your laptop or computer</b>

###### There are three sub-directories, each from a different sample: 
+ TL-IgG1
+ TL-IgG4
+ Untreated

Within each sample directory, there should be a "web_summary.html"  file and a "cloupe.cloupe" file.

##### 4.2 Open the web_summary.html files in a web browser;

Documentation of the web summary file and QC metrics:

https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-outputs-web-summary-count
https://www.10xgenomics.com/analysis-guides/quality-assessment-using-the-cell-ranger-web-summary

What do the QC metrics indicate about these samples?

##### 4.3. Start Loupe, and open the file cloupe.cloupe;

 Tutorials for Loupe:

https://www.10xgenomics.com/support/software/loupe-browser/latest/tutorials/introduction/lb-sc-interface-and-navigation

https://www.10xgenomics.com/support/software/loupe-browser/latest/tutorials/assay-analysis/lb-single-cell-gene-expression-and-flex


#### 5. Things to take into consideration when working with real data

- It could take a few hours to process each set of real data files. You will need to run cellranger with `nohup ` or in a persistent "screen" session;

- Most likely you would need to run cellranger on multiple samples. Some of the BioHPC servers (large memory gen2) have >100 cpu cores, and there is no performance benefit to run cellranger on a single sample with >32 CPU cores.  You would want to run the jobs in parallel to use the all the available CPU cores.
 
**You can do this step later during the week.** 

##### 5.1 Create a batch script file using a text editor (e.g. Notepad++ on Windows, or BBEdit on Mac, or nano on Linux), with the following content. Upload to the /workdir/$USER/cellranger directory of your assigned server


```
cellranger count --id=run_IgG1d --sample=IgG1d --transcriptome=/workdir/qisun/cellranger/refdata-gex-GRCh38-2020-A --fastqs=/workdir/qisun/cellranger/IgG1d  --localcores=8 --localmem=24

cellranger count --id=run_IgG4 --sample=IgG4 --transcriptome=/workdir/qisun/cellranger/refdata-gex-GRCh38-2020-A --fastqs=/workdir/qisun/cellranger/IgG4  --localcores=8 --localmem=24

cellranger count --id=run_UT --sample=UT --transcriptome=/workdir/qisun/cellranger/refdata-gex-GRCh38-2020-A --fastqs=/workdir/qisun/cellranger/UT  --localcores=8 --localmem=24
```

###### Notes:
+ This code uses the full path for the reference index (`--transcriptome`) and fastq location (`--fastqs=`), allowing the script to be run from any directory.
+ The output directories (named as `--id=`) will be created in the directory that you run the script.


##### 5.2 Run the script in "screen" persistent session.

Tutorial for "screen" can be found [here](https://biohpc.cornell.edu/lab/doc/Linux_exercise_part2.pdf).


There are two ways to run this script:

Run in serial (run cellranger on one sample at a time):

```
export PATH=/programs/cellranger-7.2.0:$PATH
sh run.sh
```

Run in parallel ("-j 2" : process two samples  at a time):

```
export PATH=/programs/cellranger-7.2.0:$PATH
parallel -j 2 < run.sh
```

Notes: 

- When running in parallel, make sure that the numbers set in "--localmem --localcores" multiply by number of jobs do not exceed the total RAM or CPU cores on the server you are using. Your assigned server for this workshop has 128GB RAM and 24 CPU cores. 
