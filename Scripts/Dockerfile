# This is Dockerfile used to build scrna2024 image used in Week 4
# command to build is:

# docker1 build -f Dockerfile -t scrna2024 `pwd`

#If you add commands to the end of this file, you can call 'build', and it will add a layer to previously built image

# Once you are happy with an image, save it with docker1 save -o image.tar imageID

FROM ubuntu:22.04
MAINTAINER mjhubisz@cornell.edu
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN apt install -y python3 python3-pip vim git locate
    
RUN echo "deb [trusted=yes] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y r-base r-base-dev

    
RUN pip install pip==24.0
RUN pip install scanpy jupyterlab anndata2ri


# install seurat
RUN apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev libglpk-dev
RUN Rscript -e 'install.packages("Seurat")'

RUN Rscript -e 'if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}; BiocManager::install("MAST")'

RUN Rscript -e 'BiocManager::install("rhdf5")'

RUN Rscript -e 'install.packages("argparse")'

RUN Rscript -e 'install.packages("IRkernel"); IRkernel::installspec(name="ir432", displayname="R 4.3.2", user=FALSE)'

# monocle
RUN apt-get update &&  apt install -y cmake libgdal-dev libfontconfig1-dev libcairo2-dev libharfbuzz-dev libfribidi-dev libudunits2-dev
RUN Rscript -e 'BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "lme4", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "batchelor", "HDF5Array", "terra", "ggrastr"))'
RUN Rscript -e 'install.packages("devtools"); devtools::install_github("cole-trapnell-lab/monocle3")'

RUN Rscript -e 'install.packages("hdf5r")'
RUN Rscript -e 'install.packages("R.utils"); devtools::install_github("satijalab/seurat-wrappers"); BiocManager::install("glmGamPoi")'

RUN pip install iprogress
RUN pip install Cython
RUN pip install velocyto

# I removed these tools because they were huge (~8 Gb)
#RUN pip install scvi-tools
RUN pip install scikit-misc
RUN pip install louvain 
RUN pip install harmonypy
RUN pip install magic-impute pandas==2.2.1
RUN pip install scrublet

RUN mkdir -p /root/scripts/
COPY startJupyter.sh /root/scripts/
COPY overrides.json /usr/local/share/jupyter/lab/settings/
    
ENTRYPOINT ["bash", "/root/scripts/startJupyter.sh"]
CMD ["/workdir"]
ENV SHELL=/bin/bash

RUN apt install -y git rsync
RUN apt install -y samtools

