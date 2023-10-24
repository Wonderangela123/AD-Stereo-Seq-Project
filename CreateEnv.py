## Python environment
conda create -n st python=3.8
conda activate st

## Development version
git clone -b dev https://github.com/STOmics/stereopy.git
cd stereopy
python setup.py install

conda install -c conda-forge notebook
conda install -c conda-forge scanpy
conda install -c conda-forge anndata


## R environment
conda create -n st-R
conda activate st-R

conda install -c conda-forge r-base
conda install -c conda-forge notebook
conda install -c conda-forge r-irkernel  ## use R in notebook

conda install -c conda-forge r-seurat

## install Bioconda packages in R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")

## install "tar.gz" in linux 
## download from https://bioconductor.org/packages/release/bioc/src/contrib/zellkonverter_1.10.1.tar.gz
R CMD INSTALL zellkonverter_1.10.1.tar.gz
