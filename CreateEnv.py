## Python environment
conda create -n ADstereo python=3.8
conda activate ADstereo

## Development version
git clone -b dev https://github.com/STOmics/stereopy.git
cd stereopy
python setup.py install

conda install -c conda-forge notebook


## R environment
conda create -n ADstereo-R
conda activate ADstereo-R

conda install -c conda-forge r-base
conda install -c conda-forge notebook
conda install -c conda-forge r-irkernel  ## use R in notebook

conda install -c conda-forge r-seurat

## install "SingleR" in R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")
