# Stereo-seq analysis pipeline
# Step one: Convert the stereoseq gene expression data into Seurat to perform the multiple sample integration

""" Run in python """
# Read the Stereoseq data through python and write the gene expression matrix into feather file and gene name into csv 
import stereo as st
import utils
from stereo.core.ms_data import MSData
from stereo.core.ms_pipeline import slice_generator
import pandas as pd
import pyarrow.feather as feather

samples_list = ["B01806B5", "B01806B6", "B01809A3", "B01809A4", "B01809C2", "D02175A4", "D02175A6"]
for sample in samples_list:
    print(sample)
    data_path = sample+ "/GeneExpMatrix/" + sample + ".cellbin.gef"
    st.io.read_gef_info(data_path)
    data1 = st.io.read_gef(file_path=data_path, bin_type='cell_bins')
    df = data1.to_df()
    df = df.T
    feather.write_feather(df, sample + "/matrix.feather")
    gene_list = pd.DataFrame(df.index)
    gene_list.to_csv(sample + "/gene_list.csv")
""""""

# After writing the gene expression and gene label for each individual, use the Seurat in R to read the data
# High memory needed and time consuming (finished on loni)

""" Run in R """
library(Seurat)
library(zellkonverter)
library(Matrix)

# clean the environment
rm(list = ls())
gc()

sample_list = c("B01806B5", "B01806B6", "B01809A3", "B01809A4", "B01809C2", "D02175A4", "D02175A6")
for (sample in sample_list){
    print(sample)
    # import gene expression matrix
    library(arrow)
    data = read_feather(paste0(sample,"/matrix.feather"))
    
    gene_list = read.csv(paste0(sample,"/gene_list.csv"))
    gene_list[,2]
    
    rownames(data) = gene_list[,2]
    
    # Create Seurat data object
    pbmc <- CreateSeuratObject(
      counts = data, 
      project = "pbmc3k", 
      min.cells = 0, 
      min.features = 0
      )

    # Find the mitochondria genes, cells with high mitochondira gene expression are considered as dead cells    
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    
    
    # Add the label for each sample
    pbmc$sample = sample
    pbmc$diagnosis = "Normal"
    
    # save the R RDS file for each sample    
    saveRDS(pbmc, paste0(sample, "/Seurat_gene.RDS"))
    rm(list=ls())
    gc()
}
""""""

# Perform the integration using sctransform for normalization and CCA for integration (working on loni)
""" Run in R """
library(Seurat)
library(tidyr)

# read the RDS file from each sample folder
data_list = c("B01806B5", "B01806B6", "B01809A3", "B01809A4", "B01809C2", "D02175A4", "D02175A6")
combine_data_list = list()

for (i in data_list){
    file = readRDS(paste0(i, "/Seurat_gene.RDS"))
    matrix = file@assays$RNA@counts
    matrix<-matrix[,!grepl("__index_level_0__",colnames(matrix))] # remove the extra columns
    pbmc <- CreateSeuratObject(counts = matrix, project = "pbmc3k", min.cells = 3, min.features = 200)
    pbmc$sample = i
    combine_data_list = c(combine_data_list, pbmc)
}

ifnb.list <- lapply(X = combine_data_list,FUN = SCTransform)


features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)

# save the data as RDS in the Result file
saveRDS(immune.combined.sct, "Result/immune.combined.RDS")
""""""









































