library(Seurat)
library(SeuratDisk)

setwd("/work/aliu10/AD_Stereoseq_Project/reference")
data_pub <- LoadH5Seurat("integrated.h5Seurat")
data_pub$region = "BA9"

Convert('/work/ygong/stereo_seq_public/MIT_AD.h5ad', dest = "h5seurat", overwrite = TRUE)
data_MIT = LoadH5Seurat('/work/ygong/stereo_seq_public/MIT_AD.h5seurat')
data_MIT$region = "BA10"

data_MIT_neuron <- subset(data_MIT, subset = broad.cell.type %in% c("Ex", "In"))

x <- NormalizeData(data_MIT_neuron)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)

# identify variable features for each dataset independently
my_list = list(x, data_pub)
features <- SelectIntegrationFeatures(object.list = my_list)

anchors <- FindIntegrationAnchors(object.list = my_list, anchor.features = features)
# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

png("/work/aliu10/AD_Stereoseq_Project/reference/BA9vs10.png")
DimPlot(combined, reduction = "umap", group.by = "region", raster=FALSE)
dev.off() 
