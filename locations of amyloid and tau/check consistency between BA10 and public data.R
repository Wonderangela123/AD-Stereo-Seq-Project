library(Seurat)
library(SeuratDisk)

setwd("/work/aliu10/AD_Stereoseq_Project/reference")
data_pub <- LoadH5Seurat("integrated.h5Seurat")
data_pub$region = "BA9"

setwd("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/")
Convert("B01809C2_extract.anndata.h5ad", dest = "h5seurat", overwrite = TRUE)
data = LoadH5Seurat("B01809C2_extract.anndata.h5seurat")
data$region = "BA10"

# identify variable features for each dataset independently
my_list = list(data, data_pub)
features <- SelectIntegrationFeatures(object.list = my_list)

anchors <- FindIntegrationAnchors(object.list = my_list, anchor.features = features)
# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

p1 <- DimPlot(combined, reduction = "umap", group.by = "BA10")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
