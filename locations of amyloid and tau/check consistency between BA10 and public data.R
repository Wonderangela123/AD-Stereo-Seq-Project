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

########################
Scaling features for provided objects

Finding all pairwise anchors

Running CCA

Merging objects

Finding neighborhoods

Finding anchors

	Found 32521 anchors

Filtering anchors

	Retained 731 anchors

Merging dataset 1 into 2

Extracting anchors for merged samples

Finding integration vectors

Finding integration vector weights

Integrating data
