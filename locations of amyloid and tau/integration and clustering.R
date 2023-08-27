library(Seurat)
library(SeuratDisk)

# read h5 files for each sample and merge (cases)
cases_list = list("GSM3704357_1-MAP2_filtered_feature_bc_matrix.h5", 
                  "GSM3704359_2-MAP2_filtered_feature_bc_matrix.h5", 
                  "GSM3704361_3-MAP2_filtered_feature_bc_matrix.h5", 
                  "GSM3704363_4-MAP2_filtered_feature_bc_matrix.h5", 
                  "GSM3704365_5-MAP2_filtered_feature_bc_matrix.h5", 
                  "GSM3704367_6-MAP2_filtered_feature_bc_matrix.h5", 
                  "GSM3704369_7-MAP2_filtered_feature_bc_matrix.h5",
                  "GSM3704371_8-MAP2_filtered_feature_bc_matrix.h5")

my_list_cases = list()
for (sample in cases_list){
    print(sample)
    data = Read10X_h5(file = paste0("/work/aliu10/AD_Stereoseq_Project/reference/", sample))
    object = CreateSeuratObject(data, project = "SeuratProject", assay = "RNA")
    object$diagnosis = "case"
    my_list_cases[[length(my_list_cases) + 1]] <- object
}

merge_cases = Reduce(merge, my_list_cases)

# read h5 files for each sample and merge (controls)
controls_list = c("GSM6261344_Control-1-MAP2_filtered_feature_bc_matrix.h5",
                  "GSM6261345_Control-2-MAP2_filtered_feature_bc_matrix.h5",
                  "GSM6261346_Control-3-MAP2_filtered_feature_bc_matrix.h5",
                  "GSM6261347_Control-4-MAP2_filtered_feature_bc_matrix.h5",
                  "GSM6261348_Control-5-MAP2_filtered_feature_bc_matrix.h5",
                  "GSM6261349_Control-6-MAP2_filtered_feature_bc_matrix.h5",
                  "GSM6261350_Control-7-MAP2_filtered_feature_bc_matrix.h5",
                  "GSM6261351_Control-8-MAP2_filtered_feature_bc_matrix.h5")

my_list_controls = list()
for (sample in controls_list){
    print(sample)
    data = Read10X_h5(file = paste0("/work/aliu10/AD_Stereoseq_Project/reference/", sample))
    object = CreateSeuratObject(data, project = "SeuratProject", assay = "RNA")
    object$diagnosis = "control"
    my_list_controls[[length(my_list_controls) + 1]] <- object
}

merge_controls = Reduce(merge, my_list_controls)

# normalize and identify variable features for each dataset independently
my_list = list(merge_cases, merge_controls)
neuron.list <- lapply(X = my_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = neuron.list)
neuron.list <- lapply(X = neuron.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = neuron.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"

# Run the standard workflow for clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

# save as anndata
setwd("/work/aliu10/AD_Stereoseq_Project/reference/")
SaveH5Seurat(combined, filename = "clusters.h5Seurat", overwrite = TRUE)
Convert("clusters.h5Seurat", dest = "h5ad", overwrite = TRUE)
