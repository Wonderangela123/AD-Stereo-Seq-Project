library(Seurat)
library(SeuratDisk)
library(hdf5r)

## integrate samples
samples = list("GSM3704357_1-MAP2_filtered_feature_bc_matrix.h5", 
          "GSM3704359_2-MAP2_filtered_feature_bc_matrix.h5", 
          "GSM3704361_3-MAP2_filtered_feature_bc_matrix.h5", 
          "GSM3704363_4-MAP2_filtered_feature_bc_matrix.h5", 
          "GSM3704365_5-MAP2_filtered_feature_bc_matrix.h5", 
          "GSM3704367_6-MAP2_filtered_feature_bc_matrix.h5", 
          "GSM3704369_7-MAP2_filtered_feature_bc_matrix.h5",
          "GSM3704371_8-MAP2_filtered_feature_bc_matrix.h5",
          "GSM6261344_Control-1-MAP2_filtered_feature_bc_matrix.h5",
          "GSM6261345_Control-2-MAP2_filtered_feature_bc_matrix.h5",
          "GSM6261346_Control-3-MAP2_filtered_feature_bc_matrix.h5",
          "GSM6261347_Control-4-MAP2_filtered_feature_bc_matrix.h5",
          "GSM6261348_Control-5-MAP2_filtered_feature_bc_matrix.h5",
          "GSM6261349_Control-6-MAP2_filtered_feature_bc_matrix.h5",
          "GSM6261350_Control-7-MAP2_filtered_feature_bc_matrix.h5",
          "GSM6261351_Control-8-MAP2_filtered_feature_bc_matrix.h5")

path_dir = "/work/aliu10/stereo_project/reference/"

sample_list = list()
for (sample in samples){
    data = Read10X_h5(file = paste0(path_dir, sample))
    object = CreateSeuratObject(data, project = "SeuratProject", assay = "RNA")
    object$sample = sample
    object$diagnosis = ifelse(sample %in% samples[1:8], "case", "control") # case: tau-containing neuron; control: non tau-containing neuron
    sample_list[[length(sample_list) + 1]] <- object
}

# normalize and identify variable features for each dataset independently
neuron.list <- lapply(X = sample_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = neuron.list)
neuron.list <- lapply(X = neuron.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# perform integration (rpca)
anchors <- FindIntegrationAnchors(object.list = neuron.list, anchor.features = features, reduction = "rpca")
combined <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

# # Run the standard workflow for visualization and clustering
# combined <- ScaleData(combined, verbose = FALSE)
# combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
# combined <- FindClusters(combined, resolution = 0.5)

# save as anndata
setwd("/work/aliu10/stereo_project/")
SaveH5Seurat(combined, filename = "integrated_ref.h5Seurat", overwrite = TRUE)
Convert("integrated_ref.h5Seurat", dest = "h5ad", overwrite = TRUE)


library(zellkonverter)
library(scuttle)
library(SingleR)
library(Seurat)
library(SeuratDisk)

test = readH5AD(file = '/work/aliu10/stereo_project/integrated_ref.h5ad')
ref = readH5AD(file = '/work/ygong/stereo_seq_public/MIT_AD.h5ad')

# Convert unwanted levels to NA
ref$broad.cell.type[ref$broad.cell.type %in% c("Mic", "Per", "End")] <- NA
# Drop unused levels
ref$broad.cell.type <- droplevels(as.factor(ref$broad.cell.type))

assays(test)$logcounts <- assays(test)$X
assays(ref)$counts <- assays(ref)$X
ref=logNormCounts(ref)

## Extract neurons
# Obtain annotation dictionary
# labels: cell types; assay.type.test/ref: An integer scalar or string specifying the assay of test/ref containing the relevant expression matrix.
anno_dict_cluster = SingleR(test, ref, clusters=test$seurat_clusters, labels=ref$broad.cell.type, assay.type.test="logcounts", assay.type.ref="logcounts")

# the numbers in test$seurat_clusters are 0-based index (i.e., 0, 1, 2,...,24)
test$anno_cluster <- anno_dict_cluster$labels[test$seurat_clusters + 1]

subset_test <- test[, colData(test)$anno_cluster == 'Ex']

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = assays(subset_test)$logcounts, meta.data = as.data.frame(colData(subset_test)))

SaveH5Seurat(seurat_obj, filename = "reference.h5Seurat", overwrite = TRUE)
Convert("reference.h5Seurat", dest = "h5ad", overwrite = TRUE)
