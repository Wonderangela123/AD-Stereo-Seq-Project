library(Seurat)
library(Matrix)
library(arrow)
library(scuttle)
library(SingleR)
library(spatialLIBD)
library(dplyr)
library(future)
# install the Visium prefrontal cortext data
spe <- fetch_data(type = "spe")

# extract the row data
raw_data = spe@assays@data$counts

# get the gene name
# the gene name is the reference H5 file download from the 
rownames(raw_data) = rownames(Read10X_h5("../spatial_prefrontal_ref_data/151507_filtered_feature_bc_matrix.h5"))

# get the layer info for different spot
layer_info = spe@colData$spatialLIBD

# create the seurat object for the ref data
spatial_ref = CreateSeuratObject(counts = raw_data, min.cells=0, min.features=0)

# add the layer label to each cell
spatial_ref$label = layer_info

# import the stereopy data
data_path = "../temp_result/exp_matrix.csv"
test_sample = read_feather(data_path)
gene_list = as.vector(unlist(test_sample[, "__index_level_0__"]))
test_sample = select(test_sample, -c("__index_level_0__"))
rownames(test_sample) = gene_list
spatial_data = CreateSeuratObject(counts = test_sample, min.cells=0, min.features=0)    

# perform the annotation
spatial_ref_sce = as.SingleCellExperiment(spatial_ref)
spatial_data_sce = as.SingleCellExperiment(spatial_data)

spatial_ref_sce = logNormCounts(spatial_ref_sce)
spatial_data_sce = logNormCounts(spatial_data_sce)

pred.grun <- SingleR(
    test=spatial_data_sce, 
    ref=spatial_ref_sce, 
    labels=spatial_ref_sce$label, 
    de.method="wilcox",
    num.threads = 20
)

# save the result to the temp file
annotation_result = data.frame(cbind(pred.grun$labels, pred.grun$pruned.labels))
rownames(annotation_result) = rownames(pred.grun)
write.csv(annotation_result, "../temp_result/annotation_result.csv")
