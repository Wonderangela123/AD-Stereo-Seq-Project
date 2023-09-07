library(SingleR)
library(zellkonverter) 

## read h5ad files into SingleCellExperiment objects
test = readH5AD(file = '/work/aliu10/AD_Stereoseq_Project/processed/data/integrated_subtype.anndata.h5ad')
ref = readH5AD(file = '/work/aliu10/AD_Stereoseq_Project/reference/reference.h5ad')

ref$anno_cluster_leiden <- as.character(ref$anno_cluster_leiden)
# Get a logical vector indicating which cells are "Ex" or "In" types
selected_cells <- grepl("Ex|In", ref$anno_cluster_leiden)
# Subset the SCE object using the logical vector
ref <- ref[, selected_cells]

## Obtain annotation dictionary
## labels: cell types; assay.type.test/ref: An integer scalar or string specifying the assay of test/ref containing the relevant expression matrix.
anno_dict_cluster = SingleR(test, ref, clusters=test$leiden, labels=ref$diagnosis, assay.type.test=1, assay.type.ref=1) ## cluster-level annotation

## The function "data.tl.annotation" requires "Categorical categories must be unique", hence, create a variable to paste cluster and labels together.
## cluster annotation
anno_dict_cluster$cluster = row.names(anno_dict_cluster)
anno_dict_cluster$dict = paste(anno_dict_cluster$labels, anno_dict_cluster$cluster, sep = ".")

write.table(anno_dict_cluster, file = '/work/aliu10/AD_Stereoseq_Project/processed/data/annotation_dict_subtype.txt')
