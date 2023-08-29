library(SingleR)
library(zellkonverter) 

## read h5ad files into SingleCellExperiment objects
test = readH5AD(file = '/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/B01809C2_recluster.anndata.h5ad')
ref = readH5AD(file = '/work/ygong/stereo_seq_public/MIT_AD.h5ad')

## Obtain annotation dictionary
## labels: cell types; assay.type.test/ref: An integer scalar or string specifying the assay of test/ref containing the relevant expression matrix.
anno_dict_cluster = SingleR(test, ref, clusters=test$leiden, labels=ref$broad.cell.type, assay.type.test=1, assay.type.ref=1) ## cluster-level annotation
anno_dict_cell = SingleR(test, ref, labels=ref$broad.cell.type, assay.type.test=1, assay.type.ref=1) ## cell-level annotation


## The function "data.tl.annotation" requires "Categorical categories must be unique", hence, create a variable to paste cluster and labels together.
## cluster annotation
anno_dict_cluster$cluster = row.names(anno_dict_cluster)
anno_dict_cluster$dict = paste(anno_dict_cluster$labels, anno_dict_cluster$cluster, sep = ".")

write.table(anno_dict_cluster, file = '/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/annotation_dict_cluster_recluster.txt')



# ## cell annotation
# anno_dict_cell$cell = row.names(anno_dict_cell)
# anno_dict_cell$dict = paste(anno_dict_cell$labels, anno_dict_cell$cell, sep = ".")

# write.table(anno_dict_cell, file = '/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/annotation_dict_cell_recluster.txt')
