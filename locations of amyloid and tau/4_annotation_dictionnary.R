library(SingleR)
library(zellkonverter) 

## read h5ad files into SingleCellExperiment objects
test = readH5AD(file = '/work/aliu10/AD_Stereoseq_Project/processed/data/integrated_subtype.anndata.h5ad')
ref = readH5AD(file = '/work/aliu10/AD_Stereoseq_Project/reference/reference.h5ad')

ref$anno_cluster_leiden <- as.character(ref$anno_cluster_leiden)
ref$anno_cluster_leiden[grepl("Ex", ref$anno_cluster_leiden)] <- "Ex"
ref$anno_cluster_leiden[grepl("In", ref$anno_cluster_leiden)] <- "In"
ref$anno_cluster_leiden[grepl("Opc", ref$anno_cluster_leiden)] <- "Opc"
ref$anno_cluster_leiden[grepl("Ast", ref$anno_cluster_leiden)] <- "Ast"
ref$anno_cluster_leiden[grepl("Oli", ref$anno_cluster_leiden)] <- "Oli"
ref$anno_cluster_leiden <- as.factor(ref$anno_cluster_leiden)

## Obtain annotation dictionary
## labels: cell types; assay.type.test/ref: An integer scalar or string specifying the assay of test/ref containing the relevant expression matrix.
anno_dict_cluster = SingleR(test, ref, clusters=test$leiden, labels=ref$anno_cluster_leiden, assay.type.test=1, assay.type.ref=1) ## cluster-level annotation
# anno_dict_cell = SingleR(test, ref, labels=ref$anno_cluster_leiden, assay.type.test=1, assay.type.ref=1) ## cell-level annotation


## The function "data.tl.annotation" requires "Categorical categories must be unique", hence, create a variable to paste cluster and labels together.
## cluster annotation
anno_dict_cluster$cluster = row.names(anno_dict_cluster)
anno_dict_cluster$dict = paste(anno_dict_cluster$labels, anno_dict_cluster$cluster, sep = ".")

write.table(anno_dict_cluster, file = '/work/aliu10/AD_Stereoseq_Project/processed/data/annotation_dict_subtype.txt')



# ## cell annotation
# anno_dict_cell$cell = row.names(anno_dict_cell)
# anno_dict_cell$dict = paste(anno_dict_cell$labels, anno_dict_cell$cell, sep = ".")

# write.table(anno_dict_cell, file = '/work/aliu10/AD_Stereoseq_Project/processed/data/annotation_dict_cell.txt')
