library(SingleR)
library(zellkonverter) 
library(scater)

## read h5ad files into SingleCellExperiment objects
test = readH5AD(file = '/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.anndata_with_clustering.h5ad')
ref = readH5AD(file = '/work/ygong/stereo_seq_public/MIT_AD.h5ad')

assayNames(ref)[assayNames(ref) == "X"] <- "counts"
ref = logNormCounts(ref)

anno_dict_cluster = SingleR(test, ref, clusters=test$leiden, labels=ref$broad.cell.type, assay.type.test=1, assay.type.ref=2) ## cluster-level annotation

## The function "data.tl.annotation" requires "Categorical categories must be unique", hence, create a variable to paste cluster and labels together.
## cluster annotation
anno_dict_cluster$cluster = row.names(anno_dict_cluster)
anno_dict_cluster$dict = paste(anno_dict_cluster$labels, anno_dict_cluster$cluster, sep = ".")

write.table(anno_dict_cluster, file = '/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2_annotation_dict_cluster.txt')
