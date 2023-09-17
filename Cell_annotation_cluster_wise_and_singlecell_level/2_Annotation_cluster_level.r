library(SingleR)
library(zellkonverter) 
library(scuttle)

## read h5ad files into SingleCellExperiment objects
test = readH5AD(file = 'Result/B01809C2.anndata_with_clustering.h5ad')
ref = readH5AD(file = 'Result/BA10_reference1.h5ad')

anno_dict_cell = SingleR(
    test, 
    ref, 
    labels=ref$broad.cell.type, 
    assay.type.test = 1, 
    assay.type.ref = 1,
    clusters = test$leiden
    
)
saveRDS(anno_dict_cell, "Result/cell_level_anno_cluster_wise.RDS")

# cluster wise result read
result = readRDS("Result/cell_level_anno_cluster_wise.RDS")
write.csv(result, "Result/cell_level_anno_cluster_wise.csv")

# read the cluster info
label = read.csv("Result/cell_level_anno_cluster_wise.csv")
label = data.frame(cbind(rownames(label), label$pruned.labels))
colnames(label) = c("group", "label")
label$group = as.character(label$group)

leiden_cluster = data.frame(test$leiden)
leiden_cluster$bin = colnames(test) 

combined = merge(leiden_cluster, label, by.x = "test.leiden", by.y = "group")
write.csv(combined, "Result/cell_level_anno_cluster_wise_combined.csv")


