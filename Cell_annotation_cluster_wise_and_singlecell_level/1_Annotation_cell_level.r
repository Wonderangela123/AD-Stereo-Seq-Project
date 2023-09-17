library(SingleR)
library(zellkonverter) 
library(scuttle)

## read h5ad files into SingleCellExperiment objects
test = readH5AD(file = 'Result/B01809C2.anndata.h5ad')
ref = readH5AD(file = 'Result/BA10_reference1.h5ad')

anno_dict_cell = SingleR(test, ref, labels=ref$broad.cell.type, assay.type.test = 1, assay.type.ref = 1)
saveRDS(anno_dict_cell, "Result/cell_level_anno.RDS")

# save the result
annotation = readRDS("Result/cell_level_anno.RDS")
write.csv(annotation, "Result/cell_level_anno.csv")

