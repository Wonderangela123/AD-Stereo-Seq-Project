library(Seurat)
library(SeuratDisk)

setwd("/work/aliu10/AD_Stereoseq_Project/processed/")
samples_list = c("A02092E1", "B01809C2", "B02008C6", "B02008D2", "B02009F6", "C02248B5",
                 "D02175A4", "D02175A6", "B01809A3", "B01809A4", "B01806B5", "B01806B6")

# merge cases
cases_list = list()
for (sample in samples_list[1:6]){
    print(sample)
    Convert(paste0(sample, "/GeneExpMatrix/", sample, ".h5ad"), dest = "h5seurat", overwrite = TRUE)
    data = LoadH5Seurat(paste0(sample, "/GeneExpMatrix/", sample,".h5seurat"))
    cases_list[[length(cases_list) + 1]] <- data
}

merged_cases = Reduce(merge, cases_list)

SaveH5Seurat(merged_cases, filename = "cases.h5Seurat", overwrite = TRUE)
Convert("cases.h5Seurat", dest = "h5ad", overwrite = TRUE)

# merge controls
controls_list = list()
for (sample in samples_list[7:12]){
    print(sample)
    Convert(paste0(sample, "/GeneExpMatrix/", sample, ".h5ad"), dest = "h5seurat", overwrite = TRUE)
    data = LoadH5Seurat(paste0(sample, "/GeneExpMatrix/", sample,".h5seurat"))
    controls_list[[length(controls_list) + 1]] <- data
}

merged_controls = Reduce(merge, controls_list)

SaveH5Seurat(merged_controls, filename = "controls.h5Seurat", overwrite = TRUE)
Convert("controls.h5Seurat", dest = "h5ad", overwrite = TRUE)
