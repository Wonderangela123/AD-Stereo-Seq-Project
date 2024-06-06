library(Matrix)
library(dplyr)
library(rjson)
library(Seurat)
library(ggplot2)
#library(SeuratDisk)
library(Matrix)
#library(SpaceX)
library(cowplot)
library(patchwork)
library(tidyr)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
#library(SeuratData)
library(igraph)
library(Matrix)
enableWGCNAThreads(nThreads = 48)

# read the matrix
matrix = readMM("../Result/Concentric_to_R/concentric_outer3_all.mtx")

# read the cell and gene names
meta = read.csv("../Result/Concentric_to_R/outer3_meta_data.csv")
feature_name = read.csv("../Result/Concentric_to_R/outer3_raw_feature.csv")

rownames(matrix) = meta$X
colnames(matrix) = feature_name$X

meta_info = read.csv("../Result/Concentric_to_R/outer3_meta_data.csv")

data = CreateSeuratObject(counts = t(matrix), project = "pbmc3k", min.cells = 0, min.features = 0)

# add the meta information
for (i in meta_info){
    for (i in colnames(meta_info)){
        data[[i]] = meta_info[[i]]
    }
}

meta_info = read.csv("../Result/Concentric_to_R/meta_data.csv")

meta_sub = meta_info[,c("barcode", "concentric")]

meta_all = data@meta.data
meta_all$concentric = NULL

meta_df <- merge(meta_all, meta_sub, by = "barcode", all.x = TRUE)

sub = subset(meta_df, annotation == "Ex" & levels == "severe" & concentric == "inner")

meta_df$concentric[is.na(meta_df$concentric)] <- "outer3"

for (i in meta_df){
    for (i in colnames(meta_df)){
        data[[i]] = meta_df[[i]]
    }
}

data <-  data %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

data <- SetupForWGCNA(
  data,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "vis"
)

data <- MetacellsByGroups(
  seurat_obj = data,
  group.by = c("annotation","sample_id"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'annotation' # set the Idents of the metacell seurat object
)

data <- NormalizeMetacells(data)

saveRDS(data, "../Result/WGCNA/WGCNA_inner_outer3_concentric.RDS")
##############################################################################################
data <- SetDatExpr(
  data,
  group_name = "Ex", # the name of the group of interest in the group.by column
  group.by='annotation', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data', # using normalized data
  use_metacells = TRUE
)

# Test different soft powers
data<- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

power_table <- GetPowerTable(data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
