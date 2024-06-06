library(dplyr)
library(rjson)
library(Seurat)
library(ggplot2)
#library(SeuratDisk)
library(Matrix)
#library(SpaceX)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
#library(SeuratData)

enableWGCNAThreads(nThreads = 48)

# read the count data
matrix = readMM("../Result/R_conversion/raw_matrix.mtx")

gene_name = read.csv("../Result//R_conversion//gene_names.csv")
meta = read.csv("../Result/R_conversion/meta_info.csv")

rownames(matrix) = meta$X
colnames(matrix) = gene_name$X0

# create seurat object
merge = CreateSeuratObject(counts = t(matrix), project = "pbmc3k", min.cells = 0, min.features = 0)

for (i in colnames(meta)){
    merge[[i]] = meta[[i]]
}

merge <- merge %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

merge <- FindNeighbors(merge, dims = 1:30)
merge <- FindClusters(merge,verbose = TRUE)
#merge <- RunUMAP(merge, dims = 1:30)

merge <- SetupForWGCNA(
  merge,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "vis"
)

merge <- MetacellsByGroups(
  seurat_obj = merge,
  group.by = c("annotation", "sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'annotation' # set the Idents of the metacell seurat object
)

merge <- NormalizeMetacells(merge)

merge <- SetDatExpr(
  merge,
  group_name = c("L1", "L2/3", "L4", "L5", "L6", "WM"), # the name of the group of interest in the group.by column
  group.by='annotation', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

saveRDS(merge, "../Result/Stereo/WGCNA/WGCNA_merge.RDS")
merge_WGCNA = readRDS("../Result/Stereo/WGCNA/WGCNA_merge.RDS")

# Test different soft powers
merge<- TestSoftPowers(
  merge,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(merge)

# assemble with patchwork
#wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(merge_WGCNA)

merge_WGCNA <- ConstructNetwork(
  merge, soft_power=7,
  setDatExpr=FALSE,
  tom_name = 'all', # name of the topoligical overlap matrix written to disk,
  overwrite_tom = TRUE
)

# need to run ScaleData first or else harmony throws an error:
merge_WGCNA <- ScaleData(merge_WGCNA, features=VariableFeatures(merge_WGCNA))

# compute all MEs in the full single-cell dataset
merge_WGCNA <- ModuleEigengenes(
 merge_WGCNA,
 group.by.vars="diagnosis"
)

# harmonized module eigengenes:
hMEs <- GetMEs(merge_WGCNA)

# module eigengenes:
MEs <- GetMEs(merge_WGCNA, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
merge_WGCNA <- ModuleConnectivity(
  merge_WGCNA,
  group.by = 'annotation', group_name = c('L1', "L2/3", "L4", "L5", "L6", "WM")
)

# rename the modules
merge_WGCNA <- ResetModuleNames(
  merge_WGCNA,
  new_name = "Layer_modules"
)

# get the module assignment table:
modules <- GetModules(merge_WGCNA)

# show the first 6 columns:


hub_df <- GetHubGenes(merge_WGCNA, n_hubs = 10)

merge_WGCNA <- ModuleExprScore(
  merge_WGCNA,
  n_genes = 25,
  method='Seurat'
)

saveRDS(merge_WGCNA, "../Result/Stereo/WGCNA/WGCNA_merge.RDS")
