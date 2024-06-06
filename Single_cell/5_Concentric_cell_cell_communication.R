library(CellChat)
library(Seurat)
library(Matrix)

cell_chat_prediction = function(cellchat)
{
    cellchat <- setIdent(cellchat, ident.use = "annotation") # set "labels" as default cell identity
    levels(cellchat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

    # set the ligand-receptor interaction database
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)

    dplyr::glimpse(CellChatDB$interaction)

    # use a subset of CellChatDB for cell-cell communication analysis
    #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB

    # set the used database in the object
    cellchat@DB <- CellChatDB.use

    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 8) # do parallel
    options(future.globals.maxSize = 2000 * 10000000^2)  # Setting limit to 2 GiB

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 2)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    return(cellchat)
}

# set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# read the matrix
matrix = readMM("../Result/Concentric_to_R/concentric_matrix.mtx")

# read the cell and gene names
meta = read.csv("../Result/Concentric_to_R/meta_data.csv")
feature_name = read.csv("../Result/Concentric_to_R/raw_feature.csv")

rownames(matrix) = meta$X
colnames(matrix) = rownames(feature_name)

meta_info = read.csv("../Result/Concentric_to_R/meta_data.csv")

data = CreateSeuratObject(counts = t(matrix), project = "pbmc3k", min.cells = 0, min.features = 0)

# add the meta information
for (i in meta_info){
    for (i in colnames(meta_info)){
        data[[i]] = meta_info[[i]]
    }
}

data_inner = subset(data, concentric == "inner")
data_outer = subset(data, concentric == "outer2")

# separate them into different groups
for (i in c("control", "moderate", "advanced", "severe"))
{
    tmp = subset(data_inner, subset = levels == i)
    count = tmp@assays$RNA@counts
    meta = tmp@meta.data
    cellchat <- createCellChat(object = count, meta = meta, group.by = "annotation")
    cellchat = cell_chat_prediction(cellchat)
    
    assign(paste0(i, "_inner_cell_chat"), cellchat)
    print(i)
}

# separate them into different groups
for (i in c("control", "moderate", "advanced", "severe"))
{
    tmp = subset(data_outer, subset = levels == i)
    count = tmp@assays$RNA@counts
    meta = tmp@meta.data
    cellchat <- createCellChat(object = count, meta = meta, group.by = "annotation")
    cellchat = cell_chat_prediction(cellchat)
    
    assign(paste0(i, "_outer_cell_chat"), cellchat)
    print(i)
}

object.list <- list(
    control_outer = control_outer_cell_chat, 
    moderate_outer = moderate_outer_cell_chat,
    advanced_outer = advanced_outer_cell_chat,
    severe_outer = severe_outer_cell_chat,
    control_inner = control_inner_cell_chat, 
    moderate_inner = moderate_inner_cell_chat,
    advanced_inner = advanced_inner_cell_chat,
    severe_inner = severe_inner_cell_chat
)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(control_outer_cell_chat, "../Result/cell_cell_communication/control_outer.RDS")
saveRDS(moderate_outer_cell_chat, "../Result/cell_cell_communication/moderate_outer.RDS")
saveRDS(advanced_outer_cell_chat, "../Result/cell_cell_communication/advanced_outer.RDS")
saveRDS(severe_outer_cell_chat, "../Result/cell_cell_communication/severe_outer.RDS")
saveRDS(control_inner_cell_chat, "../Result/cell_cell_communication/control_inner.RDS")
saveRDS(moderate_inner_cell_chat, "../Result/cell_cell_communication/moderate_inner.RDS")
saveRDS(advanced_inner_cell_chat, "../Result/cell_cell_communication/advanced_inner.RDS")
saveRDS(severe_inner_cell_chat, "../Result/cell_cell_communication/severe_inner.RDS")
saveRDS(cellchat, "../Result/cell_cell_communication/cellchat.RDS")


