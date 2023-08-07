import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
import stereo as st
import numpy as np
from stereo.core.stereo_exp_data import AnnBasedStereoExpData
import pyarrow.feather as feather
import rpy2
import rpy2.robjects as robjects

# to use GraphST for the spatial cluster
def run_GraphST(sample, bin_size, celltype=False, data_type = "Visium"): # set any number you want if you choose celltype = True

    # set the gpu
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    n_clusters = 7

    # load the data
    if celltype == True:
        data_path = "../{}/GeneExpMatrix/{}.{}.gef".format(sample, sample, "cellbin")
        st.io.read_gef_info(data_path)
        data1 = st.io.read_gef(file_path=data_path, bin_type = "cell_bins")
    
    else:
        data_path = "../{}/GeneExpMatrix/{}.{}.gef".format(sample, sample, "tissue")
        data1 = st.io.read_gef(file_path=data_path, bin_size = bin_size)
    
    
    # conver data into anndata
    data1.tl.cal_qc()
    data1.tl.raw_checkpoint()
    adata = st.io.stereo_to_anndata(data1,flavor='seurat',output='seurat_out.h5ad')

    # define model
    if data_type == "Stereo":
        model = GraphST.GraphST(adata, datatype='Stereo', device=device)

    elif data_type == "Visium":
        model = GraphST.GraphST(adata, datatype,device=device)

    # run model
    adata = model.train()

    from GraphST.utils import clustering
    
    tool = 'mclust' # mclust, leiden, and louvain
    
    # clustering
    from GraphST.utils import clustering
    
    if tool == 'mclust':
        clustering(adata, n_clusters, method=tool)
    elif tool in ['leiden', 'louvain']:
        clustering(adata, n_clusters, method=tool, start=0.1, end=2.0, increment=0.01)
    
    adata.write_h5ad("../Result/Graph_{}_cluster_{}_{}.h5ad".format(sample, celltype, data_type))
    pass    

if __name__ == "__main__":
    # test the function
    # data are saved in ../Result folder
    run_GraphST("B01806B5", celltype=True, bin_size = 100, data_type="Stereo")
























































    





