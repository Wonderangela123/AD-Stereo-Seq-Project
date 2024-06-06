import stereo as st
import os
import pandas as pd
import gc
import utils
from importlib import reload
reload(utils)

######### read the data as ann h5ad files
control_list = ["B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6"]
moderate_list = ["B02008D2", "B02009F6"]
advanced_list = ["B01809C2", "C02248B5"]
severe_list = ["A02092E1", "B02008C6"]

###### Control
# import the stereo data
for num,files in enumerate(control_list):
    data_path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(files, files)
    data = st.io.read_gef(file_path=data_path, bin_size=50)
    data.cells["diagnosis"] = "control"
    data.cells["levels"] = "control"
    data.cells["sample"] = files
    globals()[files] = data
    
##### moderate
# import the stereo data
for num,files in enumerate(moderate_list):
    data_path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(files, files)
    data = st.io.read_gef(file_path=data_path, bin_size=50)
    data.cells["diagnosis"] = "case"
    data.cells["levels"] = "moderate"
    data.cells["sample"] = files
    globals()[files] = data

#### advanced
for num,files in enumerate(advanced_list):
    data_path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(files, files)
    data = st.io.read_gef(file_path=data_path, bin_size=50)
    data.cells["diagnosis"] = "case"
    data.cells["levels"] = "advanced"
    data.cells["sample"] = files
    globals()[files] = data

#### severe_list
for num,files in enumerate(severe_list):
    data_path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(files, files)
    data = st.io.read_gef(file_path=data_path, bin_size=50)
    data.cells["diagnosis"] = "case"
    data.cells["levels"] = "severe"
    data.cells["sample"] = files
    globals()[files] = data

data = st.utils.data_helper.merge(B01806B5, B01806B6, B01809A3, B01809A4, D02175A4, D02175A6,
                                 A02092E1, B02009F6, C02248B5, B02008C6, B02008D2, B01809C2)
    
data.tl.cal_qc()
data.plt.violin()

data.tl.filter_cells(
        min_n_genes_by_counts=400,
        inplace=True
        )

data.tl.raw_checkpoint()
    
data.tl.normalize_total()
data.tl.log1p()    
    
data.tl.pca(use_highly_genes=False, n_pcs=50, res_key='pca')   
    
data.tl.batches_integrate(pca_res_key='pca', res_key='pca_integrated')

data.tl.neighbors(
        pca_res_key='pca_integrated',
        n_pcs=30,
        res_key='neighbors'
        )

data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')

data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden', resolution=1)
data.plt.cluster_scatter(res_key='leiden', reorganize_coordinate = 2)

# one file save
meta = data.cells.to_df()
meta.to_csv("../Result/Stereo/One_file/meta_info")
st.io.write_h5ad(
        data,
        use_raw=True,
        use_result=True,
        key_record=None,
        output = "../Result/Stereo/One_file/bin50_integrate.h5ad",
        split_batches = False
        )