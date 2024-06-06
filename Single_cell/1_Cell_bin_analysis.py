import stereo as st
import os
import gc
import numpy as np
import pandas as pd
from stereo.utils.data_helper import split
import numpy as np
import utils
from importlib import reload
import data_processing
import scanpy as sc
import anndata as ad
reload(utils)
reload(data_processing)

# load the data and convert into h5ad
case = ["A02092E1", "B02009F6", "C02248B5", "B02008C6", "B02008D2", "B01809C2"]
control = ["D02175A4", "D02175A6", "B01809A3", "B01809A4", "B01806B5", "B01806B6"]
for sample in case:
    data_path = "../Data/{}/Result/{}.cellbin.gef".format(sample, sample)
    tmp = st.io.read_gef(file_path=data_path, bin_type='cell_bins')
    tmp.tl.cal_qc()
    tmp.tl.raw_checkpoint()
    tmp.cells["sample_id"] = sample
    tmp.cells["diagnosis"] = "case"
    globals()[sample] = tmp
    st.io.stereo_to_anndata(tmp,flavor='scanpy', output = "../Result/anndata/raw_data/{}.h5ad".format(sample))

for sample in control:
    data_path = "../Data/{}/Result/{}.cellbin.gef".format(sample, sample)
    tmp = st.io.read_gef(file_path=data_path, bin_type='cell_bins')
    tmp.tl.cal_qc()
    tmp.tl.raw_checkpoint()
    tmp.cells["sample_id"] = sample
    tmp.cells["diagnosis"] = "control"
    globals()[sample] = tmp
    st.io.stereo_to_anndata(tmp,flavor='scanpy', output = "../Result/anndata/raw_data/{}.h5ad".format(sample))

######### read the data as ann h5ad files
control_list = ["B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6"]
moderate_list = ["B02008D2", "B02009F6"]
advanced_list = ["B01809C2", "C02248B5"]
severe_list = ["A02092E1", "B02008C6"]
###### Control
# import the stereo data
for num,files in enumerate(control_list):
    data_path = "../Result/anndata/raw_data/{}.h5ad".format(files)
    data = sc.read_h5ad(data_path)
    data.obs["diagnosis"] = "control"
    data.obs["sample"] = files
    data.obs["levels"] = "control"
    globals()[files] = data
    print(files)

###### moderate_list
# import the stereo data
for num,files in enumerate(moderate_list):
    data_path = "../Result/anndata/raw_data/{}.h5ad".format(files)
    data = sc.read_h5ad(data_path)
    data.obs["diagnosis"] = "case"
    data.obs["sample"] = files
    data.obs["levels"] = "moderate"
    globals()[files] = data
    print(files)

##### advanced_list
# import the stereo data
for num,files in enumerate(advanced_list):
    data_path = "../Result/anndata/raw_data/{}.h5ad".format(files)
    data = sc.read_h5ad(data_path)
    data.obs["diagnosis"] = "case"
    data.obs["sample"] = files
    data.obs["levels"] = "advanced"
    globals()[files] = data
    print(files)

##### severe_list
# import the stereo data
for num,files in enumerate(severe_list):
    data_path = "../Result/anndata/raw_data/{}.h5ad".format(files)
    data = sc.read_h5ad(data_path)
    data.obs["diagnosis"] = "case"
    data.obs["sample"] = files
    data.obs["levels"] = "severe"
    globals()[files] = data
    print(files)

merge = ad.concat(
    [B01806B5, B01806B6, B01809A3, B01809A4, 
     D02175A4, D02175A6,A02092E1, B02009F6, 
     C02248B5, B02008C6, B02008D2, B01809C2],
    join = "outer"
)

merge.obs_names_make_unique()
merge = merge[merge.obs["n_genes_by_counts"]>=150]
merge.raw = merge

merge = data_processing.data_process(merge, "sample")

merge.write_h5ad("../Result/anndata/integrated_data/integrated.h5ad")





























