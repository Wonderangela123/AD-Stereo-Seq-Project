import stereo as st
import os
import pandas as pd
import gc
import utils
from importlib import reload
reload(utils)
import scanpy as sc
import anndata as ad

def data_process(adata):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10, zero_center=False)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.external.pp.harmony_integrate(adata, "sample")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep= "X_pca_harmony")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution = 0.05)
    return adata

######### read the data as ann h5ad files
control_list = ["B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6"]
moderate_list = ["B02008D2", "B02009F6"]
advanced_list = ["B01809C2", "C02248B5"]
severe_list = ["A02092E1", "B02008C6"]

###### Control
# import the stereo data
for num,files in enumerate(control_list):
    data_path = "../Result/Anndata/Raw/Multi_file/{}_raw.h5ad".format(files)
    data = sc.read_h5ad(data_path)
    data.obs["diagnosis"] = "control"
    data.obs["levels"] = "control"
    data.obs["sample"] = files
    globals()[files] = data
    
##### moderate
# import the stereo data
for num,files in enumerate(moderate_list):
    data_path = "../Result/Anndata/Raw/Multi_file/{}_raw.h5ad".format(files)
    data = sc.read_h5ad(data_path)
    data.obs["diagnosis"] = "case"
    data.obs["levels"] = "moderate"
    data.obs["sample"] = files
    globals()[files] = data
    
#### advanced
for num,files in enumerate(advanced_list):
    data_path = "../Result/Anndata/Raw/Multi_file/{}_raw.h5ad".format(files)
    data = sc.read_h5ad(data_path)
    data.obs["diagnosis"] = "case"
    data.obs["levels"] = "advanced"
    data.obs["sample"] = files
    globals()[files] = data
    

#### severe_list
for num,files in enumerate(severe_list):
    data_path = "../Result/Anndata/Raw/Multi_file/{}_raw.h5ad".format(files)
    data = sc.read_h5ad(data_path)
    data.obs["diagnosis"] = "case"
    data.obs["levels"] = "severe"
    data.obs["sample"] = files
    globals()[files] = data

merge = ad.concat([B01806B5, B01806B6, B01809A3, B01809A4, D02175A4, D02175A6,
                  B02008D2, B02009F6, B01809C2, C02248B5, A02092E1, B02008C6],
                    join = "outer"
                 )
merge.var['mt'] = merge.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(merge, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

merge = merge[merge.obs["n_genes_by_counts"]>400]
merge = data_process(merge)

merge.write_h5ad("../Result/Anndata/New_integrate/merge.h5ad")