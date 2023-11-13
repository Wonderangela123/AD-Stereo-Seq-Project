# We will use the scanpy to analysis the data
import scanpy as sc
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import gc

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
    sc.external.pp.harmony_integrate(adata, "sample_ID")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep= "X_pca_harmony")
    sc.tl.umap(adata)
    return adata

# read the anndata
data_list = []
for file in os.listdir("../Bin50_result/anndata/raw/"):
    data = sc.read_h5ad("../Bin50_result/anndata/raw/" + file)
    data.obs["sample_ID"] = file.split(".")[0]
    data.obs["tech"] = "stereo"
    globals()[file.split(".")[0]] = data
    print(file.split(".")[0])

case = ad.concat([A02092E1, B02009F6, C02248B5, B02008C6, B02008D2, B01809C2], join = "outer")
case.obs_names_make_unique()
control = ad.concat([D02175A4, D02175A6, B01809A3, B01809A4, B01806B5, B01806B6], join = "outer")
control.obs_names_make_unique()

del A02092E1, B02009F6, C02248B5, B02008C6, B02008D2, B01809C2,D02175A4, D02175A6, B01809A3, B01809A4, B01806B5, B01806B6
gc.collect()

case.obs["diagnosis"] = "case"
control.obs["diagnosis"] = "control"

# merge the case and control data
merge = ad.concat([case, control], join = "outer")

# save the merge data
if not os.path.exists("../Bin50_result/anndata/one_raw/"):
    os.mkdir("../Bin50_result/anndata/one_raw/")
merge.write_h5ad("../Bin50_result/anndata/one_raw/raw_merged.h5ad")

del case, control
gc.collect()

merge = data_process(merge)
















