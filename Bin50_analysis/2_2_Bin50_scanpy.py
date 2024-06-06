import stereo as st
import os
import pandas as pd
import gc
import utils
from importlib import reload
reload(utils)
import pandas as pd
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
    sc.tl.leiden(adata, resolution = 0.1)
    return adata

merge = st.io.read_stereo_h5ad(
    file_path="../Result/Stereo/One_file/bin50_integrate.h5ad",
    use_raw=True,
    use_result=True
)
meta = pd.read_csv("../Result/Stereo/One_file/meta_info")
utils.cell_t0_res_key(merge)

ref = sc.read("../../MIT_snRNAseqPFC_BA10/MIT_BA10_ref.h5ad")
ref.obs["sample"] = ref.obs["projid"].astype("str")

merge_ann = st.io.stereo_to_anndata(merge,flavor='scanpy', split_batches=False)
for i in meta.columns:
    merge_ann.obs[i] = meta[i].to_list()
    

merge_ref = data_process(merge_ref)

merge_ref.write_h5ad("../Result/Anndata/Integrate_ref/merge_ref_test.h5ad")