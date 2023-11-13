import stereo as st
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import gc

# create a folder
if not os.path.exists("../Bin50_result"):
    os.mkdir("../Bin50_result")

if not os.path.exists("../Bin50_result/anndata"):
    os.mkdir("../Bin50_result/anndata")

if not os.path.exists("../Bin50_result/anndata/raw"):
    os.mkdir("../Bin50_result/anndata/raw")
    
for sample in os.listdir("../../processed_data"):
    print(sample)
    path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(sample, sample)
    
    # import the data as bin50
    tmp = st.io.read_gef(file_path=path, bin_size=50)

    # perform the QC
    tmp.tl.cal_qc()
    tmp.tl.filter_cells(
        min_gene=200,
        max_gene=3000,
        min_n_genes_by_counts=25,
        max_n_genes_by_counts=2000,
        pct_counts_mt=15,
        inplace=True
        )

    # check the raw data
    tmp.tl.raw_checkpoint()

    st.io.stereo_to_anndata(
        tmp,
        flavor='scanpy',
        output='../Bin50_result/anndata/raw/{}.h5ad'.format(sample)
    )
    globals()[sample] = tmp

# Data integration
data = st.utils.data_helper.merge(A02092E1,B01806B5,B01806B6,B01809A3,
                                  B01809A4,B01809C2,B02008C6,B02008D2,
                                  B02009F6,C02248B5,D02175A4,D02175A6)
# remove the single files and release the memory
del A02092E1,B01806B5,B01806B6,B01809A3,B01809A4,B01809C2,B02008C6,B02008D2,B02009F6,C02248B5,D02175A4,D02175A6
gc.collect()

# Calculate the quality control

# test the gene expression
#total_counts = pd.DataFrame(data.tl.result["total_counts"])
#genes_counts = pd.DataFrame(data.tl.result["n_genes_by_counts"])
#pct_mt = pd.DataFrame(data.tl.result["pct_counts_mt"])

# Check 
# total counts
#sns.violinplot(y=total_counts["total_counts"])
#plt.axhline(y=200, color='r', linestyle='--')
#plt.axhline(y=3000, color='g', linestyle='-.')
#plt.title('Total Counts')

# Genes per cell
#sns.violinplot(y=genes_counts["n_genes_by_counts"])
#plt.axhline(y=25, color='r', linestyle='--')
#plt.axhline(y=2000, color='g', linestyle='-.')
#plt.title('Genes Per Cells')

# Mitochondria per cell
#sns.violinplot(y=pct_mt["pct_counts_mt"])
#plt.axhline(y=15, color='g', linestyle='-.')
#plt.title('MT_PCT_Cell')

# filter the cells
"""
data.tl.filter_cells(
        min_gene=200,
        max_gene=3000,
        min_n_genes_by_counts=25,
        max_n_genes_by_counts=2000,
        pct_counts_mt=15,
        inplace=True
        )
"""
data.tl.raw_checkpoint()

data.tl.normalize_total()
data.tl.log1p()

# Perform the PCA
data.tl.pca(use_highly_genes=False, n_pcs=50, res_key='pca')

# Perform the harmony integration
data.tl.batches_integrate(pca_res_key='pca', res_key='pca_integrated')

data.tl.neighbors(pca_res_key='pca_integrated', n_pcs=50, res_key='neighbors_integrated')
data.tl.umap(pca_res_key='pca_integrated', neighbors_res_key='neighbors_integrated', res_key='umap_integrated')
data.plt.batches_umap(res_key='umap_integrated')

if not os.path.exists("../Bin50_result/Stereo_seq"):
    os.mkdir("../Bin50_result/Stereo_seq")    
if not os.path.exists("../Bin50_result/Stereo_seq/One_file"):
    os.mkdir("../Bin50_result/Stereo_seq/One_file")
if not os.path.exists("../Bin50_result/Stereo_seq/One_file/"):
    os.mkdir("../Bin50_result/Stereo_seq/One_file/")
if not os.path.exists("../Bin50_result/Stereo_seq/Multiple_file/"):
    os.mkdir("../Bin50_result/Stereo_seq/Multiple_file/")


st.io.write_h5ad(
        data,
        use_raw=True,
        use_result=True,
        key_record=None,
        output = "../Bin50_result/Stereo_seq/One_file/bin50_integrated.h5ad",
        split_batches = False
)

st.io.write_h5ad(
        data,
        use_raw=True,
        use_result=True,
        key_record=None,
        output = "../Bin50_result/Stereo_seq/Multiple_file/bin50_integrated.h5ad",
        split_batches = True
)











































