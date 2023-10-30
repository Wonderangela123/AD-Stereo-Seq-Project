import stereo as st
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

# extract the expression data from the stereo folders
data_list = []
for num,files in enumerate(os.listdir("../processed_data/")):
    print(files)
    data_path = "../processed_data/{}/GeneExpMatrix/{}.cellbin.gef".format(files, files)
    print(data_path)
    data_name = files.split(".")[0]
    data = st.io.read_gef(file_path=data_path, bin_type="cell_bins")
    globals()[data_name] = data 

# merge the data
data_combined = st.utils.data_helper.merge(A02092E1, B02009F6, C02248B5, B02008C6, B02008D2, B01809C2,
                                         D02175A4, D02175A6, B01809A3, B01809A4, B01806B5, B01806B6)

# calculate the total gene counts, detected genes, and percentage of mitochondria
data_combined.tl.cal_qc()

# test the gene expression
total_counts = pd.DataFrame(data_combined.tl.result["total_counts"])
genes_counts = pd.DataFrame(data_combined.tl.result["n_genes_by_counts"])
pct_mt = pd.DataFrame(data_combined.tl.result["pct_counts_mt"])

# Check 
# total counts
sns.violinplot(y=total_counts["total_counts"])
plt.axhline(y=100, color='r', linestyle='--')
plt.axhline(y=1000, color='g', linestyle='-.')
plt.title('Total Counts')

# Genes per cell
sns.violinplot(y=genes_counts["n_genes_by_counts"])
plt.axhline(y=25, color='r', linestyle='--')
plt.axhline(y=800, color='g', linestyle='-.')
plt.title('Genes Per Cells')

# Mitochondria per cell
sns.violinplot(y=pct_mt["pct_counts_mt"])
plt.axhline(y=25, color='g', linestyle='-.')
plt.title('MT_PCT_Cell')

# filter the cells
data_combined.tl.filter_cells(
        min_gene=100,
        max_gene=1000,
        min_n_genes_by_counts=25,
        max_n_genes_by_counts=800,
        pct_counts_mt=25,
        inplace=True
        )

# save the check point
data_combined.tl.raw_checkpoint()

# data normalization
data_combined.tl.normalize_total()
data_combined.tl.log1p()

# Perform the PCA
data_combined.tl.pca(use_highly_genes=False, n_pcs=50, res_key='pca')

# Perform the harmony integration
data_combined.tl.batches_integrate(pca_res_key='pca', res_key='pca_integrated')

# Calculate the neighbors
data_combined.tl.neighbors(pca_res_key='pca_integrated', n_pcs=50, res_key='neighbors_integrated')

# Calculate the spatial neighbors (for spatial clustering)
data_combined.tl.spatial_neighbors(
    neighbors_res_key = "neighbors_integrated",
    res_key = "spatial_neighbors"
)

# Run and plot UMAP
data_combined.tl.umap(pca_res_key='pca_integrated', neighbors_res_key='neighbors_integrated', res_key='umap_integrated')
data_combined.plt.batches_umap(res_key='umap_integrated', )

# Calculate the spatial leiden
data_combined.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden')
data_combined.plt.cluster_scatter(res_key='spatial_leiden')

# save as stereo ann file (one combine file)
st.io.write_h5ad(
        data_combined,
        use_raw=True,
        use_result=True,
        key_record=None,
        output = "Result/cell_bin_combined.h5ad",
        split_batches = False
        )

# write the scanpy anndata file 
adata = st.io.stereo_to_anndata(data,flavor='scanpy',output='Result/integrated_data.h5ad')


















































