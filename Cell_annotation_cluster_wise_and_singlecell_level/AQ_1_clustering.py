# Read the Stereo-seq data
import stereo as st
import scanpy as sc
import pandas as pd

data = st.io.read_gef(file_path="/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.cellbin.gef",
                      bin_type='cell_bins')


## Filtering
data.tl.filter_cells(
        min_gene=50,
        # maximum number of counts required for a cell to pass fitlering. 
        max_gene=700,
        min_n_genes_by_counts=100, 
        # maximum number of genes expressed required for a cell to pass filtering.
        max_n_genes_by_counts=1100, 
        # Remove cells that have too many mitochondrial genes expressed, without enough genes expressed, and out of count range.
        # maximum number of pct_counts_mt required for a cell to pass filtering.
        pct_counts_mt=20,
        # whether to inplace the previous data or return a new data.
        inplace=True
        )

#data.tl.raw_checkpoint() # In order to save the data and recall it conveniently, you can save the raw expression matrix.

## Normalization (a combination method of normalize_total and log1p to normalize gene expression matrix)
data.tl.normalize_total()
data.tl.log1p()

# Highly variable genes
# Identify highly variable genes in cells.(In the subsequent "data.tl.pca" method, the parameter use_highly_genes can be set as True/False.)
data.tl.highly_variable_genes(
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=2000,
            res_key='highly_variable_genes'
            )

# save raw data
data.tl.raw_checkpoint()
data.tl.raw

# for the clustering
data.plt.highly_variable_genes(res_key='highly_variable_genes')
data.tl.scale()

data.tl.pca(
        use_highly_genes=True,
        n_pcs=30,
        res_key='pca'
)

data.tl.neighbors(pca_res_key='pca', res_key='neighbors')
data.tl.umap(
        pca_res_key='pca',
        neighbors_res_key='neighbors',
        res_key='umap'
        )
# 22 clusters
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution = 0.2)
data.plt.cluster_scatter(res_key='leiden')
data.plt.umap(res_key='umap', cluster_key='leiden')

# save the raw data
raw_data = data.tl.raw

# add the clustering information on the raw data
raw_data.tl.result["leiden"] = data.tl.result["leiden"]
raw_data.tl.raw_checkpoint()

st.io.stereo_to_anndata(raw_data,
                        flavor='seurat',
                        output='/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.anndata_with_clustering.h5ad')

##########################################################################################################################################################
# After getting annotation dictionary
import pandas as pd

# read text file into pandas DataFrame
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2_annotation_dict_cluster.txt", sep=" ")

# convert annotation dictionary into list
annotation_dict_cluster = anno_dict_cluster.dict.tolist()

data.tl.annotation(
        annotation_information=annotation_dict_cluster, ## annotation dictionary
        cluster_res_key='leiden',
        res_key='anno_cluster_leiden' ## store annotation in "res_key" as a keyword
        )

data.plt.umap(res_key='umap', cluster_key='anno_cluster_leiden')
