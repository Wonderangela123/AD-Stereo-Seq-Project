import stereo as st
import pandas as pd

data = st.io.read_gef(file_path="/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.cellbin.gef",
                      bin_type='cell_bins')

# data.tl.cal_qc()
# data.plt.violin()

## Filtering
data.tl.filter_cells(
        min_gene=100,
        # maximum number of counts required for a cell to pass fitlering. 
        # max_gene=800,
        min_n_genes_by_counts=10, 
        # maximum number of genes expressed required for a cell to pass filtering.
        # max_n_genes_by_counts=600, 
        # Remove cells that have too many mitochondrial genes expressed, without enough genes expressed, and out of count range.
        # maximum number of pct_counts_mt required for a cell to pass filtering.
        pct_counts_mt=5,
        # whether to inplace the previous data or return a new data.
        inplace=True
        )

data.tl.raw_checkpoint() # In order to save the data and recall it conveniently, you can save the raw expression matrix.

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

data.tl.scale() # Scale each gene to unit variance. Clip values exceeding standard deviation 10. 


# Embedding
## PCA (Principal component analysis)
data.tl.pca(
        use_highly_genes=True,
        n_pcs=30,
        res_key='pca'
        )

## Neighborhood graph: After PCA, we compute the neighborhood graph of cells using the PCA representation of the expression matrix.
data.tl.neighbors(pca_res_key='pca', res_key='neighbors')

## UMAP: It’s strongly to suggest embedding the graph in two dimensions using UMAP.
data.tl.umap(
        pca_res_key='pca',
        neighbors_res_key='neighbors',
        res_key='umap'
        )


# Clustering (Leiden)
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=0.5)

# read text file into pandas DataFrame
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2_annotation_dict_cluster.txt", sep=" ")

# convert annotation dictionary into list
annotation_dict_cluster = anno_dict_cluster.dict.tolist()

data.tl.annotation(
        annotation_information=annotation_dict_cluster, ## annotation dictionary
        cluster_res_key='leiden',
        res_key='anno_cluster_leiden' ## store annotation in "res_key" as a keyword
        )

st.io.write_h5ad(
        data,
        use_raw=True,
        use_result=True,
        output='/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.stereo.h5ad',
        )

# Subset the cells annotated as "neuron"
id = data.tl.result['anno_cluster_leiden']['group'].str.contains('Ex|In')
id = id[id].index.tolist()

data1 = data
data1.cells.cell_name = data.cells.cell_name[id]
data1.exp_matrix = data.exp_matrix[id]
data1.position = data.position[id, :]

data1.tl.raw_checkpoint() 

data1.tl.quantile()

data1.tl.highly_variable_genes(
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=2000,
            res_key='highly_variable_genes'
            )

data1.tl.scale() # Scale each gene to unit variance. Clip values exceeding standard deviation 10. 

data1.tl.pca(
        use_highly_genes=True,
        n_pcs=30,
        res_key='pca'
        )

data1.tl.neighbors(pca_res_key='pca', res_key='neighbors')

data1.tl.umap(
        pca_res_key='pca',
        neighbors_res_key='neighbors',
        res_key='umap'
        )

data1.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=0.5)

# save StereoExpObject as AnnData in h5ad file
st.io.stereo_to_anndata(data1,
                        flavor='seurat',
                        output='/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2_subtype.anndata.h5ad')
