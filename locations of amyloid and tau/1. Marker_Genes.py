# Read the Stereo-seq data
import stereo as st

sample = "B01809C2" ## AD case
data_path = "/work/ygong/stereo_seq_public/{}/GeneExpMatrix/{}.cellbin.gef".format(sample, sample)
st.io.read_gef_info(data_path)
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins')


# Data preprocessing includes three modules: quality control, filtering and normalization.

## Quality control: This module calculates the quality distribution of original data, using three indicators: 
## 1) total_counts - the total counts per cell
## 2) n_genes_by_counts - the number of genes expressed in count maxtrix
## 3) pct_counts_mt - the percentage of counts in mitochondrial genes

## Filtering
data.tl.filter_cells(
        # minimum number of counts required for a cell to pass fitlering.
        min_gene=20, 
        # minimum number of genes expressed required for a cell to pass filtering.
        min_n_genes_by_counts=200,
        # maximum number of genes expressed required for a cell to pass filtering.
        max_n_genes_by_counts=4000, 
        # Remove cells that have too many mitochondrial genes expressed, without enough genes expressed, and out of count range.
        # maximum number of pct_counts_mt required for a cell to pass filtering.
        pct_counts_mt=10,
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
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden')


# Find Marker Genes
data.tl.find_marker_genes(
        cluster_res_key='leiden',
        method='t_test',
        use_highly_genes=False,
        use_raw=True,
        output="/work/aliu10/AD_Stereoseq_Project/processed_data/{}/Gene_markers.csv".format(sample)
        )
### filter out genes based on log fold change and fraction of genes expressing the gene within and outside each group (optional)
data.tl.filter_marker_genes(
    marker_genes_res_key='marker_genes',
    min_fold_change=1,
    min_in_group_fraction=0.25,
    max_out_group_fraction=0.5,
    res_key='marker_genes_filtered',
    output="/work/aliu10/AD_Stereoseq_Project/processed_data/{}/Gene_markers.csv".format(sample)
)

## save StereoExpObject as AnnData in h5ad file
st.io.stereo_to_anndata(data,
                        flavor='seurat',
                        output='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}.anndata.h5ad'.format(sample, sample))


## write a new h5ad with StereoExpData, if key_record = None, it will use the res_key stored in data.tl.key_record
st.io.write_h5ad(data,
                 use_raw=True,
                 use_result=True,
                 key_record=None,
                 output='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}.stereo.h5ad'.format(sample, sample))


