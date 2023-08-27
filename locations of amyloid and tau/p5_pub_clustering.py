import stereo as st

data = st.io.read_ann_h5ad(
       file_path='/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/B01809C2_extract.anndata.h5ad',
       spatial_key=None,
       )

data.tl.raw_checkpoint() # In order to save the data and recall it conveniently, you can save the raw expression matrix.

# Highly variable genes
# Identify highly variable genes in cells.(In the subsequent "data.tl.pca" method, the parameter use_highly_genes can be set as True/False.)
data.tl.highly_variable_genes(
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=5000,
            res_key='highly_variable_genes'
            )

# data.tl.scale() # Scale each gene to unit variance. Clip values exceeding standard deviation 10. 


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


# Clustering (Leiden) resolution=XXX ("Neuron" paper annotated 17 clusters)
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden')

## save StereoExpObject as AnnData in h5ad file
st.io.stereo_to_anndata(data,
                        flavor='seurat',
                        output='/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/B01809C2_extract_clusters.anndata.h5ad')
