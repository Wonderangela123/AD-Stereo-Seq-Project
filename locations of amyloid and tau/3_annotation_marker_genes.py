import stereo as st
import pandas as pd

# # read the h5ad file
# sample = "B01809C2"
# data = st.io.read_stereo_h5ad(
#         file_path='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}.stereo.h5ad'.format(sample, sample),
#         use_raw=True,
#         use_result=True
#         )

# recover the data
data = st.io.read_ann_h5ad(
       file_path='/work/aliu10/AD_Stereoseq_Project/processed/cases.h5ad',
       spatial_key=None,
       bin_type="cell_bins"
       )

data.tl.filter_cells(
        min_gene=20, 
        min_n_genes_by_counts=200,
        max_n_genes_by_counts=4000, 
        pct_counts_mt=10,
        inplace=True
        )

data.tl.raw_checkpoint() 

data.tl.normalize_total()
data.tl.log1p()

data.tl.highly_variable_genes(
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=2000,
            res_key='highly_variable_genes'
            )

data.tl.scale() # Scale each gene to unit variance. Clip values exceeding standard deviation 10. 

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

data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=0.1)
  
# read text file into pandas DataFrame
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed/annotation_dict_cluster.txt", sep=" ")

# convert annotation dictionary into list
annotation_dict_cluster = anno_dict_cluster.dict.tolist()

data.tl.annotation(
        annotation_information=annotation_dict_cluster, ## annotation dictionary
        cluster_res_key='leiden',
        res_key='anno_cluster_leiden' ## store annotation in "res_key" as a keyward
        )

# ## check if the same cell types are close, e.g., "Mic.7" is close to "Mic.16"
# data.plt.umap(res_key='umap', cluster_key='anno_cluster_leiden')

# Subset the cells annotated as "neuron"
id = data.tl.result['anno_cluster_leiden']['group'].str.contains('Ex|In')
id = id[id].index.tolist()

# # read annData object
# sample = "B01809C2"
# data1 = st.io.read_ann_h5ad(
#         file_path='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}.anndata.h5ad'.format(sample, sample),
#         spatial_key=None,
#         )

# data1.bin_type = "cell_bins"
data1 = data
data1.cells.cell_name = data1.cells.cell_name[id]
data1.exp_matrix = data1.exp_matrix[id]

# Reclustering
data1.tl.raw_checkpoint()

data1.tl.highly_variable_genes(
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=5000,
            res_key='highly_variable_genes'
            )

# data1.tl.scale() # will produce NaN

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

data1.tl.leiden(neighbors_res_key='neighbors',res_key='leiden')

# ## check if the cell types are seperate"
# data1.plt.umap(res_key='umap', cluster_key='leiden')

# # Find Marker Genes
# data1.tl.find_marker_genes(
#          cluster_res_key='leiden',
#          method='t_test',
#          use_highly_genes=True,
#          use_raw=True,
#          output="/work/aliu10/AD_Stereoseq_Project/processed_data/{}/Gene_markers_recluster.csv".format(sample)
#          )

## save StereoExpObject as AnnData in h5ad file
st.io.stereo_to_anndata(data1,
                        flavor='seurat',
                        output='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}_subtype.anndata.h5ad'.format(sample, sample))


