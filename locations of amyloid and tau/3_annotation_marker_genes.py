import stereo as st
import pandas as pd

# read the h5ad file
sample = "B01809C2"
data = st.io.read_stereo_h5ad(
        file_path='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}.stereo.h5ad'.format(sample, sample),
        use_raw=True,
        use_result=True
        )
  
# read text file into pandas DataFrame
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/annotation_dict_cluster.txt", sep=" ")

# convert annotation dictionary into list
annotation_dict_cluster = anno_dict_cluster.dict.tolist()

data.tl.annotation(
        annotation_information=annotation_dict_cluster, ## annotation dictionary
        cluster_res_key='leiden',
        res_key='anno_cluster_leiden' ## store annotation in "res_key" as a keyward
        )

# ## check if the same cell types are close, e.g., "Mic.7" is close to "Mic.16"
# data.plt.umap(res_key='umap', cluster_key='anno_cluster_leiden')

# Subset the cells annotated as "microglia"
id = data.tl.result['anno_cluster_leiden']['group'].str.contains('Ex|In')
id = id[id].index.tolist()

# read annData object
sample = "B01809C2"
data1 = st.io.read_ann_h5ad(
        file_path='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}.anndata.h5ad'.format(sample, sample),
        spatial_key=None,
        )

data1.bin_type = "cell_bins"
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

data1.tl.scale() 

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

# Find Marker Genes
data1.tl.find_marker_genes(
         cluster_res_key='leiden',
         method='t_test',
         use_highly_genes=True,
         use_raw=True,
         output="/work/aliu10/AD_Stereoseq_Project/processed_data/{}/Gene_markers_recluster.csv".format(sample)
         )

## save StereoExpObject as AnnData in h5ad file
st.io.stereo_to_anndata(data1,
                        flavor='seurat',
                        output='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}_recluster.anndata.h5ad'.format(sample, sample))


## write a new h5ad with StereoExpData, if key_record = None, it will use the res_key stored in data.tl.key_record
st.io.write_h5ad(data1,
                 use_raw=True,
                 use_result=True,
                 key_record=None,
                 output='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}_recluster.stereo.h5ad'.format(sample, sample))





# ## gene expression regulation
# markers = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/Gene_markers.csv", sep=",", index_col=0)

# ## microglia
# mic = markers.loc[:,["7_pvalues_adj", "7_log2fc", "16_pvalues_adj", "16_log2fc"]]  ## extract microglia (cluster 7 and 16)
# mic = mic.loc[mic.loc[:,"7_pvalues_adj"] < 0.05] ## set the singnificance level as 0.05
# mic = mic.loc[mic.loc[:,"16_pvalues_adj"] < 0.05]

# gene_name = pd.DataFrame(data.genes.gene_name)

# microglia = pd.merge(gene_name, mic, left_index=True, right_index=True)
# microglia.to_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/gene_expression_microglia.csv", sep = " ", index = False)


