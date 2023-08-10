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

# Subset the cells annotated as "microglia"
id = data.tl.result['anno_cluster_leiden']['group'].str.contains('Mic')
id = id[id].index.tolist()

data.cells.cell_name = data.cells.cell_name[id]
data.exp_matrix = data.exp_matrix[id]

# Reclustering
data.tl.raw_checkpoint()

data.tl.highly_variable_genes(
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=5000,
            res_key='highly_variable_genes'
            )

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

data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden')












# ## cluster-level annotations
# data.plt.cluster_scatter(res_key='anno_cluster_leiden')

# ## check if the same cell types are close, e.g., "Ex.1" is close to "Ex.3"
# data.plt.umap(res_key='umap', cluster_key='anno_cluster_leiden')



# ## gene expression regulation
# markers = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/Gene_markers.csv", sep=",", index_col=0)

# ## microglia
# mic = markers.loc[:,["7_pvalues_adj", "7_log2fc", "16_pvalues_adj", "16_log2fc"]]  ## extract microglia (cluster 7 and 16)
# mic = mic.loc[mic.loc[:,"7_pvalues_adj"] < 0.05] ## set the singnificance level as 0.05
# mic = mic.loc[mic.loc[:,"16_pvalues_adj"] < 0.05]

# gene_name = pd.DataFrame(data.genes.gene_name)

# microglia = pd.merge(gene_name, mic, left_index=True, right_index=True)
# microglia.to_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/gene_expression_microglia.csv", sep = " ", index = False)


