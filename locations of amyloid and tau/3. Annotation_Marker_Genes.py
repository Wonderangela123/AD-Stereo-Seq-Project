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

## cluster-level annotations
data.plt.cluster_scatter(res_key='anno_cluster_leiden')

## check if the same cell types are close, e.g., "Ex.1" is close to "Ex.3"
data.plt.umap(res_key='umap', cluster_key='anno_cluster_leiden')



## gene expression regulation
markers = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/Gene_markers.csv", sep=",", index_col=0)

## microglia
mic = markers.loc[:,["7_pvalues_adj", "7_log2fc", "16_pvalues_adj", "16_log2fc"]]  ## extract microglia (cluster 7 and 16)
mic = mic.loc[mic.loc[:,"7_log2fc"] * mic.loc[:,"16_log2fc"] > 0] ## keep only the genes with the same expression regulation patterns in various clusters
mic = mic.loc[mic.loc[:,"7_pvalues_adj"] < 0.05] ## set the singnificance level as 0.05
mic = mic.loc[mic.loc[:,"16_pvalues_adj"] < 0.05]

import numpy as np
mic['regulation'] = np.where(mic['7_log2fc'] > 0, "Upregulation", "Downregulation")
mic.iloc[:,-1].to_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/gene_regulation_microglia.csv", sep = " ")

