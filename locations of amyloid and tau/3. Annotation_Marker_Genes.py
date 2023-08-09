import stereo as st
import pandas as pd


# read the Anndata h5ad file
sample = "B01809C2"
ann_h5ad = '/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}.stereo.h5ad'.format(sample, sample)
data = st.io.read_ann_h5ad(
        file_path=ann_h5ad,
        spatial_key=None,
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
