import stereo as st
import pandas as pd

# read the h5ad file
sample = "B01809C2"
data = st.io.read_stereo_h5ad(
        file_path='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}_recluster.stereo.h5ad'.format(sample, sample),
        use_raw=True,
        use_result=True
        )
  
# read text file into pandas DataFrame
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/annotation_dict_cluster_recluster.txt", sep=" ")

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
        file_path='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}_recluster.anndata.h5ad'.format(sample, sample),
        spatial_key=None,
        )

data1.bin_type = "cell_bins"
data1.cells.cell_name = data1.cells.cell_name[id]
data1.exp_matrix = data1.exp_matrix[id]

# Find Marker Genes
data1.tl.find_marker_genes(
         cluster_res_key='leiden',
         method='t_test',
         use_highly_genes=True,
         use_raw=True,
         output="/work/aliu10/AD_Stereoseq_Project/processed_data/{}/Gene_markers_extract.csv".format(sample)
         )

## save StereoExpObject as AnnData in h5ad file
st.io.stereo_to_anndata(data1,
                        flavor='seurat',
                        output='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}_extract.anndata.h5ad'.format(sample, sample))


## write a new h5ad with StereoExpData, if key_record = None, it will use the res_key stored in data.tl.key_record
st.io.write_h5ad(data1,
                 use_raw=True,
                 use_result=True,
                 key_record=None,
                 output='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}_extract.stereo.h5ad'.format(sample, sample))

