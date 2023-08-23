import stereo as st
import pandas as pd

## recover the StereoExpObject
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

# Subset the cells annotated as "neuron"
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

data = data1
##################################################################

# read text file into pandas DataFrame
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/annotation_dict_cluster_recluster.txt", sep=" ")

# convert annotation dictionary into list
annotation_dict_cluster = anno_dict_cluster.dict.tolist()

data.tl.annotation(
        annotation_information=annotation_dict_cluster, ## annotation dictionary
        cluster_res_key='leiden',
        res_key='anno_cluster_leiden' ## store annotation in "res_key" as a keyward
        )

# Subset the cells annotated as "neuron"
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


## save StereoExpObject as AnnData in h5ad file
data1.tl.raw_checkpoint() # convert to AnnData should have raw data
st.io.stereo_to_anndata(data1,
                        flavor='seurat',
                        output='/work/aliu10/AD_Stereoseq_Project/processed_data/{}/{}_extract.anndata.h5ad'.format(sample, sample))
