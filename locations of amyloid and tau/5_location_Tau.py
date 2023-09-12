import stereo as st
import pandas as pd

data = st.io.read_ann_h5ad(
       file_path='/work/aliu10/AD_Stereoseq_Project/processed/data/integrated.h5ad',
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

data.tl.quantile()

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
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed/data/annotation_dict_cluster.txt", sep=" ")

# convert annotation dictionary into list
annotation_dict_cluster = anno_dict_cluster.dict.tolist()

data.tl.annotation(
        annotation_information=annotation_dict_cluster, ## annotation dictionary
        cluster_res_key='leiden',
        res_key='anno_cluster_leiden' ## store annotation in "res_key" as a keyward
        )

# save as stereodata
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

data1.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=0.1)

# read text file into pandas DataFrame
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed/data/annotation_dict_subtype.txt", sep=" ")

# convert annotation dictionary into list
annotation_dict_cluster = anno_dict_cluster.dict.tolist()

data1.tl.annotation(
        annotation_information=annotation_dict_cluster, ## annotation dictionary
        cluster_res_key='leiden',
        res_key='anno_cluster_leiden' ## store annotation in "res_key" as a keyward
        )

# read stereodata
data = st.io.read_stereo_h5ad(
        file_path='/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.stereo.h5ad',
        use_raw=True,
        use_result=True,
        )

# draw a plot with 'tau' location
import pandas as pd
import matplotlib.pyplot as plt

# Extract cell names of tau from data1
tau = data1.cells.cell_name[data1.tl.result['anno_cluster_leiden']['group'].str.contains('cases')]

# Convert the cell names from 'data' into a pandas series
cell_names_series = pd.Series(data.cells.cell_name)

# Identify which cells in 'data' match the cell names in 'tau'
id_tau = cell_names_series.isin(tau)

# If match then 'red' otherwise 'grey'
data.cells['tau'] = ['red' if cell else 'gray' for cell in id_tau]
data.cells['tau'] = data.cells['tau'].astype('category')

# plot the cells by 'tau' (red)
fig, ax = plt.subplots(figsize=(10, 10))
plt.scatter(data.position[:,0], data.position[:,1], c=data.cells['tau'], s=10)
ax.invert_yaxis() # change the direction of y axis
