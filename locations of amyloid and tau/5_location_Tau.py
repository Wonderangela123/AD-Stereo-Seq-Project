import stereo as st
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

data = st.io.read_stereo_h5ad(
        file_path='/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.stereo.h5ad',
        use_raw=True,
        use_result=True,
        )

# Subset the cells annotated as "neuron"
id = data.tl.result['anno_cluster_leiden']['group'].str.contains('Ex|In')
id = id[id].index.tolist()

data1 = data
data1.cells.cell_name = data.cells.cell_name[id]
data1.exp_matrix = data.exp_matrix[id]
data1.position = data.position[id, :]

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

data1.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=0.5)

# After getting annotation dictionary
import pandas as pd

# read text file into pandas DataFrame
anno_dict_cluster = pd.read_csv("/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2_annotation_dict_subtype.txt", sep=" ")

# convert annotation dictionary into list
annotation_dict_cluster = anno_dict_cluster.dict.tolist()

data1.tl.annotation(
        annotation_information=annotation_dict_cluster, ## annotation dictionary
        cluster_res_key='leiden',
        res_key='anno_cluster_leiden' ## store annotation in "res_key" as a keyword
        )

#######################################################
# data will be covered, so we have to load data again.
data = st.io.read_stereo_h5ad(
        file_path='/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.stereo.h5ad',
        use_raw=True,
        use_result=True,
        )

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
plt.scatter(data.position[:,0], data.position[:,1], c=data.cells['tau'], s=5)

legend_elements = [Line2D([0], [0], marker='o', color='w', label='Tau', markersize=10, markerfacecolor='red'),
                   Line2D([0], [0], marker='o', color='w', label='Others', markersize=10, markerfacecolor='grey')]

# Add the legend to the plot
ax.legend(handles=legend_elements, fontsize=14)

ax.invert_yaxis() # change the direction of y axis

fig.savefig("/work/aliu10/AD_Stereoseq_Project/processed/data/B01809C2/GeneExpMatrix/B01809C2.png")
