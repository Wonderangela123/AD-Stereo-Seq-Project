import stereo as st
from stereo.core.stereo_exp_data import AnnBasedStereoExpData
import scanpy as sc
import anndata
import os
import numpy as np

## Integration
samples = ["A02092E1", "B01809C2", "B02008C6", "B02008D2", "B02009F6", "C02248B5",
           "B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6"]

path_dir = "/work/aliu10/stereo_project/results/"

stereo_list = []
for i in samples:
    data_path = os.path.join(path_dir, "integrated_wilcoxon-{}.h5ad".format(i))
    # read data as stereobject 
    data = st.io.read_stereo_h5ad(
        file_path=data_path,
        use_raw=False, # use normalized data
        use_result=True,
        )
    data.cells['sample'] = i
    if i in samples[0:5]:
        data.cells['diagnosis'] = 'case'
    else: 
        data.cells['diagnosis'] = 'control'
    stereo_list.append(data)
        
merged_data = st.utils.data_helper.merge(stereo_list[0], stereo_list[1], stereo_list[2], stereo_list[3], stereo_list[4], stereo_list[5],
                                         stereo_list[6], stereo_list[7], stereo_list[8], stereo_list[9], stereo_list[10], stereo_list[11])

## Cell type annotation
annotation_dict = {
    '1':'excitatory neurons',
    '2':'oligodendrocytes.2',
    '3':'astrocytes.3',
    '4':'oligodendrocytes.4',
    '5':'astrocytes.5',
    '6':'oligodendrocytes.6',
    '7':'oligodendrocytes.7'
}

merged_data.tl.annotation(
               annotation_information=annotation_dict,
               cluster_res_key='spatial_leiden',
               res_key='anno_leiden'
               )

for i in range(12):
    stereo_list[i].tl.annotation(
                   annotation_information=annotation_dict,
                   cluster_res_key='spatial_leiden',
                   res_key='anno_leiden'
                   )

    conditions = [
        stereo_list[i].tl.result['anno_leiden']['group'].str.contains('astrocytes'),
        stereo_list[i].tl.result['anno_leiden']['group'].str.contains('oligodendrocytes')
    ]
    choices = ['astrocytes', 'oligodendrocytes']
    stereo_list[i].cells['anno_leiden_regroup'] = np.select(conditions, choices, default='excitatory neurons')

## Cell type annotation cluster scatter plots
# stereo_list[8].plt.cluster_scatter(res_key='anno_leiden_regroup')

## Extract neurons
id_neuron = merged_data.tl.result['anno_leiden']['group'].str.contains('excitatory neurons')
id_neuron = id_neuron[id_neuron].index.tolist()
data1 = merged_data
data1.cells.cell_name = merged_data.cells.cell_name[id_neuron]
data1.exp_matrix = merged_data.exp_matrix[id_neuron]
data1.position = merged_data.position[id_neuron, :]
data1.position_z = merged_data.position_z[id_neuron, :]

## Normalizing data should be skipped for cell subtyping because normalize total counts over all genes "per cell" such that each cell has the same total count after normalization in the last step. 
## PCA
data1.tl.pca(use_highly_genes=False, n_pcs=50, res_key='pca')
## UMAP
data1.tl.neighbors(pca_res_key='pca', n_pcs=50, res_key='neighbors')
data1.tl.spatial_neighbors(neighbors_res_key='neighbors', res_key='spatial_neighbors') # compute spatial neighbors
data1.tl.umap(pca_res_key='pca', neighbors_res_key='spatial_neighbors', res_key='spatial_umap')
## Clustering
data1.tl.leiden(neighbors_res_key='spatial_neighbors',res_key='spatial_leiden', resolution=1.2)

## Annotation
ref_file = '/work/aliu10/stereo_project/reference.h5ad'
ref = AnnBasedStereoExpData(ref_file)
data1.tl.single_r(ref_exp_data=ref, ref_use_col='diagnosis', cluster_res_key='spatial_leiden', res_key='annotation')
# data1.plt.cluster_scatter(res_key='annotation')

## Draw plots
# Extract cell names of tau from data1
tau = data1.cells.cell_name[data1.tl.result['annotation']['group'].str.contains('case')]

# Identify which cells in 'data' match the cell names in 'tau'
id_tau = np.isin(stereo_list[5].cells.cell_name, tau)
# If match then 'red' otherwise 'grey'
stereo_list[5].cells['tau'] = np.where(id_tau, 'red', 'grey')
stereo_list[5].cells['tau'] = stereo_list[5].cells['tau'].astype('category')
# plot the cells by 'tau' (red)
fig, ax = plt.subplots(figsize=(10, 10))
plt.scatter(stereo_list[5].position[:,0], stereo_list[5].position[:,1], c=stereo_list[5].cells['tau'], s=5)
legend_elements = [Line2D([0], [0], marker='o', color='w', label='Tau', markersize=10, markerfacecolor='red'),
                   Line2D([0], [0], marker='o', color='w', label='Others', markersize=10, markerfacecolor='grey')]
# Add the legend to the plot
ax.legend(handles=legend_elements, fontsize=14)
ax.invert_yaxis() # change the direction of y axis
fig.savefig(os.path.join(path_dir, "{}_tau_location.png".format(samples[5]))
