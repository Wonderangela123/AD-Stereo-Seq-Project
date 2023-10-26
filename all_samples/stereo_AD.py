import scanpy as sc
import stereo as st
import anndata
import os

## integrate samples
samples = ["A02092E1", "B01809C2", "B02008C6", "B02008D2", "B02009F6", "C02248B5",
           "B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6"]

path_dir = "/work/aliu10/AD_Stereoseq_Project/processed_data"

stereo_list = []
anndata_list = []
for i in samples:
    data_path = os.path.join(path_dir, "cell_bin_combined-{}.h5ad".format(i))
    # read data as stereobject 
    data = st.io.read_stereo_h5ad(
        file_path=data_path,
        use_raw=True,
        use_result=True,
        )
    stereo_list.append(data)
    # save data as anndata
    adata = st.io.stereo_to_anndata(data,flavor='scanpy')
    adata.obs["sample"] = i
    if i in samples[0:5]:
        adata.obs["treatment"] = "case"    # define group as "case"
    else:
        adata.obs["treatment"] = "control" # define group as "control"
    anndata_list.append(adata)

# merge data    
merged_data = st.utils.data_helper.merge(stereo_list[0], stereo_list[1], stereo_list[2], stereo_list[3], stereo_list[4], stereo_list[5],
                                         stereo_list[6], stereo_list[7], stereo_list[8], stereo_list[9], stereo_list[10], stereo_list[11])

# normalize data
merged_data.tl.normalize_total()
merged_data.tl.log1p()

# PCA
merged_data.tl.pca(use_highly_genes=False, n_pcs=50, res_key='pca')

# integrate
merged_data.tl.batches_integrate(pca_res_key='pca', res_key='pca_integrated')

# perform UMAP to check the consistency of different samples (check integration effects) 
merged_data.tl.neighbors(pca_res_key='pca_integrated', n_pcs=50, res_key='neighbors_integrated')
merged_data.tl.spatial_neighbors(neighbors_res_key='neighbors_integrated', res_key='spatial_neighbors') # compute spatial neighbors
merged_data.tl.umap(pca_res_key='pca_integrated', neighbors_res_key='spatial_neighbors', res_key='umap_integrated')

# UMAP plot
# merged_data.plt.batches_umap(res_key='umap_integrated', title ='UMAP between samples')

merged_data.tl.leiden(neighbors_res_key='spatial_neighbors',res_key='spatial_leiden')

merged_data.tl.raw_checkpoint()

# find marker genes
merged_data.tl.find_marker_genes(cluster_res_key='spatial_leiden', method='wilcoxon_test', use_highly_genes=False, use_raw=True, output = "/work/aliu10/stereo_project/results/marker_genes_wilcoxon.txt")
merged_data.tl.filter_marker_genes(marker_genes_res_key='marker_genes',
                                   min_fold_change=1,
                                #    min_in_group_fraction=0.25,
                                #    max_out_group_fraction=0.5,
                                   res_key='marker_genes_filtered',
                                   output = "/work/aliu10/stereo_project/results/marker_genes_filtered_wilcoxon.txt"
)

st.io.write_h5ad(
        merged_data,
        use_raw=True,
        use_result=True,
        key_record=None,
        output="/work/aliu10/stereo_project/results/integrated_wilcoxon.h5ad")
