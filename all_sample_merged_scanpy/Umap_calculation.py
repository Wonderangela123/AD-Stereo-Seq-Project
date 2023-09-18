import scanpy as sc
import stereo as st
import anndata

merged_data = sc.read_h5ad("Result/merged_data.h5ad")

# clustering analysis
sc.pp.neighbors(merged_data, n_neighbors=10, n_pcs=40)
sc.tl.leiden(merged_data)

# Embedding the neighborhood graph
sc.tl.paga(merged_data)
sc.pl.paga(merged_data, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(merged_data, init_pos='paga')

sc.pl.umap(merged_data, color = ["sample"])

sc.pl.umap(merged_data, color = ["sample"], save= ".pdf")















