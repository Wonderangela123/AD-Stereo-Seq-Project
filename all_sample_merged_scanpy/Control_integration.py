import scanpy as sc
import stereo as st
import anndata

# extract the expression data from the stereo data
control = ["B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6"]

raw_list = []
for i in control:
    data_path = "../processed_data/{}/GeneExpMatrix/{}.cellbin.gef".format(i, i)
    st.io.read_gef_info(data_path)
    data = st.io.read_gef(file_path=data_path, bin_type='cell_bins')
    data.tl.cal_qc()
    data.tl.filter_cells(
        min_gene=50,
        # maximum number of counts required for a cell to pass fitlering. 
        max_gene=700,
        min_n_genes_by_counts=100, 
        # maximum number of genes expressed required for a cell to pass filtering.
        max_n_genes_by_counts=1100, 
        # Remove cells that have too many mitochondrial genes expressed, without enough genes expressed, and out of count range.
        # maximum number of pct_counts_mt required for a cell to pass filtering.
        pct_counts_mt=20,
        # whether to inplace the previous data or return a new data.
        inplace=True
    )
    data.tl.raw_checkpoint()
    adata = st.io.stereo_to_anndata(data,flavor='scanpy')
    adata.obs["sample"] = i
    adata.obs["treatment"] = "control"
    raw_list.append(adata)

# Concatenate the AnnData objects
# use the outer to incoporate all the genes
merged_data = anndata.concat(raw_list, join = "outer")

sc.pl.highest_expr_genes(merged_data, n_top=20, )
sc.pp.normalize_total(merged_data, target_sum=1e4)

# log normalization
sc.pp.log1p(merged_data)

sc.pp.highly_variable_genes(merged_data, min_mean=0.0125, max_mean=3, min_disp=0.5)
merged_data.raw = merged_data

# select the highly variable genes
merged_data = merged_data[:, merged_data.var.highly_variable]

sc.pp.regress_out(merged_data, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(merged_data, max_value=10)

sc.tl.pca(merged_data, svd_solver='arpack')




















