# Read the Stereoseq data through python and write the gene expression matrix into feather file and gene name into csv 
import utils
import pandas as pd
import stereo as st
import pyarrow.feather


sample = "B01809C2"
data_path = "/work/ygong/stereo_seq_public/{}/GeneExpMatrix/{}.cellbin.gef".format(sample, sample)
st.io.read_gef_info(data_path)
data1 = st.io.read_gef(file_path=data_path, bin_type='cell_bins')
df = data1.to_df()
df = df.T
pyarrow.feather.write_feather(df, "/work/aliu10/AD_Stereoseq_Project/processed_data/{}_matrix.feather".format(sample))
gene_list = pd.DataFrame(df.index)
gene_list.to_csv("/work/aliu10/AD_Stereoseq_Project/processed_data/{}_gene_list.csv".format(sample))
