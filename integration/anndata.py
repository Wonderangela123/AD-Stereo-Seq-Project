# Read the Stereo-seq data
import stereo as st

samples_list = ["A02092E1", "B01806B5", "B01806B6", "B01809A3", "B01809A4", "B01809C2", "B02008C6", "B02008D2", "B02009F6", "C02248B5", "D02175A4", "D02175A6"]
                
for sample in samples_list:
    print(sample)
    data_path = "/work/aliu10/AD_Stereoseq_Project/processed/{}/GeneExpMatrix/{}.cellbin.gef".format(sample, sample)
    st.io.read_gef_info(data_path)
    data = st.io.read_gef(file_path=data_path, bin_type='cell_bins')
    data.tl.raw_checkpoint()
    st.io.stereo_to_anndata(data,
                            flavor='seurat',
                            output='/work/aliu10/AD_Stereoseq_Project/processed/{}/GeneExpMatrix/{}.anndata.h5ad'.format(sample, sample))
