# Read the Stereoseq data
import utils
import stereo as st


sample = "B01809C2"
data_path = "/work/ygong/stereo_seq_public/{}/GeneExpMatrix/{}.cellbin.gef".format(sample, sample)
st.io.read_gef_info(data_path)
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins')
