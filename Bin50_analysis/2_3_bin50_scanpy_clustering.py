#import stereo as st
import os
import gc
import numpy as np
import pandas as pd
#from stereo.utils.data_helper import split
import numpy as np
import scanpy as sc

merge = sc.read("../Result/Anndata/Integrate/stereo_seq_integrate.h5ad")

# perform the clustering 
sc.tl.leiden(merge, resolution = 0.1)

# save h5ad data
merge.write_h5ad("../Result/Anndata/Integrate/stereo_seq_integrate_clustering.h5ad")