import stereo as st
import os
import pandas as pd
import gc
import utils
from importlib import reload
reload(utils)
import scanpy as sc
import anndata as ad

# check the data
merge = sc.read_h5ad("../Result/Anndata/New_integrate/merge.h5ad")

sc.tl.leiden(merge, resolution = 0.15)
merge.write_h5ad("../Result/Anndata/New_integrate/merge.h5ad")