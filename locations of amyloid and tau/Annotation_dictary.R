library(SingleR)
library(anndata)  ## read_h5ad
library(zellkonverter) ## AnnData2SCE

## read h5ad files into AnnData objects
data = read_h5ad(file = '/work/aliu10/AD_Stereoseq_Project/processed_data/B01809C2/B01809C2.stereo.h5ad')
ref = read_h5ad(file = '/work/ygong/stereo_seq_public/MIT_AD.h5ad')

## convert AnnData objects into SingleCellExperiment objects
test = AnnData2SCE(data)
ref = AnnData2SCE(ref)

## Obtain annotation dictionary
## labels: cell types; assay.type.test: An integer scalar or string specifying the assay of test containing the relevant expression matrix.
annotation_dict = SingleR(test, ref, clusters=test$leiden, labels=ref$broad.cell.type, assay.type.test=1, assay.type.ref=1) 
