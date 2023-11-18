from stereo import image as im
import os
from stereo.tools.cell_cut import CellCut
from stereo.tools.cell_correct import cell_correct

def cell_segmentation(sample):
    model_dir = '../Deep_Cell_Model/Deep_Cell_Model/'
    for i in os.listdir("../../../processed_data/{}/Register/".format(sample)):
        if "regist.tif" in i:
            img_path = '../../../processed_data/{}/Register/{}'.format(sample, i)
    out_path = 'Result'
    tissue_seg_model_path = "../weight_tissue_cut_tool_220304.hdf5"
    tissue_seg_method = 1

    # Cell segmentation
    im.cell_seg_deepcell(
        model_dir, 
        img_path, 
        out_path, 
        tissue_seg_model_path=tissue_seg_model_path, 
        tissue_seg_method=1
    )

    # process the mask result for the cell bin segmentation
    cgef_out_dir = "Result"
    bgef_path = "../../../processed_data/{}/GeneExpMatrix/{}.gef".format(sample, sample)
    for i in os.listdir("Result"):
        if "mask.tif" in i:
            mask_path = "Result/{}".format(i)

    cc = CellCut(cgef_out_dir=cgef_out_dir)
    out_path = cc.cell_cut(bgef_path=bgef_path, mask_path=mask_path)
    































if __name__ == '__main__':
    cell_segmentation('D02175A4')
