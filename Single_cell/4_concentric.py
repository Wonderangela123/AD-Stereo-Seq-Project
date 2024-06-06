import stereo as st
import scanpy as sc
import anndata as ad
import pandas as pd
import utils
import numpy as np
import matplotlib.pyplot as plt

def meta_info(stereo_data, meta_data):
    stereo_data.cells["leiden"] = meta_data["leiden"]
    stereo_data.cells["batch"] = meta_data["sample_ID"]
    leiden = pd.DataFrame({"group":stereo_data.cells["leiden"]})
    leiden.index = stereo_data.cell_names
    stereo_data.tl.result["leiden"] = leiden  

def add_position(stereo_data):
    x = stereo_data.position[:,0].copy().tolist()
    y = stereo_data.position[:,1].copy().tolist()
    stereo_data.cells["x"] = x
    stereo_data.cells["y"] = y
    return stereo_data

def is_in_circle(x, y, center, radius):
    dx = x - center[0]
    dy = y - center[1]
    distances = np.sqrt(dx*dx + dy*dy)
    return distances <= radius

def draw_circle(stereo_data, center, radius):
    inside_circle = is_in_circle(
        stereo_data.cells["x"], 
        stereo_data.cells["y"], 
        center, 
        radius
    )
    loc_each_dot = pd.DataFrame({"x": stereo_data.cells["x"], "y":stereo_data.cells["y"]})
    return loc_each_dot[inside_circle] 

def count_circle_overlaps(dots, location, radius):
    x = dots["x"]
    y = dots["y"]
    count = []
    indx = []
    for i in range(len(location)):
        center = [location.iloc[i,1],location.iloc[i,2]]
        indx_tmp = str(location.iloc[i,1])+"x"+str(location.iloc[i,2])
        tmp = is_in_circle(x, y, center, radius)
        count.append(tmp)
        indx.append(indx_tmp)
    count_all = pd.DataFrame(count)
    count_all.index = indx
    return count_all

def check_circle(stereo_data, filtered_positions):
    loc_each_dot = pd.DataFrame({"x": stereo_data.cells["x"], "y":stereo_data.cells["y"]})
    plt.figure(figsize=(10, 10))
    plt.scatter(loc_each_dot["x"], loc_each_dot["y"], color='blue')
    plt.scatter(filtered_positions["x"], filtered_positions["y"], color='green')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Spatial distribution of barcodes')
    plt.gca().invert_yaxis()
    plt.show()
    
def inner_outer_points(data, sample, inner_circle, outer_circle, location):
    globals()[sample] = data.tl.filter_by_clusters(cluster_res_key = "sample", groups = sample)
    
    
    globals()["location_{}".format(sample)] = location[location["sample_ID"]==sample]
    # for the gene expression
    inner_points_list_all = []
    outter_points_spec_list_all = []
    for i in range(globals()["location_{}".format(sample)].shape[0]):
        coor = globals()["location_{}".format(sample)].iloc[i,[1,2]].tolist()
        inner_points = draw_circle(globals()[sample], [coor[0], coor[1]], radius = inner_circle)
        inner_points["Barcode"] = inner_points.index.tolist()

        # for the outer cells
        outter_points = draw_circle(globals()[sample], [coor[0], coor[1]], radius = outer_circle)
        #print(outter_points)
        outter_points["Barcode"] = outter_points.index.tolist()
        outter_points_spec = outter_points[~outter_points["Barcode"].isin(inner_points["Barcode"])]
        #print(outter_points_spec)
        count = count_circle_overlaps(outter_points_spec, globals()["location_{}".format(sample)], inner_circle)
        #print(count.sum())
        outter_points_spec_list = count.sum()[count.sum()==0].index.tolist()


        # save the inner points
        inner_points_list = inner_points["Barcode"].tolist()

        inner_points_list_all += inner_points_list
        outter_points_spec_list_all += outter_points_spec_list    
    
    # make the elements unique
    inner_points_list_all = list(set(inner_points_list_all))
    outter_points_spec_list_all = list(set(outter_points_spec_list_all))
        
    return inner_points_list_all, outter_points_spec_list_all

def concentric_preprocess(data, levels_list):
      
    inner_all = []
    outer_all = []
    for sample in globals()[levels_list]:
        inner, outer = inner_outer_points(
            data, 
            sample, 
            inner_circle=500, 
            outer_circle=1500, 
            location = location
        )
        inner_all += inner
        outer_all += outer

    all_list = inner_all + outer_all
    data_conc = data.tl.filter_cells(cell_list=all_list, inplace = False)
    return data_conc, inner_all
def concentric_with_subtype(data_conc, celltype, inner_all):    
    meta = data_conc.cells.to_df()
    inner_ref = pd.DataFrame({"barcode":inner_all, "concentric": "inner"})
    inner_ref.set_index('barcode', inplace=True)
    meta['concentric'] = meta['Unnamed: 0'].astype(str).map(inner_ref['concentric']).fillna('faraway')

    sub = meta[meta["annotation"]== celltype]

    data_sub = data.tl.filter_cells(cell_list=sub.index.tolist(), inplace = False)
    data_sub.cells["concentric"] = sub["concentric"].to_list()
    utils.cell_t0_res_key(data_sub)
    merge_raw.obs["barcode"] = merge_raw.obs.index.to_list()
    merge_subset =merge_raw[merge_raw.obs["barcode"].isin(data_sub.cells.to_df().index.to_list())]
    merge_subset.obs["concentric"] = data_sub.cells["concentric"].to_list()
    sc.tl.rank_genes_groups(merge_subset, 'concentric', method='t-test')
    return merge_subset


# data list
control_list = ["B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6"]
moderate_list = ["B02008D2", "B02009F6"]
advanced_list = ["B01809C2", "C02248B5"]
severe_list = ["A02092E1", "B02008C6"]
case_list = ["A02092E1", "B02008C6", "B02008D2", "B02009F6", "B01809C2", "C02248B5"]
all_list = ["B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6",
           "A02092E1", "B02008C6", "B02008D2", "B02009F6", "B01809C2", "C02248B5"]

merge_raw = sc.read_h5ad("../Result/anndata/integrated_data/annotation_integrate.h5ad")

# load the location
location = pd.read_csv("../../50_bin_analysis/Result/Neuron_high_mt/high_mt_Ex_location.csv")
data = st.io.read_ann_h5ad(
                file_path="../Result/anndata/integrated_data/annotation_integrate.h5ad",
                spatial_key="spatial",
              )
meta_data = pd.read_csv("../Result/anndata/integrated_data/meta.csv")
for i in meta_data.columns:
    data.cells[i] = meta_data[i].to_list()

utils.cell_t0_res_key(data)

for i in ["control_list", "moderate_list", "advanced_list", "severe_list", "case_list", "all_list"]:
    control_conc, inner_all = concentric_preprocess(data, i)
    for cell in ['Opc', 'Undefined', 'End', 'Ex', 'Mic', 'Ast', 'Inh', 'Oli']:
        control_ex = concentric_with_subtype(control_conc, cell, inner_all)
        control_ex.write_h5ad("../Result/anndata/Concentric/{}/{}.h5ad".format(i, cell))

