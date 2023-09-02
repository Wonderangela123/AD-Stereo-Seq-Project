import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import czifile
import matplotlib.pyplot as plt
from tifffile import imsave
import statistics
import pyarrow.feather as feather
import datashader as ds
from datashader import transfer_functions as tf
import numpy as np
import utils
from PIL import Image
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
import holoviews as hv
import datashader as ds
from holoviews.operation.datashader import datashade, shade, dynspread, rasterize
import czifile

# plot the gene expression on the figure
def feature_plot(stereo_data, genes, fig_size=(20,20), spot_size = 20): # fig size is the hieght and width of the figure, spot size is the size for the spot.
    gene_to_plot = genes

    position = pd.DataFrame(stereo_data.position)
    expression_data = stereo_data.to_df()

    position.index = expression_data.index
    position.columns = ["x", "y"]
    x, y = position.iloc[:,0], position.iloc[:,1]

    gene_expression = expression_data[gene_to_plot]

    plt.figure(figsize=fig_size)
    plt.scatter(x, y, c= gene_expression, cmap = "viridis", alpha = 0.7, s = spot_size)
    plt.colorbar(label='Expression level')
    plt.show()

# plot the nucleus staining image
# to check the nucleus image in the Nuclei_Images folder
def show_nucleus_image(file_path):
    # Load the .czi file
    img = czifile.imread(file_path)

    print("Shape of image: ", img.shape)
    image_2d = img[0, 0, :, :, 0]
    #image_2d = np.flip(image_2d, 1)
    # Show the image
    plt.imshow(image_2d, cmap='gray')
    plt.show()

# count the detected genes for each spot
def genes_detected(data):
    expression_matrix = data.to_df()
    gene_identified = []
    for index_name in expression_matrix.index:
        tmp = len([i for i in expression_matrix.loc[index_name] if i >0])
        gene_identified.append(tmp)

    median_value = statistics.median(gene_identified)
    mean_value = statistics.mean(gene_identified)
    max_value = max(gene_identified)
    min_value = min(gene_identified)

    result = pd.DataFrame({
    "median_value": [median_value],
    "mean_value": [mean_value],
    "max_value": [max_value],
    "min_value": [min_value]
    })
    return result

####### to draw the Concentric circle for one spot ################
####### identify the cell barcodes and gene expression in an area #############
def is_in_circle(positions, center, radius):
    # given a position (e.g., a cell's position in data.plt.cells_plotting) and radius, check if the cells are in the
    """Return a boolean array indicating which rows are within a given circle."""
    dx = positions[:, 0] - center[0]  # difference between x coordinates and circle's x center
    dy = positions[:, 1] - center[1]  # difference between y coordinates and circle's y center
    distances = np.sqrt(dx*dx + dy*dy)  # calculate distances using Pythagorean theorem
    return distances <= radius  # return boolean array

def draw_circle(stereo_data, center, radius): # to define the circle and incorporate the cells
    inside_circle = is_in_circle(stereo_data.position, center, radius)
    filtered_positions = stereo_data.position[inside_circle]
    return filtered_positions

def check_circle(stereo_data, filtered_positions): # to check the location of the circle

    plt.figure(figsize=(30,30))
    # Plot all positions in blue
    plt.scatter(stereo_data.position[:, 0], stereo_data.position[:, 1], color='blue')

    # Overlay the filtered positions in green
    plt.scatter(filtered_positions[:, 0], filtered_positions[:, 1], color='green')
    

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Spatial distribution of barcodes')

    # Invert the y-axis to flip the plot vertically
    plt.gca().invert_yaxis()

    plt.show()

def find_barcode_expression(stereo_data, filtered_positions):# extrac the gene expression matrix of the barcode selected
    result_barcode_position = []
    for i in range(len(filtered_positions)):
        barcode = np.where((stereo_data.position[:,0] == filtered_positions[i][0]) & (stereo_data.position[:,1] == filtered_positions[i][1]))
        result_barcode_position.append(barcode)
    
    flat_list = [item[0][0] for item in result_barcode_position]
    
    
    expression_matrix = stereo_data.to_df()
    
    return expression_matrix.iloc[flat_list,:]

def pseudo_image(data):
    exp_matrix = data.array2sparse()
    row_sums = exp_matrix.sum(axis=1)
    location_info = data.position
    row_sums = np.array(row_sums)
    df = pd.DataFrame({
        'x': location_info[:, 0],
        'y': location_info[:, 1],
        'signal_strength': row_sums[:, 0].flatten()
    })
    
    # verticle flip
    height = df['y'].max()
    df['y'] = -df['y'] + height

    # Create points using HoloViews
    points = hv.Points(df, kdims=['x', 'y'], vdims=['signal_strength'])

    # Apply datashading operation and dynamic map
    result = dynspread(
        datashade(
            points, 
            element_type=hv.Image, 
            aggregator=ds.mean('signal_strength'), 
            cmap=["black", "green"])).opts(width=800, height=800, bgcolor="black")
    return result

def nuclei_image_import(data_path):
    
    # Load CZI into a numpy array
    img = czifile.imread(data_path)
    image_2d = img[0, 0, :, :, 0]

    # Assuming `image_2d` is your image array with shape (18576, 18684)
    # Create indices for x and y locations
    y_indices, x_indices = np.indices(image_2d.shape)

    # Flatten the arrays
    x_flat = x_indices.flatten()
    y_flat = y_indices.flatten()
    image_flat = image_2d.flatten()

    # Create DataFrame
    df = pd.DataFrame({
        'x': x_flat,
        'y': y_flat,
        'signal_strength': image_flat
    })
    # Create points using HoloViews
    points = hv.Points(df, kdims=['x', 'y'], vdims=['signal_strength'])
    
    # Negate x and y values and then shift them by the dimensions of the image
    df['x'] = -df['x'] + image_2d.shape[1]
    df['y'] = -df['y'] + image_2d.shape[0]

    
    # Apply datashading operation and dynamic map
    result = dynspread(
                datashade(
                    points, 
                    element_type=hv.Image, 
                    aggregator=ds.mean('signal_strength'), 
                    cmap=["black", "red"])).opts(width=800, height=800, bgcolor="black")
    return result




























