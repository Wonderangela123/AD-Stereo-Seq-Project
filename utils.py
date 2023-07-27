import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import czifile
import matplotlib.pyplot as plt
from tifffile import imsave

# plot the gene expression on the figure
def feature_plot(stereo_data, genes, fig_size=(20,20), spot_size = 20):
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
def show_nucleus_image(file_path):
    # Load the .czi file
    img = czifile.imread(file_path)

    print("Shape of image: ", img.shape)
    image_2d = img[0, 0, :, :, 0]
    #image_2d = np.flip(image_2d, 1)
    # Show the image
    plt.imshow(image_2d, cmap='gray')
    plt.show()

