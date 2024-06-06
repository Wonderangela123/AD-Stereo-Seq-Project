def gene_mean_exp(table, data):
    list1 = table["genes"].tolist()
    list2 = data.gene_names.tolist()
    indices = [list2.index(item) for item in list1 if item in list2]
    table["mean_exp"] = data.exp_matrix[:, indices].mean(axis=0).tolist()[0]                  

def stereo_DE_table(data, res_key, gene_num = 40,mean_exp = False):
    import pandas as pd
    # select the group we are interested in
    first_list = list(data.tl.result["marker_genes"])
    second_list = data.tl.result[res_key].group.unique().tolist()
    selected_elements = [item for item in first_list for string in second_list if string in item]
    
    # save the DE table for each group into the DE_table_list 
    DE_table_list = []
    for name in selected_elements:
        DE_result = pd.DataFrame(data.tl.result["marker_genes"][name].copy())
        DE_result = DE_result[["pvalues_adj", "log2fc", "genes"]]
        DE_result["group"] = name
        
        if mean_exp == True:
            gene_mean_exp(DE_result, data)
        
        DE_result = DE_result.sort_values(by= ["log2fc"], ascending=[False])[0:gene_num]
        DE_table_list.append(DE_result)
    result_table = pd.concat(DE_table_list, axis=0, ignore_index=True)
        
    return result_table

# convert all the cells info into the res key
def cell_t0_res_key(data):
    import pandas as pd
    
    # add the position info into the cells
    data.cells["x"] = data.position[:,0]
    data.cells["y"] = data.position[:,1]
    
    for names in data.cells.to_df().columns:
        tmp = pd.DataFrame({"bins": data.cell_names, "group": data.cells[names]})
        data.tl.result[names] = tmp

# print the gene expression
def DE_table(data, gene_numbers):
    import pandas as pd
    result = data.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    table = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'logfoldchanges', 'pvals_adj']}).head(gene_numbers)
    # Assuming 'table' is your DataFrame
    # Separating gene names, values, and P-values
    gene_names_df = table.iloc[:, ::3]  # Select columns with gene names
    values_df = table.iloc[:, 1::3]     # Select columns with values
    p_values_df = table.iloc[:, 2::3]   # Select columns with P-values

    # Reset index to keep the index column for melting
    gene_names_df = gene_names_df.reset_index()
    values_df = values_df.reset_index()
    p_values_df = p_values_df.reset_index()

    # Melting the DataFrames
    melted_gene_names = pd.melt(gene_names_df, id_vars=['index'], var_name='Group_key', value_name='Gene_name')
    melted_values = pd.melt(values_df, id_vars=['index'], var_name='Group_key', value_name='Value')
    melted_p_values = pd.melt(p_values_df, id_vars=['index'], var_name='Group_key', value_name='P_adjusted')

    # Adjusting the 'Group_key' to have consistent group names in all DataFrames
    melted_gene_names['Group_key'] = melted_gene_names['Group_key'].str.replace('_n', '')
    melted_values['Group_key'] = melted_values['Group_key'].str.replace('_l', '')
    melted_p_values['Group_key'] = melted_p_values['Group_key'].str.replace('_p', '')

    # Combining the melted DataFrames
    long_table = pd.DataFrame({
        'index': melted_gene_names['index'],
        'Group_key': melted_gene_names['Group_key'],
        'Gene_name': melted_gene_names['Gene_name'],
        'Value': melted_values['Value'],
        'P_adjusted': melted_p_values['P_adjusted']
    })
    return long_table

def create_global_colormap(data, res_key, colormaps):
    import pandas as pd
    """
    Create a global color map for all unique values in the res_key column
    across all samples in the dataset.
    """
    unique_values = pd.unique(data.cells.to_df()[res_key])
    colormap = {uv: color for uv, color in zip(unique_values, colormaps)}
    return colormap

def scatter_layer_multi_single(ax, data, res_key, sample, colormap, angle_degrees = 0, h_flip = False, v_flip = False, save_path=None, title = True):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix
    meta = data.cells.to_df()
    meta[res_key] = data.tl.result[res_key]["group"].tolist()
    sub = meta[meta["sample"] == sample]

    sub['color'] = sub[res_key].map(colormap)

    
    # Rotation matrix
    angle_degrees = 0
    angle_radians = np.radians(angle_degrees)
    rotation_matrix = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians)],
        [np.sin(angle_radians), np.cos(angle_radians)]
    ])
    
    # Extracting original coordinates
    original_coordinates = sub[['x', 'y']].values
    
    # Rotating coordinates
    rotated_coordinates = original_coordinates.dot(rotation_matrix)
    scatter = ax.scatter(
        rotated_coordinates[:, 0], 
        rotated_coordinates[:, 1], 
        c=sub['color'], 
        cmap='viridis', 
        s = 4
    )
    
    if h_flip == True:
        current_xlim = ax.get_xlim()
        ax.set_xlim(current_xlim[::-1])
    
    if v_flip == True:
        current_ylim = ax.get_ylim()
        ax.set_ylim(current_ylim[::-1])
        
    ax.axis('off') 
    if title == True:
        ax.set_title(f'{sample}', fontweight = "bold", fontsize=20)
    

def scatter_gene_multi_single(ax,data, gene_name, sample,angle_degrees=0, h_flip = False, v_flip = False, save_path=None, title=True, side_bar = False):
    
    """
    Aim to plot the specific sample from the merge data.
    I have added the option for the figure rotation 
    """
    
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    # get location of the gene
    gene_index = list(data.gene_names).index(gene_name)
    gene_expression = data.exp_matrix[:, gene_index]

    if isinstance(gene_expression, csr_matrix):
        gene_expression = gene_expression.toarray().squeeze()

    meta = data.cells.to_df()
    meta["position_x"] = data.position[:,0]
    meta["position_y"] = data.position[:,1]
    meta["gene_exp"] = gene_expression
    sub = meta[meta["sample"]==sample]
    plt.figure(figsize=(5, 5))
    
    # check the max expression value
    exp_max = np.max(meta["gene_exp"])
    exp_min = np.min(meta["gene_exp"])
    
    # define the norm
    norm = plt.Normalize(exp_min, exp_max)
    
    angle_radians = np.radians(angle_degrees)

    # Rotation matrix
    rotation_matrix = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians)],
        [np.sin(angle_radians), np.cos(angle_radians)]
    ])

    # Extracting original coordinates
    original_coordinates = sub[['position_x', 'position_y']].values

    # Rotating coordinates
    rotated_coordinates = original_coordinates.dot(rotation_matrix)
    scatter = ax.scatter(
        rotated_coordinates[:, 0], 
        rotated_coordinates[:, 1], 
        c=sub['gene_exp'], 
        cmap='viridis', 
        s = 4,
        norm = norm
    )
    
    if h_flip == True:
        current_xlim = ax.get_xlim()
        ax.set_xlim(current_xlim[::-1])

    if v_flip == True:
        current_ylim = ax.get_ylim()
        ax.set_ylim(current_ylim[::-1])
        
    ax.axis('off') 
    
    if title == True:
        ax.set_title(f'{sample}', fontweight = "bold", fontsize=20)
    
    if side_bar == True:
        fig = ax.figure  # Get the figure to attach color bar
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.set_label('Gene Expression')
    #if save_path:
    #    plt.savefig(save_path, format= "png", dpi=600)
    
    #plt.show()


# scatter plot for the gene
def scatter_gene(data, gene_name, stereo = False):
    from scipy.sparse import csr_matrix
    import matplotlib.pyplot as plt
    
    if stereo == False:
        gene_index = list(data.var_names).index(gene_name)
        gene_expression = data.X[:, gene_index]
    else:
        gene_index = list(data.gene_names).index(gene_name)
        gene_expression = data.exp_matrix[:, gene_index]
        
    
    if isinstance(gene_expression, csr_matrix):
        gene_expression = gene_expression.toarray().squeeze()

    if stereo == False:
        plt.scatter(data.obs['x'], data.obs['y'], c=gene_expression, cmap='viridis', s = 1)
        plt.colorbar(label='Expression level of gene_name')
        plt.xlabel('X coordinate')
        plt.ylabel('Y coordinate')
        plt.title('Expression of gene_name')

        plt.gca().invert_xaxis()
        plt.show()
    else:
        plt.scatter(data.position[:,0], data.position[:,1], c=gene_expression, cmap='viridis', s = 4)
        plt.colorbar(label='Expression level of gene_name')
        plt.xlabel('X coordinate')
        plt.ylabel('Y coordinate')
        plt.title('Expression of gene_name')

        plt.gca().invert_xaxis()
        plt.show()

# only for the stereo_seq data
def multiple_scatter_gene(data, gene_name):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix
    
    gene_index = list(data.gene_names).index(gene_name)
    gene_expression = data.exp_matrix[:, gene_index]

    if isinstance(gene_expression, csr_matrix):
        gene_expression = gene_expression.toarray().squeeze()

    nrows, ncols = 3, 4
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 12))

    # Flatten the axes array for easy iteration
    axes = axes.flatten()

    # Plotting each sample in its own subplot
    for i, sample in enumerate(data.cells['sample'].unique()):
        ax = axes[i]
        sample_mask = data.cells['sample'] == sample
        sc = ax.scatter(data.position[sample_mask, 0], data.position[sample_mask, 1],
                        c=gene_expression[sample_mask], cmap='viridis', s=4)
        ax.set_title(sample)
        ax.invert_xaxis()  # Invert the x-axis if needed

    # Add a global colorbar for the figure
    plt.tight_layout()
    plt.show()    
    

# check layer
def scatter_layer(data):
    import matplotlib.pyplot as plt

    color = {"L1": "#1f77b4", "WM": "#ff7f0e", "L5": "#2ca02c", 
             "L4": "#d62728", "L2/3": "#9467bd", "L6": "#8c564b"}

    df = data.obs.copy()
    plt.figure(figsize=(10, 8))

    # Specify the order of layers for plotting
    layers_order = ["L1", "L2/3", "L4", "L5", "L6", "WM"]

    # Plot each layer in the specified order
    for layer in layers_order:
        layer_data = df[df['annotation'] == layer]
        plt.scatter(layer_data['x'], layer_data['y'], c=color[layer], label=layer, s=10)

    # Adding labels, title, and legend
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.title('Scatter Plot Colored by Annotation')
    plt.legend(title='Layer')

    # Invert the X-axis
    plt.gca().invert_xaxis()

    # Show the plot
    plt.show()

# plot the stereo seq
# just in case the original function doesn't work
def layer_scater_multi(data, res_key, dot_size):
    import matplotlib.pyplot as plt
    import pandas as pd

    # Sample barcode and clustering info
    df = pd.DataFrame(data.tl.result[res_key])

    # Extract the sample number after the hyphen
    df['sample'] = df['bins'].apply(lambda x: int(x.split('-')[-1]))
    df["sample"] = df["sample"].astype(str)

    # Sample position data
    position = data.position.copy()
    position = pd.DataFrame(position)
    position["sample"] = df["sample"]
    position["sample"] = position["sample"].astype(str)
    position["group"] = data.tl.result[res_key]["group"]

    # Colormap
    colormaps = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#FFFF33",
                 "#F781BF", "#999999", "#E5D8BD", "#B3CDE3", "#CCEBC5", "#FED9A6", "#FBB4AE",
                 "#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#FFED6F",
                 "#8DA0CB", "#E78AC3", "#E5C494", "#CCCCCC", "#FB9A99", "#E31A1C", "#CAB2D6",
                 "#6A3D9A", "#B15928"]

    # Create a color palette
    unique_groups = position['group'].unique()
    group_colors = {group: colormaps[i % len(colormaps)] for i, group in enumerate(unique_groups)}
    position['color'] = position['group'].map(group_colors)

    # Determine the number of rows and columns for the subplot grid
    n_samples = len(df['sample'].unique())
    n_cols = 3  # Adjust this as needed
    n_rows = n_samples // n_cols + (n_samples % n_cols > 0)

    # Create a figure with multiple subplots
    plt.figure(figsize=(15, 5 * n_rows))  # Adjust the figure size as needed
    for i, group in enumerate(sorted(df['sample'].unique()), 1):
        plt.subplot(n_rows, n_cols, i)
        group_positions = position[position["sample"] == group]
        plt.scatter(group_positions.iloc[:, 0], group_positions.iloc[:, 1], c=group_positions['color'], s = dot_size)
        plt.title(f'Scatter Plot for Group {group}')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')

    plt.tight_layout()
    plt.show()

# plot the violin plot
def VlnPlot(data, res_key, gene):
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    cluster = data.tl.result[res_key].copy()

    # gene expression
    gene_list = np.array(data.gene_names.copy())
    index = np.where(gene_list == gene)[0][0]

    cluster["gene_expression"] = data.exp_matrix[:,index].toarray()

    # Colormap
    colormaps = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#FFFF33",
                 "#F781BF", "#999999", "#E5D8BD", "#B3CDE3", "#CCEBC5", "#FED9A6", "#FBB4AE",
                 "#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#FFED6F",
                 "#8DA0CB", "#E78AC3", "#E5C494", "#CCCCCC", "#FB9A99", "#E31A1C", "#CAB2D6",
                 "#6A3D9A", "#B15928"]

    # Create a color palette
    unique_groups = cluster['group'].unique()
    group_colors = {group: colormaps[i % len(colormaps)] for i, group in enumerate(unique_groups)}

    plt.figure(figsize=(10, 6))  # Adjust the size as needed
    sns.violinplot(x='group', y='gene_expression', data=cluster, palette=group_colors, inner=None)

    plt.ylim(bottom=0)

    # Additional customization
    plt.title(f'Violin Plot of {gene} Expression by Group')
    plt.xlabel('Group')
    plt.ylabel('Gene Expression')
    plt.show()

# display the layer proportion        
def layer_proportion(data, res_key, group):
    
    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt

    # extract the meta_table
    meta_table = data.cells.to_df()
    
    # The colormap for the layers
    colormaps = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#FFFF33",
                 "#F781BF", "#999999", "#E5D8BD", "#B3CDE3", "#CCEBC5", "#FED9A6", "#FBB4AE",
                 "#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#FFED6F",
                 "#8DA0CB", "#E78AC3", "#E5C494", "#CCCCCC", "#FB9A99", "#E31A1C", "#CAB2D6",
                 "#6A3D9A", "#B15928"]

    if res_key == "annotation":
        # The desired layer order
        desired_order = ["L6", "L5", "L4", "L2/3", "L1"]

    # Assuming 'meta_table' is your DataFrame with the data
    # Filter the DataFrame to only include the desired annotations.
    filtered_df = meta_table[meta_table['annotation'].isin(desired_order)]

    # Create a proportion column where all values are equal to 1.
    filtered_df['proportion'] = 1

    # Pivot the DataFrame to get 'sample' as index and 'annotation' as columns.
    pivot_df = filtered_df.pivot_table(index='sample', columns='annotation', values='proportion', 
                                       aggfunc='sum', fill_value=0)

    # Ensure the columns are in the specified order and include missing columns as zeros.
    for annotation in desired_order:
        if annotation not in pivot_df:
            pivot_df[annotation] = 0
    pivot_df = pivot_df[desired_order]

    # Normalize the proportions so that each bar has the same height.
    pivot_df = pivot_df.div(pivot_df.sum(axis=1), axis=0)

    # Custom sample order provided by the user
    if group == "sample": 
        custom_order = ["D02175A4", "D02175A6", "B01809A3", "B01809A4", "B01806B5", "B01806B6", 
                        "B02008D2", "B02009F6", "B01809C2", "C02248B5", "A02092E1", "B02008C6"]
    
    if group == "diagnosis":
        custom_order = ["control", "case"]
    
    if group == "levels":
        custom_order = ["control", "moderate", "advanced", "severe"]
        
    # Reindex the DataFrame according to the custom order. Missing samples will be filled with zeros.
    pivot_df_custom_order = pivot_df.reindex(custom_order).fillna(0)

    # Generate the plot with the specified colors and custom sample order.
    plt.figure(figsize=(10, 7))
    bottom = np.zeros(len(pivot_df_custom_order))

    # Plot using the custom order.
    for idx, annotation in enumerate(desired_order):
        plt.bar(pivot_df_custom_order.index, pivot_df_custom_order[annotation], 
                bottom=bottom, color=colormaps[desired_order.index(annotation)], 
                edgecolor='white', label=annotation)
        bottom += pivot_df_custom_order[annotation].values

    # Rotate the x-axis labels vertically for better readability.
    plt.xticks(rotation=90)

    # Adding labels and title.
    plt.title('Normalized Proportion of Annotation Layers for Each Sample')
    plt.xlabel('Sample')
    plt.ylabel('Normalized Proportion')

    # Move the legend out of the plot.
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], title='Annotation', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    # Show the plot.
    plt.show()
































