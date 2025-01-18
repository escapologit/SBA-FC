import pandas as pd
import numpy as np
from nilearn import plotting, datasets
from nilearn import image
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
import matplotlib as mpl

# Load CSV files
p_values = pd.read_csv("out_yeo/group_p_values_run-01_lsd_vs_plcb.csv", header=None).values
t_values = pd.read_csv("out_yeo/group_t_values_run-01_lsd_vs_plcb.csv", header=None).values

# Load the Yeo 7-network parcellation (thick version)
yeo_atlas = datasets.fetch_atlas_yeo_2011()
yeo_atlas_7 = yeo_atlas.thick_7  # You can use either thin_7 or thick_7

region_names = ["Visual", "Somatomotor", "Dorsal Attention", "Ventral Attention", "Limbic", "Frontoparietal Control", "Default Mode Network"]

# Identify significant pairs (p-value <= 0.05) excluding NaNs
significant_pairs = []
num_regions = p_values.shape[0]

for i in range(num_regions):
    for j in range(i + 1, num_regions):  # Only check the upper triangle of the matrix
        if p_values[i, j] <= 0.05 and not np.isnan(p_values[i, j]):
            significant_pairs.append((i, j))

# Print the significant pairs
print("Significant pairs (ROIs):")
for pair in significant_pairs:
    print(f"ROI {pair[0]+1} <-> ROI {pair[1]+1}")

# Function to plot each significant pair with distinct colors and legend
def plot_significant_pairs_with_legend(yeo_atlas_filename, significant_pairs):
    # Load the Yeo atlas as an image (with 7 regions)
    atlas_img = image.load_img(yeo_atlas_filename)
    
    # Define 7 distinct colors for the 7 Yeo networks
    yeo_colors = [
        (0.0, 0.0, 1.0),  # Blue
        (0.0, 1.0, 0.0),  # Green
        (1.0, 0.0, 0.0),  # Red
        (1.0, 1.0, 0.0),  # Yellow
        (0.0, 1.0, 1.0),  # Cyan
        (1.0, 0.0, 1.0),  # Magenta
        (1.0, 0.647, 0.0)  # Orange
    ]
    
    # Create a ListedColormap with 7 colors
    cmap = mpl.colors.ListedColormap(yeo_colors)
    
    # Set up the color map for plotting
    norm = mpl.colors.BoundaryNorm(boundaries=np.arange(1, 9), ncolors=7)
    
    # Plot each significant ROI pair
    for idx, (roi1, roi2) in enumerate(significant_pairs):
        # Create a mask for the current significant ROI pair (roi1 and roi2)
        roi_mask = np.zeros_like(atlas_img.get_fdata(), dtype=int)
        
        # Mask for ROI 1
        roi_mask[atlas_img.get_fdata() == (roi1 + 1)] = (roi1 % 7) + 1  # Assign a network (1-7)
        # Mask for ROI 2
        roi_mask[atlas_img.get_fdata() == (roi2 + 1)] = (roi2 % 7) + 1  # Assign a network (1-7)
        
        roi_mask_img = image.new_img_like(atlas_img, roi_mask.astype(np.int32))  # Ensure int32 type

        # Plot the brain surface with the highlighted ROIs for the current pair
        plotting.plot_roi(roi_mask_img, title=f'Significant ROI Pair: {roi1 + 1}, {roi2 + 1}',
                  draw_cross=False, colorbar=False, cmap=cmap, display_mode='ortho', vmin=1, vmax=7)

        # Define a new axes for the color bar
        cbar_ax = plt.gca().inset_axes([0.9, 0, 0.05, 0.8])  # [x, y, width, height]

        # Add the custom color bar with 7 discrete bars (one for each region)
        cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cbar_ax, orientation='vertical')
        cbar.set_ticks(np.arange(1, 8))  # Set the ticks for each region (1-7)
        cbar.set_ticklabels(['1', '2', '3', '4', '5', '6', '7'])  # Custom labels for each region
        cbar.set_label('Yeo Network Regions')

        
        # Create a legend for the current pair
        handles = []
        labels = []
        
        # Get the color from the colormap for the current ROIs
        color_roi1 = yeo_colors[roi1 % 7]
        color_roi2 = yeo_colors[roi2 % 7]
        
        # Add legend entries with the colors matching the networks
        handle1 = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_roi1, markersize=10)
        handle2 = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_roi2, markersize=10)
        handles.extend([handle1, handle2])
        labels.extend([f"{region_names[roi1]}", f"{region_names[roi2]}"])
        
        # Add the legend to the plot
        plt.legend(handles=handles, labels=labels, loc='upper right', bbox_to_anchor=(0.8, 1.3)) #title="Significant ROIs") # bbox is x,y
        
        # Show the plot
        plt.savefig(f"out_yeo/run-01_{roi1+1}vs{roi2+1}.png", dpi=300, bbox_inches='tight')
        plotting.show()

# Plot each significant pair with distinct colors and legend
plot_significant_pairs_with_legend(yeo_atlas_7, significant_pairs)