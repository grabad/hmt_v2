"""
HMT Functions for Pre-Processing SMLM images of post-translational histone modificaitons.
Data should already be separated into separate dataframes for different histone marks (i.e. H3K27me3 and H3K27ac).
"""
import napari
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import napari
from sklearn.cluster import DBSCAN
from scipy.ndimage import gaussian_filter, binary_fill_holes, label, distance_transform_edt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.qhull import QhullError
from scipy.optimize import minimize_scalar

def filter_quantile(df, column, percent=0.05):
    """
    Filter dataframe based by removing the top and bottom percent (default 0.05, leaving the 5th-95th percentile of the data).
    
    Returns:
        Pandas DataFrame: filtered_df
    """
    lower_bound = df[column].quantile(percent)
    upper_bound = df[column].quantile(1 - percent)
    
    filtered_df = df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]
    
    return filtered_df

def filter_axial(df, padding=1):
    """
    Axial data often has many localizations at the edge of the axial range. 
    Remove points at the maximum and minimum z, with optional padding to remove points very close to the edges.

    Returns:
        Pandas DataFrame: filtered_df
    """
    lower_bound = df["z [nm]"].min() + padding
    upper_bound = df["z [nm]"].max() - padding

    filtered_df = df[(df["z [nm]"] > lower_bound) & (df["z [nm]"] < upper_bound)]

    return filtered_df

def binarize_nucleus(me3_df, ac_df, thresh, bin_size=50, sigma=4.0, num_nuclei=1, show_plots=False):
    """
    Binarize nucleus based on intensity thresholding in 2D histogram.
    
    Returns:
        tuple: (binary_mask, me3_masked, ac_masked)
    """
    comb_df = pd.concat([me3_df, ac_df])
    x = comb_df["x [nm]"].values
    y = comb_df["y [nm]"].values
    
    x_bins = np.arange(x.min(), x.max() + bin_size, bin_size)
    y_bins = np.arange(y.min(), y.max() + bin_size, bin_size)
    density_map, _, _ = np.histogram2d(x, y, bins=(x_bins, y_bins))

    smoothed_map = gaussian_filter(density_map, sigma=sigma)
    binary_mask = smoothed_map > thresh
    binary_mask = binary_fill_holes(binary_mask)
    
    labeled_mask, num_features = label(binary_mask)
    if num_features > num_nuclei:
        component_sizes = np.bincount(labeled_mask.ravel())
        component_sizes[0] = 0  # Set background size to 0 so it's ignored
        
        largest_labels = np.argsort(component_sizes)[-num_nuclei:]
        binary_mask = np.isin(labeled_mask, largest_labels)
    
    # Map original coordinates to the binary mask grid to filter localizations
    x_idx = np.clip(np.digitize(x, x_bins) - 1, 0, binary_mask.shape[0] - 1)
    y_idx = np.clip(np.digitize(y, y_bins) - 1, 0, binary_mask.shape[1] - 1)
    in_mask = binary_mask[x_idx, y_idx]
    
    me3_masked = me3_df[in_mask[:len(me3_df)]]
    ac_masked = ac_df[in_mask[len(me3_df):]]

    if show_plots:
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        axes = axes.flatten()

        extent = [x_bins.min(), x_bins.max(), y_bins.min(), y_bins.max()]
        
        vmax = np.percentile(density_map, 99)
        im = axes[0].imshow(density_map.T, origin="lower", extent=extent, vmax=vmax)
        
        divider = make_axes_locatable(axes[0])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        axes[0].set_title("Density Map")
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        axes[0].set_aspect('equal', adjustable='box')

        axes[1].imshow(binary_mask.T, cmap="gray", origin="lower", extent=extent)
        axes[1].set_title("Binary Mask")
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        axes[1].set_aspect('equal', adjustable='box')

        axes[2].scatter(x, y, s=0.05, alpha=0.3, color='blue')
        axes[2].set_title("Raw Localizations")
        axes[2].set_xlim([x_bins.min(), x_bins.max()])
        axes[2].set_ylim([y_bins.min(), y_bins.max()])
        axes[2].set_xticks([])
        axes[2].set_yticks([])
        axes[2].set_aspect('equal', adjustable='box')

        axes[3].scatter(x[in_mask], y[in_mask], s=0.05, alpha=0.3, color='blue')
        axes[3].set_title("Masked Localizations")
        axes[3].set_xlim([x_bins.min(), x_bins.max()])
        axes[3].set_ylim([y_bins.min(), y_bins.max()])
        axes[3].set_xticks([])
        axes[3].set_yticks([])
        axes[3].set_aspect('equal', adjustable='box')

        plt.tight_layout()
        plt.show()

    return binary_mask, me3_masked, ac_masked

def sample_lower_densities(me3_df, ac_df, num_samples=10, random_state=42, show_plots=False, num_plots=4):
    """
    Subsamples the dataframes to simulate lower localization densities.
    Generates `num_samples` evenly spaced fractions from 1.0 down to > 0.
    
    Returns:
        list of tuples: Each tuple contains (me3_sampled, ac_sampled, fraction)
    """
    fractions = np.linspace(1.0, 1.0/num_samples, num_samples)
    
    sampled_data = []
    for frac in fractions:
        me3_samp = me3_df.sample(frac=frac, random_state=random_state) if frac < 1.0 else me3_df.copy()
        ac_samp = ac_df.sample(frac=frac, random_state=random_state+1) if frac < 1.0 else ac_df.copy() 
        sampled_data.append((me3_samp, ac_samp, frac))
        
    if show_plots:
        num_plots = min(num_plots, len(sampled_data))
        # Grab evenly spaced indices
        indices_to_plot = np.linspace(0, len(sampled_data) - 1, num_plots, dtype=int)
        
        cols = min(4, num_plots)
        rows = int(np.ceil(num_plots / cols))
        
        fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
        if rows * cols == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
            
        for i, ax in enumerate(axes):
            if i < num_plots:
                idx = indices_to_plot[i]
                me3_samp, ac_samp, frac = sampled_data[idx]
                
                ax.scatter(me3_samp["x [nm]"], me3_samp["y [nm]"], s=0.05, alpha=0.3, color='blue')
                ax.scatter(ac_samp["x [nm]"], ac_samp["y [nm]"], s=0.05, alpha=0.3, color='blue')
                
                ax.set_title(f"Density Fraction: {frac:.2f}")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_aspect('equal', adjustable='box')
            else:
                ax.axis('off') # Hide empty subplots
            
        plt.tight_layout()
        plt.show()

    return sampled_data, fractions