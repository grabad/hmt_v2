import napari
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import napari
from sklearn.cluster import DBSCAN
from scipy.ndimage import gaussian_filter, binary_fill_holes, label, distance_transform_edt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial import ConvexHull

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


def cluster_dbscan(loc_df, eps=50, min_samples=8, n_jobs=-1):
    """
    Cluster localization data using scikit-learn's DBSCAN algorithm, taking in epsilon and min_samples parameters.
    Clusters will be assigned a numeric ID and a random RGB color, with unclustered noise labeled as -1 and set gray.

    Returns:
        Pandas DataFrame: loc_df (with cluster_label and cluster_color columns)
    """
    loc_df = loc_df.copy()

    loc_np = loc_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
    loc_clust = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=n_jobs).fit_predict(loc_np)    
    loc_df["cluster_label"] = loc_clust
    
    
    rng = np.random.default_rng(0)  # reproducible colors
    alpha = 0.5

    loc_unique_clusters = np.unique(loc_clust[loc_clust != -1])
    loc_cluster_colors = {c: np.append(rng.random(3), 1) for c in loc_unique_clusters}
    loc_colors = []

    for c in loc_clust:
        if c == -1:
            loc_colors.append((0.7, 0.7, 0.7, alpha))
        else:
            loc_colors.append(loc_cluster_colors[c])
            
    loc_df["cluster_color"] = loc_colors

    return loc_df

def plot_clusters_napari(me3_df, ac_df, scaling_factor=10):
    me3_df = me3_df.copy()
    ac_df = ac_df.copy()
    
    me3_df["plot_size"] = me3_df['sigmax [nm]'].to_numpy()/110*scaling_factor
    ac_df["plot_size"] = ac_df['sigmax [nm]'].to_numpy()/110*scaling_factor

    coords = ["x [nm]", "y [nm]", "z [nm]"]
    me3_border = [0, 0.7, 0, 1]
    ac_border = [0.7, 0, 0, 1]

    me3_nano = me3_df[me3_df["cluster_label"] != -1]
    me3_noise = me3_df[me3_df["cluster_label"] == -1]
    ac_nano = ac_df[ac_df["cluster_label"] != -1]
    ac_noise = ac_df[ac_df["cluster_label"] == -1]

    viewer = napari.Viewer()
    viewer.add_points(
        data=ac_noise[coords],
        features=ac_noise,
        face_color=ac_noise["cluster_color"].tolist(),
        border_color=ac_border, 
        border_width=0.01, 
        size=ac_noise["plot_size"],    
        name="H3K27ac Noise",
        out_of_slice_display=True
    )           
    viewer.add_points(
        data=me3_noise[coords],
        features=me3_noise,
        face_color=me3_noise["cluster_color"].tolist(), 
        border_color=me3_border, 
        border_width=0.01, 
        size=me3_noise["plot_size"],   
        name="H3K27me3 Noise",
        out_of_slice_display=True
    )                     
    viewer.add_points(
        data=ac_nano[coords],
        features=ac_nano,
        face_color=ac_nano["cluster_color"].tolist(), 
        border_color=ac_border, 
        border_width=0.01, 
        size=ac_nano["plot_size"],    
        name="H3K27ac Nanodomains",
        out_of_slice_display=True
    )           
    viewer.add_points(
        data=me3_nano[coords],
        features=me3_nano,
        face_color=me3_nano["cluster_color"].tolist(),
        border_color=me3_border, 
        border_width=0.01, 
        size=me3_nano["plot_size"],   
        name="H3K27me3 Nanodomains",
        out_of_slice_display=True
    )

    viewer.dims.order = (2, 1, 0)
    viewer.dims.ndisplay = 3
    

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

def calc_nanodomain_size(locs):
    """
    Calculate and return the size of the nanodomain using 4 methods:
    1) Area of xy convex hull
    2) Volume of xyz convex hull
    3) Average dimensions of xy bounding box
    4) Average dimensions of xyz bounding box
    
    Returns:
        tuple: (hull_size_2D, hull_size_3D, bb_size_2D, bb_size_3D)
    """
    points_2D = locs[["x [nm]", "y [nm]"]].values
    points_3D = locs[["x [nm]", "y [nm]", "z [nm]"]].values
    
    if len(points_2D) < 4:
        return (0.0, 0.0, 0.0, 0.0)
        
    try:
        hull_2D = ConvexHull(points_2D)
        hull_3D = ConvexHull(points_3D)
        
        hull_size_2D = hull_2D.volume
        hull_size_3D = hull_3D.volume
        
        bb_size_2D = ((points_2D.max(axis=0) - points_2D.min(axis=0)) + 
                      (points_2D.max(axis=1) - points_2D.min(axis=1)))/2
        
        bb_size_3D = ((points_3D.max(axis=0) - points_3D.min(axis=0)) + 
                      (points_3D.max(axis=1) - points_3D.min(axis=1)) + 
                      (points_3D.max(axis=2) - points_3D.min(axis=2)))/3

        return (hull_size_2D, hull_size_3D, bb_size_2D, bb_size_3D)
    
    except Exception:
        return (0.0, 0.0, 0.0, 0.0)
    
def calc_loc_density(locs):
    """
    Calculate localization density for an individual nanodomain (or any point cloud passed in)
    
    Returns:
        tuple: (density_2D, density_3D)
    """
    points_2D = locs[["x [nm]", "y [nm]"]].values
    points_3D = locs[["x [nm]", "y [nm]", "z [nm]"]].values

    if len(points_2D) < 4:
        return (0.0, 0.0)
        
    try:
        density_2D = len(points_2D)/ConvexHull(points_2D).volume
        density_3D = len(points_3D)/ConvexHull(points_3D).volume

        return (density_2D, density_3D)
    
    except Exception:
        return (0.0, 0.0)

def create_radial_contours(binary_mask, num_bands=100, show_plots=False):
    """
    Breaks a non-circular binary mask into radial contours (concentric bands)
    using the Euclidean Distance Transform.
    
    Returns:
        tuple: (distance_map, contour_bands)
    """
    # Calculate the distance of each pixel inside the mask to the nearest background pixel
    distance_map = distance_transform_edt(binary_mask)
    
    # Create discrete contour bands
    contour_bands = np.zeros_like(distance_map, dtype=int)
    max_dist = distance_map.max()
    
    if max_dist > 0:
        # Bin the distances into the requested number of bands
        # Band 1 is the outermost boundary, Band `num_bands` is the innermost core
        bins = np.linspace(0, max_dist, num_bands + 1)
        
        # Only assign bands to pixels inside the mask
        inside_mask = binary_mask > 0
        outside_mask = binary_mask <= 0
        
        contour_bands[inside_mask] = np.clip(np.digitize(distance_map[inside_mask], bins), 1, num_bands)
        contour_bands[inside_mask] = 100 - contour_bands[inside_mask]
        contour_bands[outside_mask] = contour_bands[outside_mask] - 1
        
    if show_plots:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        
        im1 = axes[0].imshow(distance_map.T, origin='lower', cmap='viridis')
        axes[0].set_title('Distance Transform Map')
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
        
        cmap = plt.get_cmap('tab20_r').copy()
        cmap.set_under('white')
        im2 = axes[1].imshow(contour_bands.T, origin='lower', cmap=cmap, vmax=num_bands, vmin=0)
        axes[1].set_title(f'Radial Contours ({num_bands} Bands)')
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
        
        plt.tight_layout()
        plt.show()
        
    return distance_map, contour_bands

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

    return sampled_data 

"""
TODO: 
SNCR Function (KD Tree most likely)
Density Correction
Ripley's K implementation
Final outputs
"""