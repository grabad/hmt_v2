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
    """
    lower_bound = df[column].quantile(percent)
    upper_bound = df[column].quantile(1 - percent)
    
    filtered_df = df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]
    
    return filtered_df

def filter_axial(df, padding=1):
    """
    Axial data often has many localizations at the edge of the axial range. 
    Remove points at the maximum and minimum z, with optional padding to remove points very close to the edges.
    """
    lower_bound = df["z [nm]"].min() + padding
    upper_bound = df["z [nm]"].max() - padding

    filtered_df = df[(df["z [nm]"] > lower_bound) & (df["z [nm]"] < upper_bound)]

    return filtered_df


def cluster_dbscan(loc_df, eps=50, min_samples=8, n_jobs=-1):
    """
    Cluster localization data using scikit-learn's DBSCAN algorithm, taking in epsilon and min_samples parameters.
    Clusters will be assigned a numeric ID and a random RGB color, with unclustered noise labeled as -1 and set gray.
    """
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

def plot_clusters_napari(me3_df, ac_df, scaling_factor=10):
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

def calc_area(locs):
    """
    Calculate the 2D area of a point cloud projected on the XY plane using its Convex Hull.
    """
    points = locs[["x [nm]", "y [nm]"]].values
    
    # A 2D convex hull requires at least 3 non-collinear points
    if len(points) < 3:
        return 0.0
        
    try:
        hull = ConvexHull(points)
        # Note: For a 2D ConvexHull in SciPy, .volume returns the area and .area returns the perimeter
        return hull.volume
    except Exception:
        # Handles errors if points are completely collinear/flat
        return 0.0
     

def create_radial_contours(binary_mask, num_bands=100, show_plots=False):
    """
    Breaks a non-circular binary mask into radial contours (concentric bands)
    using the Euclidean Distance Transform.
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
        contour_bands[inside_mask] = np.clip(np.digitize(distance_map[inside_mask], bins), 1, num_bands)
        
    if show_plots:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        
        im1 = axes[0].imshow(distance_map.T, origin='lower', cmap='viridis')
        axes[0].set_title('Distance Transform Map')
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
        
        im2 = axes[1].imshow(contour_bands.T, origin='lower', cmap='tab20')
        axes[1].set_title(f'Radial Contours ({num_bands} Bands)')
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
        
        plt.tight_layout()
        plt.show()
        
    return distance_map, contour_bands
