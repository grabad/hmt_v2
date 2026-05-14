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


def cluster_dbscan(loc_df, eps=50, min_samples=8, n_jobs=-1):
    """
    Cluster localization data using scikit-learn's DBSCAN algorithm, taking in epsilon and min_samples parameters.
    Clusters will be assigned a numeric ID and a random RGB color, with unclustered noise labeled as -1 and set gray.

    Returns:
        Pandas DataFrame: loc_df (with cluster_label and cluster_color columns)
    """
    loc_df = loc_df.copy()

    loc_np = loc_df[["x [nm]", "y [nm]"]].to_numpy()
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
    
    if not ac_noise.empty:
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
    
    if not me3_noise.empty:   
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
    
    if not ac_nano.empty:                     
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
        
    if not me3_nano.empty:           
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

def remove_nucleoli():

    return 0

def calc_nanodomain_size(locs):
    """
    Calculate and return the size of the nanodomain using 4 methods.
    Returns a dictionary of the results.
    """
    points_2D = locs[["x [nm]", "y [nm]"]].values
    points_3D = locs[["x [nm]", "y [nm]", "z [nm]"]].values
    
    # Define a default empty dictionary for failed/small clusters
    empty_res = {'hull_2D': 0.0, 'hull_3D': 0.0, 'bb_2D': 0.0, 'bb_3D': 0.0}
    
    if len(points_2D) < 4:
        return empty_res
        
    try:
        hull_2D = ConvexHull(points_2D)
        hull_3D = ConvexHull(points_3D)
        
        hull_size_2D = hull_2D.volume
        hull_size_3D = hull_3D.volume
        
        bb_size_2D = (points_2D.max(axis=0) - points_2D.min(axis=0)).mean()
        bb_size_3D = (points_3D.max(axis=0) - points_3D.min(axis=0)).mean()

        # Return as a dictionary mapping the option names to the values
        return {
            'hull_2D': hull_size_2D,
            'hull_3D': hull_size_3D,
            'bb_2D': bb_size_2D,
            'bb_3D': bb_size_3D
        }
    
    except Exception:
        return empty_res
    
def calc_loc_density(hull_3D):
    """
    Calculate 2D and 3D localization density from a pre-computed 3D convex hull.
    
    Returns:
        tuple: (density_2D, density_3D)
    """
    points_3D = hull_3D.points
    points_2D = points_3D[:, :2]

    if len(points_2D) < 4:
        return (0.0, 0.0)
        
    try:
        # The 3D hull is already calculated. Handle case where volume is 0.
        density_3D = len(points_3D)/hull_3D.volume if hull_3D.volume > 0 else 0.0
        
        # We still need to calculate the 2D hull for the 2D density
        hull_2D = ConvexHull(points_2D)
        density_2D = len(points_2D)/hull_2D.volume if hull_2D.volume > 0 else 0.0

        return (density_2D, density_3D)
    
    except QhullError:
        return (0.0, 0.0)
    
def calculate_nanodomain_characteristics(locs):
    """
    Calculates characteristics of nanodomains found via DBSCAN clustering, stores information into new dataframe

    Returns:
        Pandas DataFrame: cluster_df
    """

    # 1. Correctly initialize with column names and use a list to collect rows.
    column_names = [
        "label",
        "cluster_color",
        "centroid",
        "vertices",
        "hull_size_2D",
        "hull_size_3D",
        "bb_size_2D",
        "bb_size_3D",
        "density_2D",
        "density_3D",
    ]
    cluster_rows = []

    num_clusters = locs['cluster_label'].max()
    for label in range(num_clusters):
        cluster_data = locs[locs['cluster_label'] == label]

        # 3. Add checks to prevent errors on small or degenerate clusters.
        if len(cluster_data) < 4:  # Need at least 4 points for a 3D hull.
            continue
        
        try:
            nanodomain_hull = ConvexHull(cluster_data[['x [nm]', 'y [nm]', 'z [nm]']].values)
        except QhullError:
            print(f"Warning: Could not compute convex hull for cluster {label}. Skipping.")
            continue
        
        size_measures = calc_nanodomain_size(nanodomain_hull)
        density_2D, density_3D = calc_loc_density(nanodomain_hull)
        
        # 4. Create a dictionary for the new row. This avoids the length mismatch error.
        new_row_dict = {
            "label": label,
            "cluster_color": cluster_data['cluster_color'].iloc[0],
            "centroid": np.mean(nanodomain_hull.points[nanodomain_hull.vertices], axis=0),
            "vertices": nanodomain_hull.points[nanodomain_hull.vertices],
            "hull_size_2D": size_measures.hull_2D,
            "hull_size_3D": size_measures.hull_3D,
            "bb_size_2D": size_measures.bb_2D,
            "bb_size_3D": size_measures.bb_2D,
            "density_2D": density_2D,
            "density_3D": density_3D}
        
        cluster_rows.append(new_row_dict)

    # 5. Create the final DataFrame from the list of rows in one efficient operation.
    cluster_df = pd.DataFrame(cluster_rows, columns=column_names)
    return cluster_df

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

    return sampled_data, fractions

def build_probability_map(distance_map, radial_density_profile):
    """
    Maps the experimental 1D density profile to the pre-calculated 3D distance shells
    to create a 3D probability map.
    """
    prob_map = np.zeros_like(distance_map, dtype=float)
    max_dist = int(np.ceil(distance_map.max()))
    
    # Apply the experimental density profile to the distance shells
    for d in range(1, max_dist + 1):
        # Find pixels that fall within this distance band (e.g., between 1 and 2 pixels away from the edge)
        mask = (distance_map >= d) & (distance_map < d + 1)
        
        # Ensure we don't index out of bounds if the nucleus is exceptionally large
        prof_idx = min(d - 1, len(radial_density_profile) - 1)
        prob_map[mask] = radial_density_profile[prof_idx]
        
    return prob_map

def plant_seeds(prob_map, real_z_coords, x_min, y_min, px_size=50.0, scaling_factor=1.0):
    # 1. Roll randoms
    random_rolls = np.random.rand(*prob_map.shape)
    thresholds = prob_map * scaling_factor
    seed_indices = np.argwhere(random_rolls < thresholds)
    
    jitter = np.random.rand(*seed_indices.shape)
    seed_coords_2d = (seed_indices + jitter) * px_size
    
    # 2. Add real-world coordinate offsets back in
    seed_coords_2d[:, 0] += x_min
    seed_coords_2d[:, 1] += y_min
    
    # 3. Assign realistic Z-coordinates to the seeds by sampling the real data
    z_samples = np.random.choice(real_z_coords, size=len(seed_coords_2d))
    
    # 4. Combine [X, Y] with [Z] to make an (N, 3) array
    seed_coords_3d = np.column_stack((seed_coords_2d, z_samples))
    
    return seed_coords_3d

def plant_seeds_POISSON(prob_map, real_z_coords, x_min, y_min, px_size=50.0, scaling_factor=1.0):
    """
    Plants seeds based on the 3D probability map using a Poisson distribution.
    Allows high-density pixels to spawn multiple seeds.
    """
    # 1. Calculate the expected number of seeds per pixel
    expected_seeds = prob_map * scaling_factor
    
    # 2. Draw actual seed counts from a Poisson distribution
    seed_counts = np.random.poisson(expected_seeds)
    
    # 3. Extract the pixel coordinates where seeds were generated
    # If a pixel got 3 seeds, we need to duplicate its coordinates 3 times
    unique_pixel_indices = np.argwhere(seed_counts > 0)
    counts_per_pixel = seed_counts[seed_counts > 0]
    
    seed_indices = np.repeat(unique_pixel_indices, counts_per_pixel, axis=0)
    
    # 4. Add random sub-pixel jitter (0 to 1) to spread them across the 50x50nm bin
    jitter = np.random.rand(*seed_indices.shape)
    seed_coords_2d = (seed_indices + jitter) * px_size
    
    # 5. Add real-world coordinate offsets back in
    seed_coords_2d[:, 0] += x_min
    seed_coords_2d[:, 1] += y_min
    
    # 6. Assign realistic Z-coordinates to the seeds by sampling the real data
    z_samples = np.random.choice(real_z_coords, size=len(seed_coords_2d))
    
    # 7. Combine [X, Y] with [Z] to make an (N, 3) array
    seed_coords_3d = np.column_stack((seed_coords_2d, z_samples))
    
    return seed_coords_3d

def plot_seeds_napari(me3_seeds=None, ac_seeds=None, seed_size=15):
    """
    Visualizes raw (N, 3) numpy matrices of ground-truth seeds in 3D.
    
    Args:
        me3_seeds (np.ndarray, optional): (N, 3) matrix of H3K27me3 seeds. Plotted in green.
        ac_seeds (np.ndarray, optional): (N, 3) matrix of H3K27ac seeds. Plotted in red.
        seed_size (int): Visual size of the points.
    """
    viewer = napari.Viewer()
    
    # Exact colors from your clustering script
    me3_color = [0, 0.7, 0, 1]  # Green
    ac_color = [0.7, 0, 0, 1]   # Red
    
    # Plot H3K27me3 Seeds
    if me3_seeds is not None and len(me3_seeds) > 0:
        viewer.add_points(
            data=me3_seeds,
            face_color=me3_color,
            border_width=0, 
            size=seed_size,    
            name="H3K27me3 Seeds",
            out_of_slice_display=True
        )
        
    # Plot H3K27ac Seeds
    if ac_seeds is not None and len(ac_seeds) > 0:
        viewer.add_points(
            data=ac_seeds,
            face_color=ac_color,
            border_width=0, 
            size=seed_size,    
            name="H3K27ac Seeds",
            out_of_slice_display=True
        )

    # Maintain the exact spatial orientation from your cluster viewer
    viewer.dims.order = (2, 1, 0)
    viewer.dims.ndisplay = 3
    
def extract_empirical_parameters_NEW(real_df, raw_me3_df, raw_ac_df, contour_bands, bin_size=50, sdis=200, step=10):
    """
    Extracts the 1D Radial Distribution Function (RDF) from experimental data.
    This captures the spatial architecture (tight vs. sprawling) of the domains.
    """
    comb_df = pd.concat([raw_me3_df, raw_ac_df])
    x_min = comb_df["x [nm]"].min()
    y_min = comb_df["y [nm]"].min()
    
    coords_2d = real_df[["x [nm]", "y [nm]"]].to_numpy()
    
    # --- 1. Extract Radial Density Profile (for Seed Placement) ---
    x_idx = np.clip(((coords_2d[:, 0] - x_min) // bin_size).astype(int), 0, contour_bands.shape[0] - 1)
    y_idx = np.clip(((coords_2d[:, 1] - y_min) // bin_size).astype(int), 0, contour_bands.shape[1] - 1)
    loc_bands = contour_bands[x_idx, y_idx]
    
    num_bands = contour_bands.max()
    radial_density_profile = np.zeros(num_bands)
    for i in range(1, num_bands + 1):
        band_area = np.sum(contour_bands == i)
        locs_in_band = np.sum(loc_bands == i)
        if band_area > 0:
            radial_density_profile[i-1] = locs_in_band / band_area
            
    # --- 2. Extract Radial Distribution Function (RDF) ---
    tree = cKDTree(coords_2d)
    radii = np.arange(step, sdis + step, step)
    
    # Calculate average density of the whole nucleus to normalize the RDF
    # This prevents 'scaling' issues between different cells
    nuc_area_nm2 = np.sum(contour_bands > 0) * (bin_size**2)
    avg_density = len(coords_2d) / nuc_area_nm2
    
    rdf_profile = []
    prev_counts = np.zeros(len(coords_2d))
    
    for r in radii:
        current_counts = tree.query_ball_point(coords_2d, r, return_length=True)
        ring_counts = current_counts - prev_counts
        prev_counts = current_counts
        
        # Mean neighbors in this ring across all points
        mean_neighbors = np.mean(ring_counts)
        
        # Normalize by the area of the ring to get local density
        ring_area = np.pi * (r**2 - (r-step)**2)
        local_density = mean_neighbors / ring_area
        
        # RDF = Local Density / Global Average Density
        # Values > 1 indicate clustering/enrichment at this radius
        rdf_profile.append(local_density / avg_density)
            
    z_degradation = np.full(len(radii), 1.5)
            
    return radial_density_profile, np.array(rdf_profile), z_degradation, x_min, y_min
    
def spawn_nanodomains_NEW(seed_coords, rdf_profile, z_degradation, step=10, base_intensity=0.5):
    """
    Spawns localizations using the RDF profile. 
    Higher base_intensity = denser nanodomains.
    """
    all_points = []
    radii = np.arange(step, (len(rdf_profile) * step) + step, step)
    
    for cluster_id, seed in enumerate(seed_coords):
        # Add the seed itself
        all_points.append([seed[0], seed[1], seed[2], cluster_id])
        
        for i, r in enumerate(radii):
            # The RDF value determines how many points to spawn at this distance
            # We scale it so that 'enrichment' creates multiple points
            num_to_spawn = int(np.round(rdf_profile[i] * base_intensity))
            
            if num_to_spawn <= 0:
                continue
                
            for _ in range(num_to_spawn):
                # Randomize radius within the 10nm bin
                radius = r - (np.random.rand() * step)
                
                # Standard 3D projection
                theta = np.random.rand() * 2 * np.pi 
                v = np.random.rand()
                phi = np.arccos(2 * v - 1) 
                
                dx = radius * np.sin(phi) * np.cos(theta)
                dy = radius * np.sin(phi) * np.sin(theta)
                dz = radius * np.cos(phi) * z_degradation[i]
                
                all_points.append([seed[0] + dx, seed[1] + dy, seed[2] + dz, cluster_id])
                
    df = pd.DataFrame(all_points, columns=["x [nm]", "y [nm]", "z [nm]", "cluster_label"])
    
    # Assign reproducible colors for visualization
    rng = np.random.default_rng(0)
    unique_clusters = df["cluster_label"].unique()
    cluster_colors = {c: np.append(rng.random(3), 1) for c in unique_clusters}
    df["cluster_color"] = [cluster_colors[c] for c in df["cluster_label"]]
    df['sigmax [nm]'] = 110.0 
    
    return df
    
def spawn_nanodomains_OLD(seed_coords, empirical_hists, z_degradation, step=10):
    """
    Uses Inverse Transform Sampling to spawn secondary localizations around seeds,
    matching the experimental radial distribution functions. 
    Tags all localizations with a cluster_label and assigns a random RGBA color.
    """
    all_points = []
    
    # Dynamically calculate the max distance based on the provided histograms
    # If empirical_hists has 15 bins, sdis becomes 150. If 20 bins, sdis becomes 200.
    sdis = len(empirical_hists) * step
    
    # Use enumerate to assign a unique cluster_id (0, 1, 2...) to each seed
    for cluster_id, seed in enumerate(seed_coords):
        
        # Add the seed itself to the list
        all_points.append([seed[0], seed[1], seed[2], cluster_id])
        
        # Step outward in increments
        for r_step in range(step, sdis + step, step):
            bin_idx = (r_step // step) - 1
            
            hist = empirical_hists[bin_idx]
            if np.sum(hist) == 0:
                continue
                
            pdf = hist / np.sum(hist)
            cdf = np.cumsum(pdf)
            
            roll = np.random.rand()
            num_to_spawn = np.searchsorted(cdf, roll) 
            
            if num_to_spawn == 0:
                continue
                
            # Project the spawned points in 3D space
            for _ in range(num_to_spawn):
                radius = r_step + (np.random.rand() * step)
                
                theta = np.random.rand() * 2 * np.pi 
                v = np.random.rand()
                phi = np.arccos(2 * v - 1) 
                
                dx = radius * np.sin(phi) * np.cos(theta)
                dy = radius * np.sin(phi) * np.sin(theta)
                
                dz = radius * np.cos(phi) * z_degradation[bin_idx]
                
                all_points.append([seed[0] + dx, seed[1] + dy, seed[2] + dz, cluster_id])
                
    # Create the base DataFrame
    df = pd.DataFrame(all_points, columns=["x [nm]", "y [nm]", "z [nm]", "cluster_label"])
    
    # ==========================================
    # --- COLOR ASSIGNMENT ---
    # ==========================================
    rng = np.random.default_rng(0)  # Reproducible colors
    unique_clusters = df["cluster_label"].unique()
    
    # Create a dictionary mapping each cluster ID to a random color array
    cluster_colors = {c: np.append(rng.random(3), 1) for c in unique_clusters}
    
    df["cluster_color"] = [cluster_colors[c] for c in df["cluster_label"]]
    df['sigmax [nm]'] = 110.0 
    
    return df
 
def extract_empirical_parameters_OLD(real_df, raw_me3_df, raw_ac_df, contour_bands, bin_size=50, sdis=200, step=10):
    comb_df = pd.concat([raw_me3_df, raw_ac_df])
    x_min = comb_df["x [nm]"].min()
    y_min = comb_df["y [nm]"].min()
    
    coords_2d = real_df[["x [nm]", "y [nm]"]].to_numpy()
    
    # --- 1. Extract Radial Density Profile ---
    x_idx = np.clip(((coords_2d[:, 0] - x_min) // bin_size).astype(int), 0, contour_bands.shape[0] - 1)
    y_idx = np.clip(((coords_2d[:, 1] - y_min) // bin_size).astype(int), 0, contour_bands.shape[1] - 1)
    
    loc_bands = contour_bands[x_idx, y_idx]
    
    num_bands = contour_bands.max()
    radial_density_profile = np.zeros(num_bands)
    for i in range(1, num_bands + 1):
        band_area = np.sum(contour_bands == i)
        locs_in_band = np.sum(loc_bands == i)
        if band_area > 0:
            radial_density_profile[i-1] = locs_in_band / band_area
            
    # --- 2. Extract Empirical Histograms (KD-Tree) ---
    tree = cKDTree(coords_2d)
    radii = np.arange(step, sdis + step, step)
    
    empirical_hists = []
    prev_counts = np.zeros(len(coords_2d))
    
    for r in radii:
        # Get the cumulative number of neighbors for EACH individual point
        current_counts = tree.query_ball_point(coords_2d, r, return_length=True)
        
        # Isolate the neighbors strictly within this specific 10nm ring
        ring_counts_per_point = current_counts - prev_counts
        prev_counts = current_counts
        
        # Build the histogram: How many points have 0 neighbors? 1 neighbor? 2?
        max_neighbors = int(ring_counts_per_point.max())
        if max_neighbors == 0:
            empirical_hists.append(np.array([1.0])) # 100% chance of 0 neighbors
        else:
            hist, _ = np.histogram(ring_counts_per_point, bins=np.arange(max_neighbors + 2))
            empirical_hists.append(hist)
            
    z_degradation = np.full(len(radii), 1.5)
            
    # Return x_min and y_min so we can align the generated seeds to real physical space
    return radial_density_profile, empirical_hists, z_degradation, x_min, y_min

def epsilon_cost(eps, coords, gt_target, min_samples=8, size_metric='geometric_2D'):
    """
    Evaluates how close the average DBSCAN cluster size is to the ground truth target.
    Uses robust geometric spanning (90th percentile) to handle sparse, 
    jagged data caused by low labeling efficiencies.
    """
    # 1. Run DBSCAN on the 2D plane
    labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(coords[:, :2])
    valid_mask = labels != -1
    
    # Severe penalty if DBSCAN fails to find ANY valid clusters
    if not np.any(valid_mask):
        return 99999.0 
        
    df = pd.DataFrame(coords[valid_mask], columns=["x [nm]", "y [nm]", "z [nm]"])
    df['label'] = labels[valid_mask]
    
    sizes = []
    
    for cluster_id, group_df in df.groupby('label'):
        cluster_coords = group_df[['x [nm]', 'y [nm]']].values
        
        # PENALTY 1: Too few points to form a biologically meaningful shape
        if len(cluster_coords) < 10:
            sizes.append(20.0) 
            continue
            
        # ====================================================
        # ROBUST METHOD: 90th Percentile Geometric Span
        # ====================================================
        if size_metric == 'geometric_2D':
            # 1. Find the center of mass
            centroid = np.mean(cluster_coords, axis=0)
            
            # 2. Calculate the distance of every point to the center
            distances = np.linalg.norm(cluster_coords - centroid, axis=1)
            
            # 3. Use the 90th percentile distance to define the boundary
            radius = np.percentile(distances, 75)
            diameter = 2 * radius
            
            # PENALTY 2: Continent Check
            # If the domain is larger than 800nm, it's definitively merged clusters
            if diameter > 800.0:
                sizes.append(1000.0)
            else:
                sizes.append(diameter)
                
        # ====================================================
        # FALLBACK / LEGACY LOGIC
        # ====================================================
        else:
            try:
                results = calc_nanodomain_size(group_df)
                sizes.append(results.get(size_metric, results.get('hull_3D', 0.0)))
            except NameError:
                sizes.append(0.0)
            
    # Penalty if all identified clusters were completely invalid/penalized
    if not sizes:
        return 99999.0
        
    # Calculate final absolute error cost
    mean_size = np.mean(sizes)
    return abs(mean_size - gt_target)

def optimize_epsilon_hybrid(fraction_dfs, gt_targets, calibration_area_nm2, min_samples=8, size_metric='geometric_2D', show_plot=True):
    """
    Finds the optimal DBSCAN epsilon using a Coarse Grid Search followed by 
    a bounded scalar minimization for high-precision refinement.
    """
    densities = []
    best_epsilons = []
    
    print(f"Starting HYBRID (Grid + Minimize) search over {len(fraction_dfs)} fractions...")
    
    for i, df in enumerate(fraction_dfs):
        coords = df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
        if len(coords) < min_samples:
            continue
            
        density = len(coords) / calibration_area_nm2
        current_target = gt_targets[i]
        
        # ==========================================
        # STAGE 1: COARSE GRID SEARCH (15 nm steps)
        # ==========================================
        # We use wide steps to quickly map the landscape
        coarse_grid = np.arange(10, 225, 15.0)
        best_coarse_eps = 10.0
        best_coarse_cost = 99999.0
        
        for test_eps in coarse_grid:
            cost = epsilon_cost(test_eps, coords, current_target, min_samples, size_metric)
            if cost < best_coarse_cost:
                best_coarse_cost = cost
                best_coarse_eps = test_eps
            
            # Stop if we hit a massive continent penalty after finding a good candidate
            if cost > 500.0 and best_coarse_cost < 100.0:
                break

        # ==========================================
        # STAGE 2: FINE MINIMIZATION
        # ==========================================
        # Define a 30nm search window around our coarse winner
        search_min = max(10.0, best_coarse_eps - 15.0)
        search_max = min(225.0, best_coarse_eps + 15.0)
        
        res = minimize_scalar(
            epsilon_cost,
            bounds=(search_min, search_max),
            args=(coords, current_target, min_samples, size_metric),
            method='bounded',
            options={'xatol': 0.2} # High precision refinement
        )
        
        # --- ACCEPTANCE ---
        if res.success and res.fun < 50.0:
            densities.append(density)
            best_epsilons.append(res.x)
            print(f"Fraction {i+1}: Density {density:.6e} | Target {current_target:.1f}nm -> Best Eps {res.x:.2f}nm (Error: {res.fun:.2f}nm)")
        else:
            print(f"Fraction {i+1}: Density {density:.6e} | FAILED (Error: {res.fun if res.success else 'Opt Failed'})")

    if len(densities) < 2:
        print("\nError: Not enough valid data points for power law.")
        return 0, 0, densities, best_epsilons
            
    # --- FIT THE POWER LAW ---
    log_dens = np.log10(densities)
    log_eps = np.log10(best_epsilons)
    par = np.polyfit(log_dens, log_eps, 1)
    
    exp = par[0]
    coef = 10**(par[1])
    r_squared = np.corrcoef(log_dens, log_eps)[0, 1]**2
    
    print(f"\nHybrid Calibration Complete! Eps = {coef:.2f} * Density^{exp:.2f} (R^2 = {r_squared:.4f})")
    
    if show_plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
            
        # Subplot 1: Log-Log Scale
        axes[0].scatter(log_dens, log_eps, color='blue', label='Optimized Epsilons')
        fit_line = par[0] * log_dens + par[1]
        axes[0].plot(log_dens, fit_line, color='red', linestyle='--', label=f'Fit: eps = {coef:.2f} * dens^{exp:.2f}')
        axes[0].set_xlabel('Log10(Density [locs/nm$^2$])')
        axes[0].set_ylabel('Log10(Epsilon [nm])')
        axes[0].set_title('Log-Log Scale')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Subplot 2: Standard Scale
        axes[1].scatter(densities, best_epsilons, color='blue', label='Optimized Epsilons')
        
        # Generate smooth curve for standard scale plot
        dens_vals = np.linspace(min(densities), max(densities), 100)
        eps_vals = coef * (dens_vals ** exp)
        axes[1].plot(dens_vals, eps_vals, color='red', linestyle='--', label=f'Fit: eps = {coef:.2f} * dens^{exp:.2f}')
        axes[1].set_xlabel('Density [locs/nm$^2$]')
        axes[1].set_ylabel('Epsilon [nm]')
        axes[1].set_title('Standard Scale')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.suptitle(f'DBSCAN Epsilon Density Calibration\n$R^2$ = {r_squared:.4f}')
        plt.tight_layout()
        plt.show()
        
    return coef, exp, densities, best_epsilons

def epsilon_cost_OLD(eps, coords, gt_target, min_samples=8, size_metric='radial_1D', step=10, max_radius=200):
    """
    Evaluates how close the average DBSCAN cluster size is to the ground truth target.
    Handles sparse/fragmented data safely by actively penalizing tiny fragments 
    and massive merged continents.
    """
    # 1. Run DBSCAN on the 2D plane
    labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(coords[:, :2])
    valid_mask = labels != -1
    
    # Severe penalty if DBSCAN fails to find ANY valid clusters
    if not np.any(valid_mask):
        return 99999.0 
        
    df = pd.DataFrame(coords[valid_mask], columns=["x [nm]", "y [nm]", "z [nm]"])
    df['label'] = labels[valid_mask]
    
    sizes = []
    
    for cluster_id, group_df in df.groupby('label'):
        
        if size_metric == 'radial_1D':
            cluster_coords = group_df[['x [nm]', 'y [nm]']].values
            
            # PENALTY 1 (Noise Fragment): Too few points to be a real domain
            if len(cluster_coords) < 10:
                sizes.append(20.0) 
                continue
                
            tree = cKDTree(cluster_coords)
            radii = np.arange(step, max_radius + step, step)
            
            prev_counts = np.zeros(len(cluster_coords))
            ring_counts = []
            
            for r in radii:
                current_counts = tree.query_ball_point(cluster_coords, r, return_length=True)
                ring_counts.append(current_counts - prev_counts)
                prev_counts = current_counts
                
            between_matrix = np.column_stack(ring_counts)
            mcount = np.mean(between_matrix, axis=0)
            
            ring_areas = np.pi * (radii**2) - np.pi * ((radii - step)**2)
            dens = mcount / ring_areas
            dens = np.nan_to_num(dens, nan=0.0, posinf=0.0, neginf=0.0)
            
            raw_max = np.max(dens)
            raw_min = np.min(dens)
            
            # PENALTY 1B (Flatline Noise)
            if raw_max <= 0:
                sizes.append(20.0)
                continue
                
            # PENALTY 2 (Continent Check): Outer edge density is still > 50% of the peak.
            # This is a massive merged cluster. Assign a massive size.
            if raw_min > (0.5 * raw_max):
                sizes.append(1000.0)
                continue
                
            dens_norm = dens - raw_min
            max_density = np.max(dens_norm)
            
            # PENALTY 1C (Normalized Flatline)
            if max_density <= 0:
                sizes.append(20.0)
                continue
                
            # Finally, find the 20% drop-off threshold
            threshold = 0.2 * max_density
            drop_indices = np.where(dens_norm < threshold)[0]
            
            if len(drop_indices) > 0:
                width = 2 * ((drop_indices[0] + 1) * step)
                sizes.append(width)
            else:
                # SAFEGUARD: Valid cluster that decays very slowly. Cap at window size.
                sizes.append(2 * max_radius)

        else:
            # Fallback for volumetric sizing (hull_3D, bb_2D, etc.)
            # Assumes calc_nanodomain_size is available in your namespace
            try:
                results = calc_nanodomain_size(group_df)
                sizes.append(results.get(size_metric, results.get('hull_3D', 0.0)))
            except NameError:
                print("Warning: calc_nanodomain_size not defined.")
                sizes.append(0.0)
            
    # Penalty if all identified clusters were completely invalid
    if not sizes:
        return 99999.0
        
    # Calculate final absolute error cost
    mean_size = np.mean(sizes)
    return abs(mean_size - gt_target)

def optimize_epsilon_grid_OLD(fraction_dfs, gt_targets, calibration_area_nm2, min_samples=8, size_metric='radial_1D', show_plot=True, debug=False):
    """
    Finds the optimal DBSCAN epsilon using a robust Grid Search.
    Bypasses local-minima traps, handles flat plateaus by choosing the tightest epsilon, 
    and fits a log-log power law curve to the valid results.
    """
    densities = []
    best_epsilons = []
    
    # Define our Grid: Test every 2 nm from 10 nm to 220 nm
    eps_grid = np.arange(20, 71, 0.5)
    
    print(f"Starting grid search over {len(fraction_dfs)} density fractions...")
    
    for i, df in enumerate(fraction_dfs):
        coords = df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
        
        # Skip if there aren't even enough points to form a single cluster
        if len(coords) < min_samples:
            print(f"Fraction {i+1}: Skipped (Too few total localizations)")
            continue
            
        # Calculate spatial density
        density = len(coords) / calibration_area_nm2
        current_target = gt_targets[i]
        
        best_eps = None
        best_cost = 99999.0
        
        # --- THE GRID SEARCH ---
        for test_eps in eps_grid:
            cost = epsilon_cost(test_eps, coords, current_target, min_samples, size_metric)
            if debug: print(f"-> Test Eps: {test_eps} | Cost: {cost}")
            # Use strict less-than (<) to handle flat plateaus. 
            # This ensures we keep the FIRST (smallest) epsilon that hits the lowest cost.
            if cost < best_cost:
                best_cost = cost
                best_eps = test_eps
                
        # --- ACCEPTANCE CRITERIA ---
        # 1. We found an epsilon
        # 2. The error was within 50 nm of our target
        # 3. We didn't just hit the artificial ceiling of our grid
        if best_eps is not None and best_cost < 50.0 and best_eps < eps_grid[-1]:
            densities.append(density)
            best_epsilons.append(best_eps)
            print(f"Fraction {i+1}: Density {density:.6e} | Target {current_target:.1f}nm -> Best Eps {best_eps:.1f}nm (Error: {best_cost:.1f}nm)")
        else:
            print(f"Fraction {i+1}: Density {density:.6e} | FAILED (Best Error: {best_cost:.1f}nm at Eps {best_eps})")

    # Verify we have enough points to draw a regression line
    if len(densities) < 2:
        print("\nError: Not enough valid density states to calculate a power law.")
        return 0, 0, densities, best_epsilons
            
    # --- FIT THE POWER LAW ---
    log_dens = np.log10(densities)
    log_eps = np.log10(best_epsilons)
    par = np.polyfit(log_dens, log_eps, 1)
    
    exp = par[0]
    coef = 10**(par[1])
    
    correlation_matrix = np.corrcoef(log_dens, log_eps)
    r_squared = correlation_matrix[0, 1]**2
    print(f"\nCalibration Complete! Formula: Eps = {coef:.2f} * Density^{exp:.2f} (R^2 = {r_squared:.4f})")
    
    # --- PLOT THE RESULTS ---
    if show_plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        
        # Subplot 1: Log-Log Scale
        axes[0].scatter(log_dens, log_eps, color='blue', label='Optimized Epsilons')
        fit_line = par[0] * log_dens + par[1]
        axes[0].plot(log_dens, fit_line, color='red', linestyle='--', label=f'Fit: eps = {coef:.2f} * dens^{exp:.2f}')
        axes[0].set_xlabel('Log10(Density [locs/nm$^2$])')
        axes[0].set_ylabel('Log10(Epsilon [nm])')
        axes[0].set_title('Log-Log Scale')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Subplot 2: Standard Scale
        axes[1].scatter(densities, best_epsilons, color='blue', label='Optimized Epsilons')
        
        # Generate smooth curve for standard scale plot
        dens_vals = np.linspace(min(densities), max(densities), 100)
        eps_vals = coef * (dens_vals ** exp)
        axes[1].plot(dens_vals, eps_vals, color='red', linestyle='--', label=f'Fit: eps = {coef:.2f} * dens^{exp:.2f}')
        axes[1].set_xlabel('Density [locs/nm$^2$]')
        axes[1].set_ylabel('Epsilon [nm]')
        axes[1].set_title('Standard Scale')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        
        plt.suptitle(f'DBSCAN Epsilon Density Calibration\n$R^2$ = {r_squared:.4f}')
        plt.tight_layout()
        plt.show()
        
    return coef, exp, densities, best_epsilons

def calculate_simulated_gt_radial(simulated_df, step=10, max_radius=200):
    """
    Calculates the average 1D radial diameter of simulated nanodomains.
    This applies the radial density drop-off method (20% threshold) to 
    isolated clusters, avoiding cross-talk.
    
    Args:
        simulated_df (pd.DataFrame): DataFrame containing 'x [nm]', 'y [nm]', 
                                     and 'cluster_label'.
        step (int): Ring increment in nm (Default: 10).
        max_radius (int): Maximum search radius in nm (Default: 200).
        
    Returns:
        float: The average 1D ground-truth diameter of the clusters (in nm).
    """
    cluster_widths = []
    
    # Iterate through each isolated ground-truth cluster
    for cluster_id, group in simulated_df.groupby('cluster_label'):
        
        # Skip noise if it exists in the dataframe
        if cluster_id == -1:
            continue
            
        coords = group[['x [nm]', 'y [nm]']].values
        
        # Require a minimum number of points to build a valid density profile
        if len(coords) < 10:
            continue
            
        # 1. Build the KD-Tree for fast neighbor lookup
        tree = cKDTree(coords)
        radii = np.arange(step, max_radius + step, step)
        
        prev_counts = np.zeros(len(coords))
        ring_counts_list = []
        
        # 2. Count neighbors in expanding rings (Pair-Correlation)
        for r in radii:
            # Query the total neighbors within radius 'r'
            current_counts = tree.query_ball_point(coords, r, return_length=True)
            
            # Subtract the inner circle to get the count just for this specific ring
            ring_counts = current_counts - prev_counts
            ring_counts_list.append(ring_counts)
            
            prev_counts = current_counts
            
        # Shape: (N_points, N_rings)
        between_matrix = np.column_stack(ring_counts_list)
        
        # 3. Calculate Spatial Density
        # Average neighbors in each ring across the whole cluster
        mcount = np.mean(between_matrix, axis=0)
        
        # Calculate the mathematical area of each 2D ring
        ring_areas = np.pi * (radii**2) - np.pi * ((radii - step)**2)
        
        # Density = points / nm^2
        dens = mcount / ring_areas
        
        # Normalize the density baseline (Matches MATLAB logic)
        dens = dens - np.min(dens)
        
        # 4. Find the 20% Drop-off Threshold
        threshold = 0.2 * np.max(dens)
        
        # Find the first ring index where density falls below the threshold
        drop_indices = np.where(dens < threshold)[0]
        
        if len(drop_indices) > 0:
            # Convert the index back to a physical radius
            # Adding 1 accounts for 0-based Python indexing (Index 0 is the 10nm ring)
            radius = (drop_indices[0] + 1) * step
            
            # Calculate 1D diameter
            width = 2 * radius
            cluster_widths.append(width)
            
    # Return the mean diameter across all valid measured clusters
    if len(cluster_widths) > 0:
        return np.mean(cluster_widths)
    else:
        return 0.0

def calculate_geometric_gt_target(master_df):
    """
    Calculates the 90th-percentile geometric diameter of the 
    simulated domains using their ground-truth labels.
    """
    gt_sizes = []
    
    # Use the 'cluster_label' assigned during the spawn_nanodomains step
    for _, group in master_df.groupby('cluster_label'):
        # Ignore background noise (label -1)
        if _ == -1 or len(group) < 10:
            continue
            
        coords = group[['x [nm]', 'y [nm]']].values
        centroid = np.mean(coords, axis=0)
        distances = np.linalg.norm(coords - centroid, axis=1)
        
        # 90th percentile radius -> Diameter
        radius = np.percentile(distances, 75)
        gt_sizes.append(2 * radius)
        
    return np.mean(gt_sizes)

"""
TODO: 
SNCR Function (KD Tree most likely)
Density Correction
Ripley's K implementation
Final outputs
"""