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

    return sampled_data 

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

def spawn_nanodomains(seed_coords, empirical_hists, z_degradation, sdis=200, step=10):
    """
    Uses Inverse Transform Sampling to spawn secondary localizations around seeds,
    matching the experimental radial distribution functions. 
    Tags all localizations with a cluster_label and assigns a random RGBA color.
    """
    all_points = []
    
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

def extract_empirical_parameters(real_df, raw_me3_df, raw_ac_df, contour_bands, bin_size=50, sdis=200, step=10):
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

def epsilon_cost(eps, coords, gt_target, min_samples=8, size_metric='hull_3D'):
    labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(coords[:, :2])
    valid_mask = labels != -1
    
    if not np.any(valid_mask):
        return 99999.0 
        
    df = pd.DataFrame(coords[valid_mask], columns=["x [nm]", "y [nm]", "z [nm]"])
    df['label'] = labels[valid_mask]
    
    sizes = []
    for _, group_df in df.groupby('label'):
        results = calc_nanodomain_size(group_df)
        
        # FIX 1: Match the uppercase 'D' in 'hull_3D'
        # Also using size_metric directly as the fallback key
        sizes.append(results.get(size_metric, results['hull_3D']))
        
    mean_size = np.mean(sizes)
    return abs(mean_size - gt_target)


def optimize_epsilon_power_law(fraction_dfs, gt_targets, nuclear_area_nm2, min_samples=8, size_metric='bb_2D', show_plot=True):
    """
    Runs the epsilon search across all density fractions and fits the log-log power law.
    Returns the coefficients needed to dynamically cluster new cells.
    """
    densities = []
    best_epsilons = []
    
    for i, df in enumerate(fraction_dfs):
        coords = df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
        if len(coords) < min_samples:
            continue
            
        density = len(coords) / nuclear_area_nm2
        
        # FIX 2: Pass only the SPECIFIC target for this fraction (gt_targets[i])
        # instead of the whole list
        current_target = gt_targets[i]
        
        res = minimize_scalar(
            epsilon_cost,
            bounds=(10, 200),
            args=(coords, current_target, min_samples, size_metric),
            method='bounded'
        )
        
        if res.success:
            densities.append(density)
            best_epsilons.append(res.x)
            print(f"Fraction {i+1}: Density {density:.6f}\n", f"Target: {current_target:.2f} -> Ideal Eps: {res.x:.2f} nm", sep="   ")
            
    # Fit the linear regression to the log-log relationship
    log_dens = np.log10(densities)
    log_eps = np.log10(best_epsilons)
    par = np.polyfit(log_dens, log_eps, 1)
    
    exp = par[0]
    coef = 10**(par[1])
    
    correlation_matrix = np.corrcoef(log_dens, log_eps)
    r_squared = correlation_matrix[0, 1]**2
    
    if show_plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        
        # Subplot 1: Log-Log Scale
        axes[0].scatter(log_dens, log_eps, color='blue', label='Optimized Epsilons')
        fit_line = par[0] * log_dens + par[1]
        axes[0].plot(log_dens, fit_line, color='red', linestyle='--', label=f'Fit: eps = {coef:.2f} * dens^{exp:.2f}')
        axes[0].set_xlabel('Log10(Density [locs/nm^2])')
        axes[0].set_ylabel('Log10(Epsilon [nm])')
        axes[0].set_title('Log-Log Scale')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Subplot 2: Standard Scale
        axes[1].scatter(densities, best_epsilons, color='blue', label='Optimized Epsilons')
        dens_vals = np.linspace(min(densities), max(densities), 100)
        eps_vals = coef * (dens_vals ** exp)
        axes[1].plot(dens_vals, eps_vals, color='red', linestyle='--', label=f'Fit: eps = {coef:.2f} * dens^{exp:.2f}')
        axes[1].set_xlabel('Density [locs/nm^2]')
        axes[1].set_ylabel('Epsilon [nm]')
        axes[1].set_title('Standard Scale')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        
        plt.suptitle(f'DBSCAN Epsilon Density Calibration\n$R^2$ = {r_squared:.4f}')
        plt.tight_layout()
        plt.show()
        
    return coef, exp, densities, best_epsilons


"""
TODO: 
SNCR Function (KD Tree most likely)
Density Correction
Ripley's K implementation
Final outputs
"""