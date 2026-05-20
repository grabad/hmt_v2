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