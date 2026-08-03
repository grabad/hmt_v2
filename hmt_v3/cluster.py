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

def cluster_dbscan(loc_df, eps=50, min_samples=8, use_z=False, n_jobs=-1):
    """
    Cluster localization data using scikit-learn's DBSCAN algorithm, taking in epsilon and min_samples parameters.
    Clusters will be assigned a numeric ID and a random RGB color, with unclustered noise labeled as -1 and set gray.

    Args:
        loc_df:       DataFrame with [x [nm], y [nm]] (and [z [nm]] if use_z=True)
        eps:          DBSCAN neighbourhood radius in nm
        min_samples:  DBSCAN min_samples parameter
        use_z:        if True, cluster in 3D (XYZ); if False (default), 2D (XY only) —
                      must match whatever dimensionality eps/min_samples were chosen for
        n_jobs:       forwarded to DBSCAN

    Returns:
        Pandas DataFrame: loc_df (with cluster_label and cluster_color columns)
    """
    loc_df = loc_df.copy()

    cols = ["x [nm]", "y [nm]", "z [nm]"] if use_z else ["x [nm]", "y [nm]"]
    loc_np = loc_df[cols].to_numpy()
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

def cluster_dbscan_density_corrected(loc_df, calibration, area_nm2, use_z=False,
                                     min_eps=10.0, min_min_samples=3, n_jobs=-1):
    """
    Cluster localization data with DBSCAN, using an (eps, min_samples) pair derived
    from a density calibration instead of a single fixed value applied irrespective
    of local labelling density. See simulate.optimize_dbscan_size_matching, which
    fits eps = eps_coef * density^eps_exp and min_samples = min_samples_slope *
    density + min_samples_intercept from simulated ground truth.

    Computes loc_df's own localization density (len(loc_df) / area_nm2) and
    evaluates those fits at that density, so the clustering parameters used here
    scale with how densely THIS dataset happens to be labelled — unlike a fixed
    ("arbitrary") eps/min_samples pair, which over- or under-merges nanodomains
    as density varies and biases any downstream size measurement.

    Args:
        loc_df:           DataFrame with [x [nm], y [nm]] (and [z [nm]] if use_z=True)
        calibration:       dict returned by simulate.optimize_dbscan_size_matching;
                          must have converged (at least 2 fitted densities)
        area_nm2:          imaged area in nm^2 used to compute loc_df's density —
                          must be measured the same way as the calibration's
                          calibration_area_nm2 (e.g. nucleus mask area, not a
                          bounding box) for the density values to be comparable
        use_z:             forwarded to cluster_dbscan; should match the use_z the
                          calibration was fit with
        min_eps:           floor on the corrected eps (nm), guards against
                          extrapolating the power-law fit to an unusably small value
        min_min_samples:   floor on the corrected min_samples
        n_jobs:            forwarded to DBSCAN

    Returns:
        tuple: (clustered_df, params) where params is a dict with the 'density'
               loc_df was measured at and the 'eps' / 'min_samples' DBSCAN used
    """
    if len(calibration.get('densities', [])) < 2:
        raise ValueError(
            "calibration did not converge (fewer than 2 fitted densities); "
            "its eps/min_samples fits are not usable.")

    density = len(loc_df) / area_nm2

    eps = calibration['eps_coef'] * density ** calibration['eps_exp']
    min_samples = calibration['min_samples_slope'] * density + calibration['min_samples_intercept']

    eps = max(eps, min_eps)
    min_samples = max(round(min_samples), min_min_samples)

    clustered_df = cluster_dbscan(loc_df, eps=eps, min_samples=min_samples, use_z=use_z, n_jobs=n_jobs)
    params = {'density': density, 'eps': eps, 'min_samples': min_samples}
    return clustered_df, params

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