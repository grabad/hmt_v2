"""
DBSCAN clustering of localizations into nanodomains, and per-nanodomain
size / density measurements (manuscript Fig. 2d/e, Fig. 6g/h).

Use cluster_dbscan_density_corrected with a calibration from
calibrate.optimize_dbscan_size_matching so eps / min_samples follow each cell's
own labelling density; cluster_dbscan is the fixed-parameter primitive.
Nanodomain sizes from calculate_nanodomain_characteristics use the same
geometric_diameter_2D statistic the calibration targets.
"""
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, QhullError
from sklearn.cluster import DBSCAN


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
    of local labelling density. See calibrate.optimize_dbscan_size_matching, which
    fits eps = eps_coef * density^eps_exp and min_samples = min_samples_slope *
    density + min_samples_intercept from simulated ground truth.

    Computes loc_df's own localization density (len(loc_df) / area_nm2) and
    evaluates those fits at that density, so the clustering parameters used here
    scale with how densely THIS dataset happens to be labelled — unlike a fixed
    ("arbitrary") eps/min_samples pair, which over- or under-merges nanodomains
    as density varies and biases any downstream size measurement.

    Args:
        loc_df:           DataFrame with [x [nm], y [nm]] (and [z [nm]] if use_z=True)
        calibration:       dict returned by calibrate.optimize_dbscan_size_matching;
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


def calculate_nanodomain_characteristics(locs, percentile=75.0):
    """
    Calculates characteristics of nanodomains found via DBSCAN clustering, one row per cluster.

    Size measures:
        bb_size_2D:            mean of the xy bounding-box width and height (nm) -- the
                               manuscript's "nanodomain size"
        bb_size_3D:            mean of the xyz bounding-box extents (nm)
        geometric_diameter_2D: 2x the `percentile`-th percentile of xy point-to-centroid
                               distance (nm) -- the same statistic simulate.
                               calculate_geometric_gt_target / domain_size_cost calibrate
                               against, so it is directly comparable to a calibration's gt_target
        hull_size_2D / _3D:    convex hull area (nm^2) / volume (nm^3)
        density_2D / _3D:      localizations per hull area / volume

    Noise (cluster_label == -1) and clusters with fewer than 4 localizations are skipped.
    3D measures are NaN when the 3D hull is degenerate (e.g. every point at one z, which
    happens with integer-binned z); the 2D measures are still reported for that cluster.

    Args:
        locs:       DataFrame with [x [nm], y [nm], z [nm], cluster_label]
                    (and optionally cluster_color), e.g. from cluster_dbscan
        percentile: percentile for geometric_diameter_2D (match the calibration's)

    Returns:
        Pandas DataFrame: cluster_df
    """
    column_names = [
        "label", "n_locs", "cluster_color",
        "centroid", "centroid_x", "centroid_y", "centroid_z", "vertices",
        "hull_size_2D", "hull_size_3D",
        "bb_size_2D", "bb_size_3D", "geometric_diameter_2D",
        "density_2D", "density_3D",
    ]
    has_color = "cluster_color" in locs.columns
    cluster_rows = []

    clustered = locs[locs["cluster_label"] != -1]
    for label, cluster_data in clustered.groupby("cluster_label"):
        if len(cluster_data) < 4:
            continue

        points_3D = cluster_data[["x [nm]", "y [nm]", "z [nm]"]].to_numpy(dtype=float)
        points_2D = points_3D[:, :2]

        try:
            hull_2D = ConvexHull(points_2D)
        except QhullError:
            continue  # collinear in xy -- no meaningful size

        try:
            hull_3D = ConvexHull(points_3D)
            hull_size_3D = hull_3D.volume
            vertices = points_3D[hull_3D.vertices]
        except QhullError:
            hull_size_3D = np.nan
            vertices = points_3D[hull_2D.vertices]

        extents = points_3D.max(axis=0) - points_3D.min(axis=0)
        centroid = points_3D.mean(axis=0)
        radial_2D = np.linalg.norm(points_2D - centroid[:2], axis=1)
        n = len(points_3D)

        cluster_rows.append({
            "label": label,
            "n_locs": n,
            "cluster_color": cluster_data["cluster_color"].iloc[0] if has_color else None,
            "centroid": centroid,
            "centroid_x": centroid[0],
            "centroid_y": centroid[1],
            "centroid_z": centroid[2],
            "vertices": vertices,
            "hull_size_2D": hull_2D.volume,
            "hull_size_3D": hull_size_3D,
            "bb_size_2D": extents[:2].mean(),
            "bb_size_3D": extents.mean(),
            "geometric_diameter_2D": 2.0 * np.percentile(radial_2D, percentile),
            "density_2D": n / hull_2D.volume if hull_2D.volume > 0 else np.nan,
            "density_3D": n / hull_size_3D if hull_size_3D > 0 else np.nan,
        })

    return pd.DataFrame(cluster_rows, columns=column_names)
