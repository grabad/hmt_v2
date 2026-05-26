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

_default_bins = np.arange(100)
_default_rdf = _default_bins * np.exp(-0.5 * (_default_bins / 15.0) ** 2)
_default_adf = np.exp(-0.5 * (_default_bins / 15.0) ** 2)

def assign_contour_bands(df, contour_bands, bin_size=50):
    """
    Assigns each localization to a radial contour band from create_radial_contours.
    Adds a 'contour_band' column to the DataFrame (0 = outside nucleus).

    Args:
        df:             DataFrame with columns [x [nm], y [nm]]
        contour_bands:  2D int array from create_radial_contours
        bin_size:       pixel size in nm (must match what was used to build contour_bands)

    Returns:
        df with a new 'contour_band' column
    """
    df = df.copy()
    coords = df[["x [nm]", "y [nm]"]].to_numpy()

    x_min = coords[:, 0].min()
    y_min = coords[:, 1].min()

    x_idx = np.clip(((coords[:, 0] - x_min) // bin_size).astype(int), 0, contour_bands.shape[0] - 1)
    y_idx = np.clip(((coords[:, 1] - y_min) // bin_size).astype(int), 0, contour_bands.shape[1] - 1)

    df["contour_band"] = contour_bands[x_idx, y_idx]
    return df


def extract_empirical_parameters(df, sdis=200, step=10):
    """
    Extracts the RDF and ADF from experimental localization data,
    ready to pass directly into spawn_nanodomains.

    RDF: mean 2D radial neighbor density in concentric XY rings at each bin.
         Captures the lateral clustering geometry of the nanodomains.
    ADF: histogram of pairwise |z-offsets| between neighbors within sdis.
         Captures the axial spread independently of the radial distribution.

    Pass in localizations from a single channel (e.g. only H3K27me3) so the
    statistics reflect within-population clustering, not cross-channel mixing.

    Args:
        df:   DataFrame with columns [x [nm], y [nm], z [nm]]
        sdis: max search radius in nm
        step: bin width in nm

    Returns:
        rdf: 1D array of length sdis//step — radial density weights per bin
        adf: 1D array of length sdis//step — axial density weights per bin
    """
    coords = df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
    coords_2d = coords[:, :2]
    radii = np.arange(step, sdis + step, step)

    # --- RDF: mean 2D ring neighbor density ---
    tree_2d = cKDTree(coords_2d)
    rdf = np.zeros(len(radii))
    prev_counts = np.zeros(len(coords))

    for i, r in enumerate(radii):
        current_counts = np.array(
            tree_2d.query_ball_point(coords_2d, r, return_length=True), dtype=float
        ) - 1  # subtract self
        ring_counts = current_counts - prev_counts
        prev_counts = current_counts

        ring_area = np.pi * (r**2 - (r - step)**2)
        rdf[i] = np.mean(ring_counts) / ring_area

    # --- ADF: pairwise |z-offset| histogram ---
    tree_3d = cKDTree(coords)
    pairs = tree_3d.query_pairs(sdis, output_type='ndarray')
    dz = np.abs(coords[pairs[:, 0], 2] - coords[pairs[:, 1], 2])
    adf, _ = np.histogram(dz, bins=np.arange(0, sdis + step, step))

    # Subtract background (mean of last 10% of bins) and clip to zero
    n_tail = max(1, len(rdf) // 10)
    baseline = rdf[-n_tail:].mean()
    rdf = np.clip(rdf - baseline, 0, None)

    return rdf, adf.astype(float)

def extract_n_locs_from_rdf(rdf, step=10):
    """
    Estimates the expected number of localizations per nanodomain from the RDF.

    Integrates rdf * ring_area across all bins to get the mean number of
    neighbors within the search radius, which reflects domain density.
    Use this to set n_locs in spawn_nanodomains from real data.

    Args:
        rdf:  1D RDF array from extract_empirical_parameters
        step: bin width in nm (must match extraction)

    Returns:
        float: expected number of localizations per domain
    """
    rdf = np.asarray(rdf, dtype=float)
    radii = (np.arange(len(rdf)) + 1) * step
    ring_areas = np.pi * (radii**2 - (radii - step)**2)
    
    return float(np.sum(rdf * ring_areas))

def spawn_nanodomains(seeds, rdf=_default_rdf, adf=_default_adf, n_locs=50, step=10, start_label=0):
    """
    Spawns localizations around each seed using Poisson-sampled counts.

    n_locs is the mean number of localizations per seed — each seed draws its
    actual count from Poisson(n_locs), giving realistic domain-to-domain
    variability. Use expected_locs_from_rdf() to derive n_locs from real data.

    Args:
        seeds:       (N, 3) array or list of [x, y, z] seed coordinates in nm
        rdf:         1D array of weights for radial distance bins
        adf:         1D array of weights for axial offset bins
        n_locs:      mean number of localizations per seed
        step:        bin width in nm
        start_label: first cluster_label value (increment to avoid collisions
                     when concatenating multiple spawn_nanodomains calls)

    Returns:
        pd.DataFrame with columns [x [nm], y [nm], z [nm], cluster_label]
    """
    seeds = np.asarray(seeds, dtype=float)
    rdf = np.asarray(rdf, dtype=float)
    adf = np.asarray(adf, dtype=float)

    rdf_prob = rdf / rdf.sum()
    adf_prob = adf / adf.sum()

    all_coords = []

    for cluster_id, seed in enumerate(seeds, start=start_label):
        count = np.random.poisson(n_locs)
        if count == 0:
            continue

        r_bins   = np.random.choice(len(rdf), size=count, p=rdf_prob)
        z_bins   = np.random.choice(len(adf), size=count, p=adf_prob)
        r        = (r_bins + np.random.rand(count)) * step
        z_offset = (z_bins + np.random.rand(count)) * step
        theta    = np.random.rand(count) * 2 * np.pi
        z_sign   = np.random.choice([-1, 1], size=count)

        coords = np.column_stack([
            seed[0] + r * np.cos(theta),
            seed[1] + r * np.sin(theta),
            seed[2] + z_sign * z_offset,
            np.full(count, cluster_id),
        ])
        all_coords.append(coords)

    data = np.vstack(all_coords)
    return pd.DataFrame(data, columns=["x [nm]", "y [nm]", "z [nm]", "cluster_label"])