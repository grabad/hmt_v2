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


def extract_radial_density_profile(df, contour_bands, x_min, y_min, bin_size=50):
    """
    Extracts localization density per contour band from experimental data.

    For each band, divides the number of localizations in that band by the
    number of pixels in that band (a proxy for band area). The result is the
    band_density_profile argument for place_seeds — it encodes where in the
    nucleus the mark is enriched (peripheral vs. interior).

    Band numbering from create_radial_contours: 99 = outermost periphery,
    0 = innermost core, -1 = outside nucleus.

    Pass localizations from a single channel so the profile reflects only
    that mark's spatial distribution.

    Args:
        df:            DataFrame with [x [nm], y [nm]] from a single channel
        contour_bands: 2D int array from create_radial_contours
        x_min:         x origin in nm (min x used to build the nucleus mask)
        y_min:         y origin in nm
        bin_size:      pixel size in nm (must match what was used to build contour_bands)

    Returns:
        profile: 1D float array of length (max_band + 1), density per band
    """
    coords = df[["x [nm]", "y [nm]"]].to_numpy()
    x_idx = np.clip(((coords[:, 0] - x_min) / bin_size).astype(int), 0, contour_bands.shape[0] - 1)
    y_idx = np.clip(((coords[:, 1] - y_min) / bin_size).astype(int), 0, contour_bands.shape[1] - 1)
    loc_bands = contour_bands[x_idx, y_idx]

    num_bands = int(contour_bands.max()) + 1
    profile = np.zeros(num_bands, dtype=float)
    for b in range(num_bands):
        n_pixels = int(np.sum(contour_bands == b))
        n_locs = int(np.sum(loc_bands == b))
        if n_pixels > 0:
            profile[b] = n_locs / n_pixels

    return profile


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


def estimate_noise_fraction(df, contour_bands, x_min, y_min, rdf=None, sdis=200, step=10,
                            bin_size=50, enrichment_threshold=1.2, verbose=True):
    """
    Estimates the fraction of noise localizations using local density thresholding.

    For each point, counts 2D neighbors within a search radius and compares to the
    expected count under a purely random (uniform Poisson) process: rho_total * pi * r².
    Points below enrichment_threshold × that expectation are classified as noise.

    If rdf is provided, the search radius is derived automatically as the radius where
    the cumulative RDF integral reaches 90% of its total — the natural domain scale.

    Choosing enrichment_threshold:
      In SMLM data, background density is often comparable to domain density, so the
      domain signal only modestly exceeds background at the domain scale. A threshold
      of 1.1–1.3 (10–30% above random) is typically appropriate. Threshold=2.0 requires
      the domain signal to double the background, which is too strict for most SMLM marks.
      If f_noise is near 1.0, lower the threshold. If it is near 0, raise it.
      A biologically reasonable range for SMLM chromatin marks is 0.4–0.8.

    Args:
        df:                   DataFrame with [x [nm], y [nm]]
        contour_bands:        2D int array from create_radial_contours
        x_min:                x origin in nm used to build the nucleus mask
        y_min:                y origin in nm
        rdf:                  corrected RDF from extract_empirical_parameters (recommended)
        sdis:                 neighbor search radius in nm, used only if rdf is None
        step:                 RDF bin width in nm (must match extraction step)
        bin_size:             pixel size in nm (must match contour_bands)
        enrichment_threshold: multiples of random expectation below which a point is noise
        verbose:              if True, prints the auto-derived domain scale and threshold

    Returns:
        f_noise: float in [0, 1]
    """
    if rdf is not None:
        rdf_arr = np.asarray(rdf, dtype=float)
        radii = (np.arange(len(rdf_arr)) + 1) * step
        ring_areas = np.pi * (radii**2 - (radii - step)**2)
        cumulative = np.cumsum(rdf_arr * ring_areas)
        total = cumulative[-1]
        if total > 0:
            sdis = float(radii[cumulative >= 0.9 * total][0])

    coords_2d = df[["x [nm]", "y [nm]"]].to_numpy()
    A_nucleus = float(np.sum(contour_bands >= 0)) * bin_size**2
    rho_total = len(df) / A_nucleus
    count_threshold = enrichment_threshold * rho_total * np.pi * sdis**2

    if verbose:
        print(f"  domain scale: {sdis:.0f} nm | "
              f"expected random count: {rho_total * np.pi * sdis**2:.1f} | "
              f"threshold: {count_threshold:.1f} neighbors")

    tree = cKDTree(coords_2d)
    n_neighbors = np.array(tree.query_ball_point(coords_2d, sdis, return_length=True), dtype=float) - 1

    return float(np.mean(n_neighbors < count_threshold))


def make_grid_seeds(n_rows=10, n_cols=10, spacing=500.0, z_values=None, origin=(0.0, 0.0, 0.0)):
    """
    Generates a regular rectangular grid of seed coordinates for sensitivity analysis.

    Seeds are placed on a uniform grid with a fixed inter-seed spacing. Z coordinates
    are either fixed at a constant value or randomly sampled from a provided array
    (matching the behaviour of place_seeds).

    Args:
        n_rows:    number of rows in the grid
        n_cols:    number of columns in the grid
        spacing:   centre-to-centre distance between adjacent seeds in nm
        z_values:  scalar z (nm) applied to every seed, or 1-D array to sample from.
                   Defaults to 0.0 when None.
        origin:    (x0, y0, z0) offset of the grid's bottom-left corner in nm

    Returns:
        np.ndarray of shape (n_rows * n_cols, 3) with [x, y, z] in nm
    """
    x0, y0, z0 = origin

    col_idx, row_idx = np.meshgrid(np.arange(n_cols), np.arange(n_rows))
    xs = x0 + col_idx.ravel() * spacing
    ys = y0 + row_idx.ravel() * spacing

    n_seeds = n_rows * n_cols
    if z_values is None:
        zs = np.full(n_seeds, float(z0))
    elif np.ndim(z_values) == 0:
        zs = np.full(n_seeds, float(z_values))
    else:
        zs = np.random.choice(np.asarray(z_values, dtype=float), size=n_seeds)

    return np.column_stack([xs, ys, zs])


def place_seeds(contour_bands, band_density_profile, real_z_coords, x_min, y_min, px_size=50.0, scaling_factor=1.0):
    """
    Places nanodomain seeds within the nucleus using a radially-weighted Poisson process.

    Each pixel inside the nucleus draws its seed count independently from
    Poisson(density * scaling_factor), where density comes from band_density_profile.
    This correctly allows high-density pixels to host multiple seeds and zero-density
    pixels to host none. Sub-pixel jitter prevents seeds from landing on a pixel-center grid.

    Z-coordinates are sampled from the empirical distribution of real_z_coords rather
    than inferred from the 2D mask, which captures the true axial extent of the nucleus.

    Args:
        contour_bands:        2D int array from create_radial_contours (-1 = outside nucleus)
        band_density_profile: 1D float array from extract_radial_density_profile
        real_z_coords:        1D array of z values from real localizations (nm), sampled for z
        x_min:                x origin in nm used to build the nucleus mask
        y_min:                y origin in nm
        px_size:              pixel size in nm (must match bin_size in contour_bands)
        scaling_factor:       multiplier on seed density; tune to match desired domain count

    Returns:
        np.ndarray of shape (N, 3) — seed [x, y, z] coordinates in nm
    """
    band_density_profile = np.asarray(band_density_profile, dtype=float)

    inside = contour_bands >= 0
    prob_map = np.zeros_like(contour_bands, dtype=float)
    prob_map[inside] = band_density_profile[contour_bands[inside]]

    seed_counts = np.random.poisson(prob_map * scaling_factor)
    pixel_indices = np.argwhere(seed_counts > 0)
    if len(pixel_indices) == 0:
        return np.empty((0, 3), dtype=float)

    counts = seed_counts[seed_counts > 0]
    seed_indices = np.repeat(pixel_indices, counts, axis=0)

    jitter = np.random.rand(*seed_indices.shape)
    xy = (seed_indices + jitter) * px_size
    xy[:, 0] += x_min
    xy[:, 1] += y_min

    z = np.random.choice(np.asarray(real_z_coords, dtype=float), size=len(xy))

    return np.column_stack([xy, z])


def add_sim_noise(sim_df, noise_fraction=0.1):
    """
    Adds uniformly distributed noise localizations within the bounding box of a simulation.

    The noise count scales with the simulation size, making it straightforward to
    sweep noise levels for sensitivity analysis without requiring a nucleus mask.
    All noise points are assigned cluster_label == -1.

    Args:
        sim_df:         DataFrame from spawn_nanodomains with [x [nm], y [nm], z [nm], cluster_label]
        noise_fraction: fraction of len(sim_df) to generate as noise (e.g. 0.2 = 20%)

    Returns:
        pd.DataFrame with columns [x [nm], y [nm], z [nm], cluster_label],
        containing only the noise rows (cluster_label == -1 for all)
    """
    n_noise = int(noise_fraction * len(sim_df))
    if n_noise == 0:
        return pd.DataFrame(columns=sim_df.columns)

    xs = np.random.uniform(sim_df["x [nm]"].min(), sim_df["x [nm]"].max(), size=n_noise)
    ys = np.random.uniform(sim_df["y [nm]"].min(), sim_df["y [nm]"].max(), size=n_noise)
    zs = np.random.uniform(sim_df["z [nm]"].min(), sim_df["z [nm]"].max(), size=n_noise)

    return pd.DataFrame({
        "x [nm]": xs,
        "y [nm]": ys,
        "z [nm]": zs,
        "cluster_label": np.full(n_noise, -1, dtype=float),
    })


def add_noise_locs(n_noise, contour_bands, real_z_coords, x_min, y_min, px_size=50.0):
    """
    Adds uniformly distributed background localizations within the nucleus.

    These represent the non-domain fraction of the data — points that are not part
    of any nanodomain and were removed from the RDF by baseline subtraction.
    Estimate n_noise as int(f_noise * len(real_df)), where f_noise comes from DBSCAN
    (fraction of points with label == -1).

    Args:
        n_noise:       number of noise localizations to generate
        contour_bands: 2D int array from create_radial_contours (-1 = outside nucleus)
        real_z_coords: 1D array of z values from real localizations (nm), sampled for z
        x_min:         x origin in nm used to build the nucleus mask
        y_min:         y origin in nm
        px_size:       pixel size in nm

    Returns:
        pd.DataFrame with columns [x [nm], y [nm], z [nm], cluster_label]
        cluster_label is -1 for all noise points
    """
    inside_indices = np.argwhere(contour_bands >= 0)
    chosen = inside_indices[np.random.randint(0, len(inside_indices), size=n_noise)]

    jitter = np.random.rand(*chosen.shape)
    xy = (chosen + jitter) * px_size
    xy[:, 0] += x_min
    xy[:, 1] += y_min

    z = np.random.choice(np.asarray(real_z_coords, dtype=float), size=n_noise)

    return pd.DataFrame(
        np.column_stack([xy, z, np.full(n_noise, -1)]),
        columns=["x [nm]", "y [nm]", "z [nm]", "cluster_label"],
    )


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