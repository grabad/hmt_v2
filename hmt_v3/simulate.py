import napari
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import napari
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score
from scipy.ndimage import gaussian_filter, binary_fill_holes, label, distance_transform_edt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.qhull import QhullError
from scipy.optimize import minimize_scalar
from scipy.signal import fftconvolve

_default_bins = np.arange(100)
# rdf is a per-ring density (locs/nm^2); spawn_nanodomains weights it by ring area.
_default_rdf = np.exp(-0.5 * (_default_bins / 15.0) ** 2)
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


def add_noise_locs(n_noise, contour_bands, real_z_coords, x_min, y_min, px_size=50.0,
                   band_weights=None):
    """
    Adds background localizations within the nucleus.

    These represent the non-domain fraction of the data — points that are not part
    of any nanodomain and were removed from the RDF by baseline subtraction.
    Estimate n_noise as int(f_noise * len(real_df)), where f_noise comes from
    estimate_noise_fraction.

    By default noise is uniform over every pixel inside the mask. Because
    binarize_nucleus fills holes, that includes signal-free regions such as
    nucleoli, so uniform noise fills places the mark never occupies while the
    profile-weighted domains do not. Pass band_weights (e.g. the radial density
    profile, or a 0/1 support mask over bands) to restrict/shape the noise to the
    same spatial support as the domains.

    Args:
        n_noise:       number of noise localizations to generate
        contour_bands: 2D int array from create_radial_contours (-1 = outside nucleus)
        real_z_coords: 1D array of z values from real localizations (nm), sampled for z
        x_min:         x origin in nm used to build the nucleus mask
        y_min:         y origin in nm
        px_size:       pixel size in nm
        band_weights:  optional 1D array indexed by contour band; a pixel's
                       sampling weight is band_weights[its band]. None = uniform.

    Returns:
        pd.DataFrame with columns [x [nm], y [nm], z [nm], cluster_label]
        cluster_label is -1 for all noise points
    """
    inside_indices = np.argwhere(contour_bands >= 0)

    if band_weights is None:
        chosen = inside_indices[np.random.randint(0, len(inside_indices), size=n_noise)]
    else:
        band_weights = np.asarray(band_weights, dtype=float)
        w = band_weights[contour_bands[inside_indices[:, 0], inside_indices[:, 1]]]
        if w.sum() <= 0:
            chosen = inside_indices[np.random.randint(0, len(inside_indices), size=n_noise)]
        else:
            sel = np.random.choice(len(inside_indices), size=n_noise, p=w / w.sum())
            chosen = inside_indices[sel]

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
        rdf:         1D per-ring neighbour density (locs/nm²), as returned by
                     extract_empirical_parameters; weighted by ring area internally
        adf:         1D histogram of |z-offset| counts for axial offset bins
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

    # rdf is a neighbour DENSITY (locs/nm^2) per ring, matching
    # extract_empirical_parameters. The number of points at radius r is
    # density * ring area, so weight by ring area before sampling the radial bin.
    # Sampling on the raw density over-concentrates points near the centre and
    # makes domains too clustered at short range.
    ring_areas = np.pi * step ** 2 * (2 * np.arange(len(rdf)) + 1)
    rdf_prob = rdf * ring_areas
    rdf_prob = rdf_prob / rdf_prob.sum()
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


def epsilon_cost_cluster(eps, coords, true_labels, min_samples=8, use_z=False):
    """
    Cost function for epsilon optimisation: returns 1 - ARI.

    Runs DBSCAN on XY (use_z=False) or XYZ (use_z=True) coords and compares
    predicted labels against ground-truth cluster_label values using ARI.
    ARI = 1 means perfect agreement; returning 1 - ARI converts it to a
    cost that scalar minimisers can drive toward zero.

    Noise points (label == -1 in either array) are treated as their own
    cluster by ARI, which penalises both over-splitting and over-merging.

    Args:
        eps:         DBSCAN neighbourhood radius in nm
        coords:      (N, 3) array of [x, y, z] in nm
        true_labels: (N,) array of ground-truth cluster IDs (cluster_label column)
        min_samples: DBSCAN min_samples parameter
        use_z:       if True, cluster in 3D (XYZ); if False, cluster in 2D (XY only)

    Returns:
        float in [0, 2]: 1 - ARI (lower is better)
    """
    fit_coords = coords[:, :3] if use_z else coords[:, :2]
    labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(fit_coords)
    if np.all(labels == -1):
        return 1.0
    return 1.0 - adjusted_rand_score(true_labels, labels)


def optimize_epsilon(fraction_dfs, calibration_area_nm2, min_samples=8, use_z=False, show_plot=True):
    """
    Finds the optimal DBSCAN epsilon for each density fraction by maximising
    correct cluster assignment (ARI against ground-truth cluster_label).

    Each DataFrame in fraction_dfs must have columns
    [x [nm], y [nm], z [nm], cluster_label], where cluster_label is the
    ground-truth assigned by spawn_nanodomains (use -1 for noise from
    add_noise_locs). The function sweeps epsilon via a coarse grid then
    refines with bounded scalar minimisation, accepts the result if
    ARI > 0.3, and fits a power law eps = coef * density^exp in log-log
    space across all accepted fractions.

    Args:
        fraction_dfs:          list of DataFrames, one per density fraction
        calibration_area_nm2:  nucleus area in nm² used to compute density
        min_samples:           DBSCAN min_samples parameter
        use_z:                 if True, cluster in 3D (XYZ); if False, cluster in 2D (XY only)
        show_plot:             if True, plots log-log and linear density vs epsilon

    Returns:
        coef (float), exp (float), densities (list), best_epsilons (list)
        where eps = coef * density^exp
    """
    densities = []
    best_epsilons = []

    dims = "3D (XYZ)" if use_z else "2D (XY)"
    print(f"Starting ARI-based epsilon optimisation over {len(fraction_dfs)} fractions [{dims}]...")

    for i, df in enumerate(fraction_dfs):
        coords = df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
        true_labels = df["cluster_label"].to_numpy()

        if len(coords) < min_samples:
            print(f"Fraction {i+1}: Skipped (too few localisations)")
            continue

        density = len(coords) / calibration_area_nm2

        # --- Stage 1: coarse grid search ---
        coarse_grid = np.arange(10, 225, 15.0)
        best_coarse_eps = coarse_grid[0]
        best_coarse_cost = 2.0

        for test_eps in coarse_grid:
            cost = epsilon_cost_cluster(test_eps, coords, true_labels, min_samples, use_z)
            if cost < best_coarse_cost:
                best_coarse_cost = cost
                best_coarse_eps = test_eps

        # --- Stage 2: bounded scalar minimisation within ±15 nm of coarse winner ---
        search_min = max(10.0, best_coarse_eps - 15.0)
        search_max = min(225.0, best_coarse_eps + 15.0)

        res = minimize_scalar(
            epsilon_cost_cluster,
            bounds=(search_min, search_max),
            args=(coords, true_labels, min_samples, use_z),
            method='bounded',
            options={'xatol': 0.2},
        )

        if res.success and (1.0 - res.fun) > 0.3:
            ari = 1.0 - res.fun
            densities.append(density)
            best_epsilons.append(res.x)
            print(f"Fraction {i+1}: Density {density:.6e} -> Eps {res.x:.2f} nm  (ARI {ari:.4f})")
        else:
            ari_str = f"{1.0 - res.fun:.4f}" if res.success else "opt failed"
            print(f"Fraction {i+1}: Density {density:.6e} | FAILED (ARI {ari_str})")

    if len(densities) < 2:
        print("\nError: not enough valid fractions for power-law fit.")
        return 0, 0, densities, best_epsilons

    # --- Fit power law in log-log space ---
    log_dens = np.log10(densities)
    log_eps = np.log10(best_epsilons)
    par = np.polyfit(log_dens, log_eps, 1)

    exp = par[0]
    coef = 10 ** par[1]
    r_squared = np.corrcoef(log_dens, log_eps)[0, 1] ** 2

    print(f"\nCalibration complete: Eps = {coef:.2f} * Density^{exp:.2f}  (R² = {r_squared:.4f})")

    if show_plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        fit_line = par[0] * np.array(log_dens) + par[1]
        axes[0].scatter(log_dens, log_eps, color='blue', label='Optimised epsilons')
        axes[0].plot(log_dens, fit_line, color='red', linestyle='--',
                     label=f'Fit: eps = {coef:.2f} · dens^{exp:.2f}')
        axes[0].set_xlabel('Log10(Density [locs/nm²])')
        axes[0].set_ylabel('Log10(Epsilon [nm])')
        axes[0].set_title('Log-Log Scale')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        dens_arr = np.array(densities)
        dens_smooth = np.linspace(dens_arr.min(), dens_arr.max(), 200)
        axes[1].scatter(dens_arr, best_epsilons, color='blue', label='Optimised epsilons')
        axes[1].plot(dens_smooth, coef * dens_smooth ** exp, color='red', linestyle='--',
                     label=f'Fit: eps = {coef:.2f} · dens^{exp:.2f}')
        axes[1].set_xlabel('Density [locs/nm²]')
        axes[1].set_ylabel('Epsilon [nm]')
        axes[1].set_title('Standard Scale')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.suptitle(f'DBSCAN Epsilon Calibration (ARI-optimised)\nR² = {r_squared:.4f}')
        plt.tight_layout()
        plt.show()

    return coef, exp, densities, best_epsilons


def measure_seed_spacing(df, eps=100.0, min_samples=8, use_z=False):
    """
    Measures the inter-domain spatial arrangement from real localizations.

    NOTE: this is density-dependent and is kept for reference only. It clusters
    with DBSCAN, so the domain count and centroid positions shift with labelling
    density (mainly through min_samples). For a density-independent measurement,
    use measure_pair_correlation (Ripley K/L and g(r)), which generate_nucleus
    now uses instead.

    Clusters localizations with DBSCAN to locate nanodomain centres, then
    computes the nearest-neighbour distance (NND) distribution between those
    centres. This captures how domains are arranged relative to one another
    (regular/dispersed vs. clustered) — the structure that a purely radial,
    independent-Poisson seed placement (place_seeds) does not reproduce.

    Pass in a single channel so the centres reflect that mark's domains only.
    A reasonable eps is the domain scale (e.g. the radius where the cumulative
    RDF reaches ~90% of its total); generate_nucleus derives this automatically.

    Args:
        df:          DataFrame with [x [nm], y [nm], z [nm]] from a single channel
        eps:         DBSCAN neighbourhood radius in nm (the domain scale)
        min_samples: DBSCAN min_samples parameter
        use_z:       cluster in 3D (XYZ) if True, else 2D (XY)

    Returns:
        centers:   (M, 3) array of domain centroid [x, y, z] in nm
        nnd:       (M,) array of 2D nearest-neighbour distances between centres (nm),
                   empty if fewer than two domains are found
        n_domains: int, number of domains found
    """
    coords = df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
    fit_coords = coords if use_z else coords[:, :2]
    labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(fit_coords)

    ids = np.unique(labels[labels >= 0])
    if len(ids) == 0:
        return np.empty((0, 3), dtype=float), np.array([]), 0

    centers = np.array([coords[labels == k].mean(axis=0) for k in ids])
    n_domains = len(centers)

    if n_domains < 2:
        return centers, np.array([]), n_domains

    tree = cKDTree(centers[:, :2])
    d, _ = tree.query(centers[:, :2], k=2)  # column 0 is the self-match
    nnd = d[:, 1]

    return centers, nnd, n_domains


def place_seeds_matched(contour_bands, band_density_profile, real_z_coords,
                        x_min, y_min, n_seeds, target_nnd=None, min_dist=None,
                        px_size=50.0, batch=2048, max_attempts_per_seed=2000, rng=None):
    """
    Places seeds matching BOTH the radial density profile and the inter-domain
    spacing distribution of the real nucleus.

    Candidates are proposed proportional to the radial band density — so the
    nuclear radial distribution is preserved exactly as in place_seeds — then
    accepted through a sequential-inhibition rule so the resulting
    nearest-neighbour distances reproduce the measured domain arrangement:

      * target_nnd given (empirical NND array from measure_seed_spacing): each
        candidate must lie at least a randomly drawn spacing t ~ target_nnd from
        every already-placed seed. Bootstrapping the required spacing from the
        empirical NND reproduces the full distribution, capturing dispersion or
        clustering as observed — not just a single exclusion radius.
      * only min_dist given: a fixed hard-core radius is enforced (classic
        Matern / simple sequential inhibition, purely dispersive).
      * neither given: reduces to an independent radial Poisson placement.

    Z-coordinates are sampled from real_z_coords, matching place_seeds.

    Args:
        contour_bands:        2D int array from create_radial_contours (-1 outside)
        band_density_profile: 1D float array from extract_radial_density_profile
        real_z_coords:        1D array of real z values (nm), sampled for seed z
        x_min, y_min:         mask origin in nm (min x/y used to build the mask)
        n_seeds:              number of seeds to place (e.g. n_domains from real data)
        target_nnd:           empirical NND array; matches the full spacing distribution
        min_dist:             fixed hard-core spacing in nm; used only if target_nnd is None
        px_size:              pixel size in nm (must match contour_bands)
        batch:                number of candidates proposed per vectorized draw
        max_attempts_per_seed:candidate budget scales as this * n_seeds
        rng:                  optional np.random.Generator for reproducibility

    Returns:
        (N, 3) array of seed [x, y, z] in nm. N <= n_seeds if the inhibition
        constraint cannot pack n_seeds within the attempt budget.
    """
    rng = np.random.default_rng() if rng is None else rng
    band_density_profile = np.asarray(band_density_profile, dtype=float)
    real_z_coords = np.asarray(real_z_coords, dtype=float)

    inside = contour_bands >= 0
    prob_map = np.zeros(contour_bands.shape, dtype=float)
    prob_map[inside] = band_density_profile[contour_bands[inside]]
    if prob_map.sum() <= 0 or n_seeds < 1:
        return np.empty((0, 3), dtype=float)

    flat_prob = prob_map.ravel()
    flat_idx = np.flatnonzero(flat_prob)
    flat_p = flat_prob[flat_idx]
    flat_p /= flat_p.sum()
    n_cols = contour_bands.shape[1]

    target_nnd = None if target_nnd is None else np.asarray(target_nnd, dtype=float)
    use_target = target_nnd is not None and len(target_nnd) > 0

    placed_xy = np.empty((n_seeds, 2), dtype=float)
    placed_z = np.empty(n_seeds, dtype=float)
    n_placed = 0
    attempts = 0
    attempt_budget = max_attempts_per_seed * n_seeds

    while n_placed < n_seeds and attempts < attempt_budget:
        choice = rng.choice(flat_idx, size=batch, p=flat_p)
        rows = choice // n_cols
        cols = choice % n_cols
        jitter = rng.random((batch, 2))
        cand_x = (rows + jitter[:, 0]) * px_size + x_min
        cand_y = (cols + jitter[:, 1]) * px_size + y_min

        for i in range(batch):
            attempts += 1

            if n_placed == 0:
                accept = True
            else:
                dx = placed_xy[:n_placed, 0] - cand_x[i]
                dy = placed_xy[:n_placed, 1] - cand_y[i]
                nearest = np.sqrt((dx * dx + dy * dy).min())
                if use_target:
                    required = target_nnd[rng.integers(len(target_nnd))]
                elif min_dist is not None:
                    required = min_dist
                else:
                    required = 0.0
                accept = nearest >= required

            if accept:
                placed_xy[n_placed] = (cand_x[i], cand_y[i])
                placed_z[n_placed] = real_z_coords[rng.integers(len(real_z_coords))]
                n_placed += 1
                if n_placed >= n_seeds:
                    break

    return np.column_stack([placed_xy[:n_placed], placed_z[:n_placed]])


def measure_pair_correlation(df, binary_mask, x_min, y_min, r_max=1000.0, dr=20.0,
                             px_size=50.0):
    """
    Measures the density-independent spatial clustering of localizations via
    Ripley's K/L functions and the pair correlation function g(r).

    Unlike DBSCAN-based domain finding, these statistics are computed directly
    from the localizations and are invariant to random under-labelling: under
    independent thinning (the model of reduced staining efficiency) K(r), L(r)
    and g(r) are unchanged in expectation, because each point's neighbour count
    and the intensity lambda both scale with the retained fraction and cancel.
    This makes them the correct tool for comparing domain structure across cells
    of differing staining density, and the density-safe replacement for
    measure_seed_spacing.

    Edge effects from the irregular nuclear boundary are corrected per point by
    normalising each point's neighbour count to the fraction of its search disk
    that actually falls inside the nucleus mask (local-area / translation
    correction), estimated by convolving the mask with disk kernels.

    Interpretation:
        L(r) - r > 0  -> clustering at scale r (excess neighbours vs. CSR)
        L(r) - r < 0  -> regularity / exclusion at scale r
        g(r) > 1      -> pair density above complete spatial randomness at r
        g(r) ~ 1      -> uncorrelated at r (returns to 1 beyond the domain scale)

    Pass a single channel so the statistics reflect that mark only.

    Args:
        df:           DataFrame with [x [nm], y [nm]] from a single channel
        binary_mask:  2D bool/int nucleus mask (True inside); contour_bands >= 0 works
        x_min, y_min: mask origin in nm (min x/y used to build the mask)
        r_max:        maximum radius in nm
        dr:           radial bin width in nm
        px_size:      pixel size in nm (must match binary_mask)

    Returns:
        dict with keys:
            'r':      (M,) radii in nm
            'K':      (M,) Ripley's K(r)          [CSR = pi r^2]
            'L':      (M,) Ripley's L(r)=sqrt(K/pi) [CSR = r]
            'g':      (M,) pair correlation g(r)   [CSR = 1]
            'lambda': intensity (localizations / nm^2)
            'area':   nucleus area in nm^2
    """
    mask = np.asarray(binary_mask) > 0
    area = float(mask.sum()) * px_size ** 2

    coords = df[["x [nm]", "y [nm]"]].to_numpy()
    n = len(coords)
    lam = n / area

    radii = np.arange(dr, r_max + dr, dr)

    # Map each point to its mask pixel (same convention as extract_radial_density_profile)
    x_idx = np.clip(((coords[:, 0] - x_min) / px_size).astype(int), 0, mask.shape[0] - 1)
    y_idx = np.clip(((coords[:, 1] - y_min) / px_size).astype(int), 0, mask.shape[1] - 1)

    tree = cKDTree(coords)
    mask_f = mask.astype(float)

    K = np.empty(len(radii))
    for k, r in enumerate(radii):
        # observed neighbours within r (exclude self)
        counts = np.asarray(tree.query_ball_point(coords, r, return_length=True), dtype=float) - 1

        # fraction of each point's search disk that lies inside the mask
        rpix = r / px_size
        krad = int(np.ceil(rpix))
        yy, xx = np.mgrid[-krad:krad + 1, -krad:krad + 1]
        disk = ((xx ** 2 + yy ** 2) <= rpix ** 2).astype(float)
        inside_area = fftconvolve(mask_f, disk, mode="same") * px_size ** 2  # nm^2 inside mask
        frac = np.clip(inside_area[x_idx, y_idx] / (np.pi * r ** 2), 1e-3, 1.0)

        # translation-corrected Ripley's K
        K[k] = (area / n ** 2) * np.sum(counts / frac)

    L = np.sqrt(np.clip(K, 0, None) / np.pi)

    # g(r) = K'(r) / (2 pi r), by finite difference
    r_mid = radii - dr / 2
    dK = np.diff(np.concatenate([[0.0], K]))
    g = dK / (2 * np.pi * r_mid * dr)

    return {"r": radii, "K": K, "L": L, "g": g, "lambda": lam, "area": area}


def estimate_n_domains_from_g(pcf, r_domain=None):
    """
    Estimates the number of nanodomains from the pair correlation function,
    density-independently, under a Neyman-Scott (Poisson-cluster) model.

    For a Poisson cluster process with parent (domain) intensity kappa, the pair
    correlation is g(r) = 1 + (1/kappa) f(r), where f integrates to 1 over the
    plane. Hence kappa = 1 / integral[(g(r) - 1) 2 pi r dr] over the within-domain
    range, and the domain count is kappa * area. Because g(r) is thinning-invariant
    (see measure_pair_correlation), this count does not depend on labelling density
    — the property DBSCAN-based counting lacks.

    Assumes domains are sparse and their centres are approximately Poisson. If
    domains are strongly ordered, g(r) dips below 1 between domains and this will
    over-count; check L(r) - r for a dispersion trough (see match dispersion in
    generate_nucleus).

    Args:
        pcf:      dict from measure_pair_correlation
        r_domain: integration cut-off in nm (where g returns to ~1). If None,
                  the first radius at which g(r) <= 1 is used.

    Returns:
        n_domains: estimated domain count (int, >= 1)
        kappa:     estimated domain intensity (domains / nm^2)
        r_domain:  cut-off radius used (nm)
    """
    r = pcf["r"]
    g = pcf["g"]
    area = pcf["area"]

    if r_domain is None:
        below = np.where(g <= 1.0)[0]
        r_domain = float(r[below[0]]) if len(below) else float(r[-1])

    sel = r <= r_domain
    y = (g[sel] - 1.0) * 2 * np.pi * r[sel]
    x = r[sel]
    integral = float(np.sum((y[1:] + y[:-1]) / 2.0 * np.diff(x))) if len(x) > 1 else 0.0
    if integral <= 0:
        return 1, 1.0 / area, r_domain

    kappa = 1.0 / integral
    n_domains = max(1, int(round(kappa * area)))
    return n_domains, kappa, r_domain


def estimate_noise_fraction_from_pcf(pcf_real, pcf_reference, r_fit_max=None, excess_frac=0.2):
    """
    Estimates the background/noise fraction density-independently, from the pair
    correlation function.

    Adding a fraction (1 - phi) of uniform background to a clustered (domain)
    process scales its correlation excess by phi^2:
        g_total(r) - 1 = phi^2 (g_domain(r) - 1),   phi = domain fraction.
    If the reference is a noise-free domain simulation built from ALL real
    localizations (domain fraction 1), the Neyman-Scott amplitude (which scales as
    1 / n_domains) cancels one factor of phi, so within the domain scale the real
    data's excess is phi times the reference's:
        (g_real(r) - 1) = phi * (g_reference(r) - 1).
    Hence phi is the least-squares slope of (g_real - 1) vs (g_reference - 1), and
    the noise fraction is 1 - phi. Because g(r) is thinning-invariant, so is this
    estimate — unlike the local-density heuristic in estimate_noise_fraction.

    IMPORTANT: the fit must be restricted to the WITHIN-domain range. The reference
    places domain centres as radial Poisson, so it has no inter-domain clustering;
    real data usually does (multi-scale / domain-of-domains structure). Beyond the
    domain scale real is MORE clustered than the reference, which would push phi >= 1
    and collapse the estimate to zero noise. Using the differential g(r) - 1 (which
    returns to 0 beyond the domain) plus an r_fit_max cutoff keeps the fit local.

    Args:
        pcf_real:      dict from measure_pair_correlation (real cell)
        pcf_reference: dict from measure_pair_correlation of a noise-free domain
                       simulation built from all localizations (domain fraction 1)
        r_fit_max:     restrict the fit to r <= r_fit_max (nm); pass the domain scale.
                       None = no radial cap (rely on excess_frac alone).
        excess_frac:   fit only bins where the reference excess exceeds this fraction
                       of its maximum (restricts to the clustered range)

    Returns:
        f_noise: estimated noise fraction in [0, 1]
    """
    r = pcf_real["r"]
    a = pcf_real["g"] - 1.0        # differential excess of the real data
    b = pcf_reference["g"] - 1.0   # differential excess of the reference

    sel = b > excess_frac * np.nanmax(b)
    if r_fit_max is not None:
        sel = sel & (r <= r_fit_max)
    if not np.any(sel) or np.nansum(b[sel] ** 2) <= 0:
        return 0.0

    phi = float(np.nansum(a[sel] * b[sel]) / np.nansum(b[sel] ** 2))  # LS slope through origin
    phi = min(max(phi, 0.0), 1.0)
    return 1.0 - phi


def fit_thomas_parameters(pcf_real, domain_scale, f_noise, r_max_fit=None):
    """
    Fits a Thomas (Neyman-Scott) cluster process to the inter-domain part of the
    real pair correlation, for clustered seed placement.

    Model — the pair correlation of seeds (domain centres) drawn as offspring of
    Poisson parents with a Gaussian dispersal kernel:
        g_seed(r) - 1 = exp(-r^2 / (4 sigma^2)) / (4 pi kappa_p sigma^2),
    where kappa_p is the parent intensity and sigma the offspring kernel width.

    The measured localization g(r) mixes the within-domain profile (below
    domain_scale) with this inter-domain term (above it), and uniform noise dilutes
    the whole excess by phi^2 (phi = 1 - f_noise). So the fit uses
    (g_real(r) - 1) / phi^2 for r > domain_scale, which isolates g_seed(r) - 1.
    Taking logs makes it linear in r^2: ln(excess) = ln(A) - r^2 / (4 sigma^2), so a
    straight-line fit recovers sigma (slope) and A = 1/(4 pi kappa_p sigma^2)
    (intercept). g(r) is thinning-invariant, so the fit is density-independent.

    Args:
        pcf_real:     dict from measure_pair_correlation (real cell)
        domain_scale: within-domain cutoff in nm; the fit uses r > domain_scale
        f_noise:      noise fraction (to undo the phi^2 dilution)
        r_max_fit:    upper radius for the fit in nm (default: last pcf radius)

    Returns:
        kappa_p: parent intensity (parents / nm^2), or None if no clustered signal
        sigma:   offspring kernel width in nm, or None
        ok:      True if a decaying clustered signal was fit; if False, callers should
                 fall back to Poisson placement
    """
    r = pcf_real["r"]
    phi2 = max((1.0 - f_noise) ** 2, 1e-6)
    excess = (pcf_real["g"] - 1.0) / phi2

    hi = r[-1] if r_max_fit is None else r_max_fit
    sel = (r > domain_scale) & (r <= hi) & (excess > 0)
    if np.count_nonzero(sel) < 2:
        return None, None, False

    slope, intercept = np.polyfit(r[sel] ** 2, np.log(excess[sel]), 1)
    if slope >= 0:  # excess grows with r: not a decaying cluster kernel
        return None, None, False

    sigma = float(np.sqrt(-1.0 / (4.0 * slope)))
    A = float(np.exp(intercept))
    if A <= 0 or not np.isfinite(sigma) or sigma <= 0:
        return None, None, False

    kappa_p = 1.0 / (4.0 * np.pi * A * sigma ** 2)
    return kappa_p, sigma, True


def place_seeds_thomas(contour_bands, band_density_profile, real_z_coords, x_min, y_min,
                       n_seeds, kappa_p, sigma, px_size=50.0, rng=None):
    """
    Places nanodomain seeds as a Thomas (Neyman-Scott) cluster process, reproducing
    inter-domain clustering while still following the radial density profile.

    Parents are placed with the radial weighting (place_seeds_matched); each parent
    scatters offspring domain centres with a Gaussian(sigma) displacement in xy.
    Offspring falling outside the nucleus mask are discarded, so the returned count
    can be slightly below n_seeds. z is sampled from real_z_coords.

    Fit kappa_p and sigma with fit_thomas_parameters. Use this in place of
    place_seeds_matched when the real cell shows inter-domain clustering (elevated
    L(r) - r beyond the domain scale).

    Args:
        contour_bands:        2D int array from create_radial_contours (-1 outside)
        band_density_profile: 1D float array from extract_radial_density_profile
        real_z_coords:        1D array of real z values (nm), sampled for seed z
        x_min, y_min:         mask origin in nm
        n_seeds:              target number of offspring seeds (e.g. n_domains)
        kappa_p:              parent intensity (parents / nm^2)
        sigma:                offspring kernel width in nm
        px_size:              pixel size in nm (must match contour_bands)
        rng:                  optional np.random.Generator for reproducibility

    Returns:
        (N, 3) array of seed [x, y, z] in nm (N <= n_seeds after mask rejection)
    """
    rng = np.random.default_rng() if rng is None else rng
    real_z_coords = np.asarray(real_z_coords, dtype=float)
    mask = contour_bands >= 0
    area = float(mask.sum()) * px_size ** 2

    n_parents = max(1, int(round(kappa_p * area)))
    parents = place_seeds_matched(contour_bands, band_density_profile, real_z_coords,
                                  x_min, y_min, n_seeds=n_parents, px_size=px_size, rng=rng)
    if len(parents) == 0:
        return np.empty((0, 3), dtype=float)

    mean_off = n_seeds / len(parents)
    counts = rng.poisson(mean_off, size=len(parents))
    total = int(counts.sum())
    if total == 0:
        return np.empty((0, 3), dtype=float)

    centers = np.repeat(parents[:, :2], counts, axis=0)
    xy = centers + rng.normal(0.0, sigma, size=(total, 2))

    xi = ((xy[:, 0] - x_min) / px_size).astype(int)
    yi = ((xy[:, 1] - y_min) / px_size).astype(int)
    within = (xi >= 0) & (xi < mask.shape[0]) & (yi >= 0) & (yi < mask.shape[1])
    inside = np.zeros(total, dtype=bool)
    inside[within] = mask[xi[within], yi[within]]
    xy = xy[inside]
    if len(xy) == 0:
        return np.empty((0, 3), dtype=float)

    z = rng.choice(real_z_coords, size=len(xy))
    return np.column_stack([xy, z])


def generate_nucleus(df, contour_bands, x_min, y_min, min_samples=8, sdis=500,
                     step=10, px_size=50.0, match_spacing=True, noise_mode="support",
                     noise_fraction=None, domain_scale=None, seed_placement="thomas",
                     rng=None, verbose=True):
    """
    Reconstructs a simulated single-channel nucleus that matches a real cell in
    three respects, using only density-independent measurements:

      1. the nuclear radial distribution of localizations (extract_radial_density_profile),
      2. the within-domain point distribution (RDF/ADF from extract_empirical_parameters), and
      3. the inter-domain arrangement, characterised by the pair correlation
         function / Ripley's L (measure_pair_correlation) rather than DBSCAN, so
         nothing here depends on labelling density.

    Domains are built first; uniform background noise is generated separately so
    the two components can be inspected on their own.

    The domain count is estimated from g(r) under a Poisson-cluster model
    (estimate_n_domains_from_g) — thinning-invariant, unlike a DBSCAN count. The
    per-domain localization count is set so the domain-only total equals
    (1 - f_noise) times the real localization count; noise adds the remaining
    f_noise fraction, so the noisy simulation matches the real cell size. Domain
    centres are placed by an inhibition process (place_seeds_matched) whose
    hard-core radius is read from any dispersion trough in the real L(r) - r; if
    the real cell shows no dispersion, centres are placed as radial Poisson.

    Args:
        df:            DataFrame with [x [nm], y [nm], z [nm]] from a single channel
        contour_bands: 2D int array from create_radial_contours (mask = bands >= 0)
        x_min, y_min:  mask origin in nm. MUST be the origin the mask grid was built
                       on — the min x/y of the PRE-mask dataframes passed to
                       binarize_nucleus, not the masked outputs. Use
                       preprocess.mask_origin(me3_df, ac_df) on the pre-mask inputs
                       to get it; a wrong origin shifts localizations into the wrong
                       bands and empties the core.
        min_samples:   unused; retained for backward compatibility
        sdis:          max search radius in nm for RDF/ADF extraction
        step:          bin width in nm
        px_size:       pixel size in nm (must match contour_bands)
        match_spacing: if True, enforce a hard-core spacing derived from L(r) - r;
                       if False, place centres as radial Poisson
        noise_mode:    where background noise may fall.
                       "support" (default): uniform over bands that contain real
                          signal (profile > 0), so noise avoids signal-free holes
                          such as filled nucleoli — the domains and noise then share
                          the same radial support.
                       "profile": proportional to the radial density profile, so
                          noise follows the mark's spatial distribution.
                       "uniform": uniform over the whole mask (original behaviour;
                          fills the empty centre).
        noise_fraction:explicit noise fraction in [0, 1]. If None (default), it is
                       estimated from g(r) via estimate_noise_fraction_from_pcf
                       (density-independent). Set it to override when you want to
                       dial the background up or down.
        domain_scale:  explicit domain scale in nm (the within-domain cutoff radius;
                       ~half the domain diameter). If None (default), it is read from
                       g(r), which can underestimate; set it to a known/measured
                       domain radius (e.g. ~100 nm for 200 nm domains).
        seed_placement:how domain centres are arranged.
                       "thomas" (default): Neyman-Scott cluster process fit to the
                          inter-domain g(r) (fit_thomas_parameters), reproducing
                          domain-of-domains clustering. Falls back to Poisson if no
                          clustered signal is found.
                       "inhibition": radial placement with a hard-core spacing read
                          from an L(r) - r dispersion trough (needs match_spacing).
                       "poisson": independent radial placement.
        rng:           optional np.random.Generator for reproducibility
        verbose:       print a short summary of the fitted parameters

    Returns:
        domains_df: DataFrame [x, y, z, cluster_label] — domains only, no noise
        noisy_df:   DataFrame [x, y, z, cluster_label] — domains + noise (noise label -1)
        info:       dict of measured/derived parameters, seeds, the fitted noise
                    fractions (f_noise, f_noise_from_g, f_noise_heuristic), and the
                    pair-correlation curves pcf_real, pcf_domains, pcf_noisy,
                    pcf_reference. Compare pcf_noisy against pcf_real to judge noise.
    """
    rng = np.random.default_rng() if rng is None else rng
    mask = contour_bands >= 0

    rdf, adf = extract_empirical_parameters(df, sdis=sdis, step=step)
    profile = extract_radial_density_profile(df, contour_bands, x_min, y_min, bin_size=px_size)

    # Density-independent inter-domain measurement (for arrangement + validation).
    pcf_real = measure_pair_correlation(df, mask, x_min, y_min, r_max=sdis, dr=2 * step, px_size=px_size)

    # Domain scale: the radius where the pair correlation returns to ~1, i.e. the
    # separation beyond which localizations are no longer clustered — roughly the
    # domain size. The g(r) first-crossing estimate is fragile (noisy finite
    # difference, can latch onto an early dip); pass domain_scale to set it directly
    # from a known/robustly-measured domain radius.
    n_domains_g, kappa, r_domain = estimate_n_domains_from_g(pcf_real, r_domain=domain_scale)

    # Within-domain RDF: keep only bins inside the domain scale. Beyond r_domain the
    # neighbour signal is inter-domain arrangement, not the domain's own profile;
    # leaving it in lets the ring-area weighting in spawn_nanodomains inflate each
    # domain to hundreds of nm (domains then merge and the noise estimate collapses).
    domain_bins = max(1, int(round(r_domain / step)))
    rdf_dom = rdf.copy()
    rdf_dom[domain_bins:] = 0.0

    # Localizations per domain from the within-domain RDF integral. Domain count =
    # total domain localizations / per-domain localizations; both scale with
    # labelling density so the count is density-independent (no DBSCAN).
    n_locs = max(extract_n_locs_from_rdf(rdf_dom, step=step), 1.0)

    # Reference: a noise-free domain field built from ALL localizations (domain
    # fraction 1). Its pair correlation is the pure-domain g(r) used to read the
    # noise fraction off the real data via the superposition scaling.
    n_ref = max(1, int(round(len(df) / n_locs)))
    ref_seeds = place_seeds_matched(contour_bands, profile, df["z [nm]"].values,
                                    x_min, y_min, n_seeds=n_ref, px_size=px_size, rng=rng)
    ref_domains = spawn_nanodomains(ref_seeds, rdf=rdf_dom, adf=adf, n_locs=n_locs, step=step)
    pcf_ref = measure_pair_correlation(ref_domains, mask, x_min, y_min,
                                       r_max=sdis, dr=2 * step, px_size=px_size)

    # Noise fraction: g(r)-based (density-independent) by default; the local-density
    # heuristic is kept for comparison. An explicit noise_fraction overrides both.
    f_noise_g = estimate_noise_fraction_from_pcf(pcf_real, pcf_ref, r_fit_max=r_domain)
    f_noise_heuristic = estimate_noise_fraction(df, contour_bands, x_min, y_min,
                                                rdf=rdf, step=step, bin_size=px_size, verbose=False)
    f_noise = noise_fraction if noise_fraction is not None else f_noise_g

    n_domains = max(1, int(round((1.0 - f_noise) * len(df) / n_locs)))

    # Read a hard-core spacing from a dispersion trough in L(r) - r beyond the
    # domain scale (L - r < 0 means regularity). Used only in inhibition mode.
    min_dist = None
    if match_spacing:
        beyond = pcf_real["r"] > r_domain
        Lr = pcf_real["L"] - pcf_real["r"]
        if np.any(beyond) and np.nanmin(Lr[beyond]) < 0:
            min_dist = float(pcf_real["r"][beyond][np.nanargmin(Lr[beyond])])

    # Seed placement.
    #   "thomas":     Neyman-Scott cluster process fit to the inter-domain g(r),
    #                 reproducing domain-of-domains clustering (falls back to Poisson
    #                 if no clustered signal is found).
    #   "inhibition": radial placement with a hard-core spacing from L(r) - r.
    #   "poisson":    independent radial placement.
    kappa_p = sigma = None
    placement_used = seed_placement
    if seed_placement == "thomas":
        kappa_p, sigma, ok = fit_thomas_parameters(pcf_real, r_domain, f_noise, r_max_fit=sdis)
        if ok:
            seeds = place_seeds_thomas(contour_bands, profile, df["z [nm]"].values,
                                       x_min, y_min, n_seeds=n_domains,
                                       kappa_p=kappa_p, sigma=sigma, px_size=px_size, rng=rng)
        else:
            placement_used = "poisson (thomas fit failed)"
            seeds = place_seeds_matched(contour_bands, profile, df["z [nm]"].values,
                                        x_min, y_min, n_seeds=n_domains, px_size=px_size, rng=rng)
    else:
        seeds = place_seeds_matched(contour_bands, profile, df["z [nm]"].values,
                                    x_min, y_min, n_seeds=n_domains,
                                    min_dist=(min_dist if seed_placement == "inhibition" else None),
                                    px_size=px_size, rng=rng)
    if len(seeds) == 0:
        raise ValueError("Seed placement produced no seeds; check the radial profile.")

    domains_df = spawn_nanodomains(seeds, rdf=rdf_dom, adf=adf, n_locs=n_locs, step=step)

    if noise_mode == "uniform":
        noise_weights = None
    elif noise_mode == "profile":
        noise_weights = profile
    else:  # "support"
        noise_weights = (profile > 0).astype(float)

    n_noise = int(f_noise * len(df))
    noise_df = add_noise_locs(n_noise, contour_bands, df["z [nm]"].values,
                              x_min, y_min, px_size=px_size, band_weights=noise_weights)
    noisy_df = pd.concat([domains_df, noise_df], ignore_index=True)

    # Validate: pair correlation of the domain-only and the full (noisy) simulation.
    # pcf_noisy vs pcf_real is the apples-to-apples comparison (both include noise).
    pcf_domains = measure_pair_correlation(domains_df, mask, x_min, y_min,
                                           r_max=sdis, dr=2 * step, px_size=px_size)
    pcf_noisy = measure_pair_correlation(noisy_df, mask, x_min, y_min,
                                         r_max=sdis, dr=2 * step, px_size=px_size)

    info = {
        "n_domains": n_domains,
        "n_domains_g_estimate": n_domains_g,
        "n_seeds_placed": len(seeds),
        "kappa": kappa,
        "r_domain_nm": r_domain,
        "min_dist_nm": min_dist,
        "placement": placement_used,
        "thomas_kappa_p": kappa_p,
        "thomas_sigma_nm": sigma,
        "f_noise": f_noise,
        "f_noise_from_g": f_noise_g,
        "f_noise_heuristic": f_noise_heuristic,
        "n_locs_per_domain": n_locs,
        "seeds": seeds,
        "n_real": len(df),
        "n_domain_locs": len(domains_df),
        "n_noise": n_noise,
        "pcf_real": pcf_real,
        "pcf_domains": pcf_domains,
        "pcf_noisy": pcf_noisy,
        "pcf_reference": pcf_ref,
    }

    if verbose:
        src = "override" if noise_fraction is not None else "g(r)"
        thomas = f" | thomas sigma={sigma:.0f} nm, parents~{int(round(kappa_p*float(mask.sum())*px_size**2))}" if sigma else ""
        print(f"  domains (locs/domain ratio): {n_domains}  | g(r) cross-check: {n_domains_g}  | domain scale ~{r_domain:.0f} nm")
        print(f"  seed placement: {placement_used}  | seeds placed: {len(seeds)}{thomas}")
        print(f"  noise fraction: {f_noise:.2f} ({src})  | g(r): {f_noise_g:.2f}  heuristic: {f_noise_heuristic:.2f}")
        print(f"  real locs: {len(df):,} | domain-only: {len(domains_df):,} | with noise: {len(noisy_df):,}")

    return domains_df, noisy_df, info


def sample_conditional_on_z(query_z, reference_df, column, k=10, mode="sample", rng=None):
    """
    Empirically samples `column` from reference_df conditioned on z, via a
    k-nearest-neighbor bootstrap in z [nm]. Lets simulated values inherit
    whatever z-dependence exists in the real data (e.g. uncertainty widening
    away from the focal plane) without fitting a parametric model.

    Args:
        query_z:       1D array of z [nm] values to condition on, one per simulated loc
        reference_df:  real localization DataFrame carrying "z [nm]" and `column`
                       (pass the channel-matched df, e.g. me3_raw_df for me3 sim rows)
        column:        column to sample, e.g. "uncertainty [nm]" or "centroid [nm]"
        k:             number of nearest real neighbors (in z) to draw from
        mode:          "sample" draws a random neighbor's value per query (adds
                       realistic scatter); "median" returns the neighbors' median
                       (a single stable value per query point)
        rng:           np.random.Generator, required when mode="sample"

    Returns:
        1D np.ndarray, len(query_z), of sampled/aggregated values
    """
    query_z = np.asarray(query_z, dtype=float)
    ref_z = reference_df["z [nm]"].to_numpy(dtype=float)
    ref_vals = reference_df[column].to_numpy(dtype=float)

    k = min(k, len(ref_z))
    tree = cKDTree(ref_z.reshape(-1, 1))
    _, idx = tree.query(query_z.reshape(-1, 1), k=k)
    idx = idx.reshape(len(query_z), k)

    neighbor_vals = ref_vals[idx]  # (n_query, k)

    if mode == "median":
        return np.median(neighbor_vals, axis=1)
    elif mode == "sample":
        if rng is None:
            raise ValueError("rng is required when mode='sample'")
        choice = rng.integers(0, k, size=len(query_z))
        return neighbor_vals[np.arange(len(query_z)), choice]
    else:
        raise ValueError(f"unknown mode: {mode!r}")


def export_to_thunderstorm(sim_df, reference_df, filename, centroid=None, frame=1, k=10, mode="sample", rng=None, seed=None):
    """
    Formats a simulated localization DataFrame and writes it to a ThunderSTORM-
    importable CSV (for ImageJ).

    ThunderSTORM's importer needs at minimum id/frame/x/y[/z], and its default
    Gaussian rendering needs "uncertainty [nm]" for weighting, so uncertainty
    is sampled from `reference_df` via sample_conditional_on_z, letting
    simulated rows inherit the real channel's z-dependence instead of getting
    a flat constant. centroid [nm] (this pipeline's z-calibration field) is
    instead set to a single fixed value per channel, since me3/ac each have
    one already-established centroid rather than a per-row distribution.
    Pass the reference_df matching the channel the sim came from (me3_raw_df
    for a me3-derived sim_df, ac_raw_df for ac) — run this separately per
    channel rather than mixing channels into one call.

    Args:
        sim_df:        simulated DataFrame with [x [nm], y [nm], z [nm], cluster_label]
        reference_df:  real raw localization DataFrame for the SAME channel
                       (must carry "z [nm]", "uncertainty [nm]"; also
                       "centroid [nm]" if centroid is left as None)
        filename:      CSV path to write the result to
        centroid:      single centroid [nm] value to assign to every row of
                       this channel. None defaults to reference_df["centroid [nm]"].median()
        frame:         constant frame index (this pipeline has no time axis)
        k, mode:       passed through to sample_conditional_on_z (uncertainty only)
        rng:           np.random.Generator for the uncertainty bootstrap. None
                       creates one internally (seeded with `seed` if given)
        seed:          seed for the auto-created rng when rng is None; ignored
                       if rng is passed explicitly

    Returns:
        pd.DataFrame with ThunderSTORM-style columns, matching what was written
    """
    z = sim_df["z [nm]"].to_numpy(dtype=float)

    if rng is None:
        rng = np.random.default_rng(seed)

    if centroid is None:
        centroid = reference_df["centroid [nm]"].median()

    out = pd.DataFrame({
        "id": np.arange(1, len(sim_df) + 1, dtype=float),
        "frame": np.full(len(sim_df), frame, dtype=float),
        "x [nm]": sim_df["x [nm]"].to_numpy(dtype=float),
        "y [nm]": sim_df["y [nm]"].to_numpy(dtype=float),
        "z [nm]": z,
        "uncertainty [nm]": sample_conditional_on_z(z, reference_df, "uncertainty [nm]", k=k, mode=mode, rng=rng),
        "centroid [nm]": np.full(len(sim_df), centroid, dtype=float),
        "cluster_label": sim_df["cluster_label"].to_numpy(),
    })

    out.to_csv(filename, index=False)

    return out