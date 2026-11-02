"""
Synthetic single-channel nuclei with known (ground-truth) nanodomains, for
calibrating and validating DBSCAN (see calibrate).

The main entry point is generate_nucleus, which rebuilds a real cell from the
statistics in `measure`: seeds (domain centres) follow the real radial density
profile and inter-domain g(r); localizations are spawned around each seed from
the real within-domain RDF / ADF; uniform background noise makes up the
estimated noise fraction. The building blocks are exposed for custom scenes:

    seed z sampling   build_z_band_pools
    seed placement    place_seeds (radial Poisson), place_seeds_matched
                      (radial + spacing inhibition), place_seeds_thomas
                      (radial + Thomas cluster process), make_grid_seeds
    localizations     spawn_nanodomains, add_noise_locs, add_sim_noise
    scenes            generate_nucleus, simulate_single_domain_scene

Simulated DataFrames carry "x [nm]", "y [nm]", "z [nm]" and a ground-truth
"cluster_label" (domain id; -1 = background noise).

Randomness: generate_nucleus / place_seeds_matched / place_seeds_thomas take an
`rng` (np.random.Generator), but spawn_nanodomains, add_noise_locs,
add_sim_noise, place_seeds and make_grid_seeds draw from the global np.random
state, so seed np.random as well when you need a reproducible nucleus.
"""
import numpy as np
import pandas as pd

from .measure import (calculate_geometric_gt_target, estimate_n_domains_from_g, estimate_noise_fraction,
                      estimate_noise_fraction_from_pcf, extract_empirical_parameters, extract_n_locs_from_rdf,
                      extract_radial_density_profile, fit_thomas_parameters, measure_pair_correlation)

_default_bins = np.arange(100)
# rdf is a per-ring density (locs/nm^2); spawn_nanodomains weights it by ring area.
_default_rdf = np.exp(-0.5 * (_default_bins / 15.0) ** 2)
_default_adf = np.exp(-0.5 * (_default_bins / 15.0) ** 2)


# --- Seed z sampling -------------------------------------------------------------

def build_z_band_pools(df, contour_bands, x_min, y_min, bin_size=50.0, min_band_n=20):
    """
    Groups real z [nm] values by the XY contour_band their (x, y) falls in, so
    seed placement can bootstrap z conditioned on radial position instead of
    from the nucleus-wide marginal.

    place_seeds / place_seeds_matched / place_seeds_thomas previously drew every
    seed's z from the whole-nucleus pool regardless of where in xy it landed —
    seeds spanned the real axial range, but a seed's z carried no information
    about its xy position, discarding whatever xy/z correlation the real
    nucleus has (e.g. the nucleus is thinner near its rim than at its centre).
    Passing this function's output to those functions instead ties z to the
    seed's own band.

    Bands with fewer than min_band_n real points fall back to the nucleus-wide
    pool (returned as global_z), since bootstrapping from only a handful of
    points gives an unstably peaked z distribution for that band — this mainly
    matters for sparse interior bands (e.g. inside a nucleolus).

    Args:
        df:            DataFrame with [x [nm], y [nm], z [nm]] from a single channel
        contour_bands: 2D int array from create_radial_contours (-1 = outside nucleus)
        x_min, y_min:  mask origin in nm used to build contour_bands (the SAME
                       origin passed to extract_radial_density_profile / generate_nucleus
                       — not necessarily df's own min x/y; see preprocess.mask_origin)
        bin_size:      pixel size in nm (must match contour_bands)
        min_band_n:    minimum real points required to use a band's own z pool

    Returns:
        (z_by_band, global_z):
            z_by_band: dict {band: np.ndarray of z [nm]}, only for bands with
                       >= min_band_n real points
            global_z:  full-nucleus z [nm] array, the fallback pool for every
                       other band
    """
    coords = df[["x [nm]", "y [nm]"]].to_numpy()
    x_idx = np.clip(((coords[:, 0] - x_min) / bin_size).astype(int), 0, contour_bands.shape[0] - 1)
    y_idx = np.clip(((coords[:, 1] - y_min) / bin_size).astype(int), 0, contour_bands.shape[1] - 1)
    bands = contour_bands[x_idx, y_idx]

    global_z = df["z [nm]"].to_numpy(dtype=float)

    z_by_band = {}
    for b in np.unique(bands[bands >= 0]):
        z_vals = global_z[bands == b]
        if len(z_vals) >= min_band_n:
            z_by_band[int(b)] = z_vals

    return z_by_band, global_z


def _sample_banded_z(bands, z_by_band, global_z, rng=np.random):
    """
    Samples one z value per entry in `bands`, bootstrapped from the real z pool
    of that entry's own contour_band (z_by_band), falling back to global_z for
    bands without a dedicated pool (see build_z_band_pools).

    `rng` may be the np.random module itself or an np.random.Generator instance
    — both expose a compatible .choice(pool, size=n), so callers that use the
    legacy module-level np.random state and callers that thread through a
    Generator can share this helper.
    """
    bands = np.asarray(bands)
    z = np.empty(len(bands), dtype=float)
    for b in np.unique(bands):
        pool = z_by_band.get(int(b), global_z)
        sel = bands == b
        z[sel] = rng.choice(pool, size=int(sel.sum()))
    return z


# --- Seed placement --------------------------------------------------------------

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


def place_seeds(contour_bands, band_density_profile, z_by_band, global_z, x_min, y_min, px_size=50.0, scaling_factor=1.0):
    """
    Places nanodomain seeds within the nucleus using a radially-weighted Poisson process.

    Each pixel inside the nucleus draws its seed count independently from
    Poisson(density * scaling_factor), where density comes from band_density_profile.
    This correctly allows high-density pixels to host multiple seeds and zero-density
    pixels to host none. Sub-pixel jitter prevents seeds from landing on a pixel-center grid.

    Z-coordinates are bootstrapped from the real localizations in that seed's own
    contour_band (z_by_band, from build_z_band_pools), so a seed's axial position is
    tied to its radial position instead of drawn from the whole-nucleus marginal.
    Bands without enough real points to bootstrap from fall back to global_z.

    Args:
        contour_bands:        2D int array from create_radial_contours (-1 = outside nucleus)
        band_density_profile: 1D float array from extract_radial_density_profile
        z_by_band, global_z:  from build_z_band_pools(df, contour_bands, x_min, y_min, ...)
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

    bands = contour_bands[seed_indices[:, 0], seed_indices[:, 1]]
    z = _sample_banded_z(bands, z_by_band, global_z, rng=np.random)

    return np.column_stack([xy, z])


def place_seeds_matched(contour_bands, band_density_profile, z_by_band, global_z,
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

    Z-coordinates are bootstrapped per-seed from that seed's own contour_band
    (z_by_band, from build_z_band_pools), matching place_seeds.

    Args:
        contour_bands:        2D int array from create_radial_contours (-1 outside)
        band_density_profile: 1D float array from extract_radial_density_profile
        z_by_band, global_z:  from build_z_band_pools(df, contour_bands, x_min, y_min, ...)
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
                band = contour_bands[rows[i], cols[i]]
                pool = z_by_band.get(int(band), global_z)
                placed_z[n_placed] = pool[rng.integers(len(pool))]
                n_placed += 1
                if n_placed >= n_seeds:
                    break

    return np.column_stack([placed_xy[:n_placed], placed_z[:n_placed]])


def place_seeds_thomas(contour_bands, band_density_profile, z_by_band, global_z, x_min, y_min,
                       n_seeds, kappa_p, sigma, px_size=50.0, rng=None):
    """
    Places nanodomain seeds as a Thomas (Neyman-Scott) cluster process, reproducing
    inter-domain clustering while still following the radial density profile.

    Parents are placed with the radial weighting (place_seeds_matched); each parent
    scatters offspring domain centres with a Gaussian(sigma) displacement in xy.
    Offspring falling outside the nucleus mask are discarded, so the returned count
    can be slightly below n_seeds. z is bootstrapped per offspring from that
    offspring's own contour_band (z_by_band, from build_z_band_pools) — not
    inherited from the parent, since the Gaussian xy jitter can move an offspring
    into a different band than its parent.

    Fit kappa_p and sigma with fit_thomas_parameters. Use this in place of
    place_seeds_matched when the real cell shows inter-domain clustering (elevated
    L(r) - r beyond the domain scale).

    Args:
        contour_bands:        2D int array from create_radial_contours (-1 outside)
        band_density_profile: 1D float array from extract_radial_density_profile
        z_by_band, global_z:  from build_z_band_pools(df, contour_bands, x_min, y_min, ...)
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
    mask = contour_bands >= 0
    area = float(mask.sum()) * px_size ** 2

    n_parents = max(1, int(round(kappa_p * area)))
    parents = place_seeds_matched(contour_bands, band_density_profile, z_by_band, global_z,
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

    bands = contour_bands[xi[inside], yi[inside]]
    z = _sample_banded_z(bands, z_by_band, global_z, rng=rng)
    return np.column_stack([xy, z])


# --- Nanodomains and background --------------------------------------------------

def spawn_nanodomains(seeds, rdf=_default_rdf, adf=_default_adf, n_locs=50, step=10, start_label=0):
    """
    Spawns localizations around each seed using Poisson-sampled counts.

    n_locs is the mean number of localizations per seed — each seed draws its
    actual count from Poisson(n_locs), giving realistic domain-to-domain
    variability. Use measure.extract_n_locs_from_rdf to derive n_locs from real
    data. Draws from the global np.random state.

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


# --- Whole-nucleus and single-domain scenes --------------------------------------

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

    Localizations per domain (n_locs) come from the within-domain RDF integral
    (extract_n_locs_from_rdf) and the domain count is (1 - f_noise) * len(df) /
    n_locs, so the domain-only total is (1 - f_noise) of the real count and noise
    adds the remaining f_noise; both counts scale with labelling density, so the
    domain count does not (no DBSCAN involved). estimate_n_domains_from_g, a
    Poisson-cluster fit to g(r), supplies the domain scale (unless domain_scale
    is given) and an independent count cross-check (info["n_domains_g_estimate"]).
    Domain centres are arranged according to seed_placement (default: a Thomas
    cluster process fit to the inter-domain g(r)).

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
        rng:           optional np.random.Generator for seed placement. Not fully
                       reproducible on its own: spawn_nanodomains and
                       add_noise_locs draw from the global np.random state, so
                       also call np.random.seed(...) for a repeatable nucleus.
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

    # Bootstrap seed z per-band rather than from the whole-nucleus marginal, so a
    # seed's axial position is tied to its radial position (see build_z_band_pools).
    z_by_band, global_z = build_z_band_pools(df, contour_bands, x_min, y_min, bin_size=px_size)

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
    ref_seeds = place_seeds_matched(contour_bands, profile, z_by_band, global_z,
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
            seeds = place_seeds_thomas(contour_bands, profile, z_by_band, global_z,
                                       x_min, y_min, n_seeds=n_domains,
                                       kappa_p=kappa_p, sigma=sigma, px_size=px_size, rng=rng)
        else:
            placement_used = "poisson (thomas fit failed)"
            seeds = place_seeds_matched(contour_bands, profile, z_by_band, global_z,
                                        x_min, y_min, n_seeds=n_domains, px_size=px_size, rng=rng)
    else:
        seeds = place_seeds_matched(contour_bands, profile, z_by_band, global_z,
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


def simulate_single_domain_scene(domain_diameter_nm=150.0, n_domain_locs=300,
                                 noise_density=5e-6, scene_size_nm=2000.0,
                                 z_extent_nm=200.0, rng=None):
    """
    Generates ONE synthetic nanodomain (a uniform-density 3D ball of
    localizations) embedded in uniform background noise across a square scene.

    generate_nucleus's dense field of many, possibly-overlapping domains has no
    unambiguous per-domain ground truth once domains merge (see
    optimize_dbscan_size_matching's docstring). A single isolated domain sidesteps
    that: there is exactly one true cluster, so its size is known exactly and any
    DBSCAN result can be checked against it directly, rather than only against a
    population-mean target.

    The domain's true diameter is measured with the same statistic
    calculate_geometric_gt_target / domain_size_cost use elsewhere (2x the 75th
    percentile of xy point-to-centroid distance), so it's directly comparable to
    whatever DBSCAN's predicted cluster sizes turn out to be, regardless of
    exactly how domain_diameter_nm maps to the sampled point cloud.

    Args:
        domain_diameter_nm: approximate diameter (nm) of the domain's sampling
                            ball; the actual measured true diameter (returned)
                            will be close to but not exactly this value
        n_domain_locs:      number of localizations in the domain at full
                            (density_fraction=1.0) density
        noise_density:      background noise localization density (locs/nm²) at
                            full density, applied over the whole scene area
        scene_size_nm:      side length (nm) of the square scene; the domain is
                            centered in it
        z_extent_nm:        full z-thickness (nm) of the scene; noise fills this
                            range uniformly, matching the domain's own z-extent
        rng:                np.random.Generator; a fresh default_rng() if None

    Returns:
        tuple: (scene_df, true_diameter_nm, scene_area_nm2)
            scene_df:         DataFrame [x [nm], y [nm], z [nm], true_label]
                              (true_label 1 = domain, 0 = background noise)
            true_diameter_nm: the domain's measured geometric_2D diameter (nm)
            scene_area_nm2:   scene_size_nm ** 2, for computing localization
                              density consistently with the calibration
    """
    if rng is None:
        rng = np.random.default_rng()

    # --- domain: uniform points within a 3D ball (r ~ U^(1/3) for uniform volume density) ---
    domain_radius = domain_diameter_nm / 2.0
    directions = rng.normal(size=(n_domain_locs, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    radii = domain_radius * rng.random(n_domain_locs) ** (1.0 / 3.0)
    domain_offsets = directions * radii[:, None]

    center = np.array([scene_size_nm / 2.0, scene_size_nm / 2.0, z_extent_nm / 2.0])
    domain_xyz = center + domain_offsets

    # --- background noise: uniform over the whole scene ---
    n_noise_locs = int(round(noise_density * scene_size_nm ** 2))
    noise_xyz = rng.random((n_noise_locs, 3)) * [scene_size_nm, scene_size_nm, z_extent_nm]

    scene_df = pd.DataFrame(
        np.vstack([domain_xyz, noise_xyz]),
        columns=["x [nm]", "y [nm]", "z [nm]"]
    )
    scene_df["true_label"] = np.concatenate([
        np.ones(n_domain_locs, dtype=int),
        np.zeros(n_noise_locs, dtype=int),
    ])

    domain_only = scene_df[scene_df["true_label"] == 1].copy()
    domain_only["cluster_label"] = 0
    true_diameter_nm = calculate_geometric_gt_target(domain_only, min_points=1, percentile=75.0)

    scene_area_nm2 = scene_size_nm ** 2
    return scene_df, true_diameter_nm, scene_area_nm2
