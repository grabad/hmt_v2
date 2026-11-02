"""
Statistics measured from a single channel of real (or simulated) localizations,
used both to parametrise the nucleus simulation (simulate.generate_nucleus) and
to check that a simulation reproduces the real cell.

    radial profile     extract_radial_density_profile: density per contour band
    within-domain      extract_empirical_parameters (RDF / ADF),
                       extract_n_locs_from_rdf, estimate_noise_fraction
    pair correlation   measure_pair_correlation (Ripley's K / L and g(r), edge-
                       corrected against the nucleus mask), plus the g(r)-based
                       estimators estimate_n_domains_from_g,
                       estimate_noise_fraction_from_pcf, fit_thomas_parameters
    domain size        calculate_geometric_gt_target (2 x 75th-percentile
                       radius, the size statistic shared with calibrate / cluster)

The pair-correlation statistics are invariant to random under-labelling, so
they compare domain structure across cells of different staining density;
DBSCAN-based measurements (measure_seed_spacing) are not.

Functions that map points onto the mask take the mask origin (x_min, y_min)
from preprocess.mask_origin; see preprocess for the grid conventions.
"""
import numpy as np
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN


# --- Radial (nuclear-scale) distribution -----------------------------------------

def extract_radial_density_profile(df, contour_bands, x_min, y_min, bin_size=50):
    """
    Extracts localization density per contour band from experimental data.

    For each band, divides the number of localizations in that band by the
    number of pixels in that band (a proxy for band area). The result is the
    band_density_profile argument for place_seeds — it encodes where in the
    nucleus the mark is enriched (peripheral vs. interior).

    Band numbering from create_radial_contours: 0 = innermost core,
    num_bands - 1 = outermost periphery, -1 = outside nucleus.

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


# --- Within-domain point distribution --------------------------------------------

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


# --- Pair correlation / Ripley's K (density-independent) -------------------------

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
    correction), estimated by convolving the mask with disk kernels. Points
    outside the mask are left out (as centres, as neighbours and from lambda),
    so a simulation whose domains spill over the mask edge is measured the same
    way as masked real data.

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
            'lambda': intensity of the in-mask localizations (localizations / nm^2)
            'area':   nucleus area in nm^2
    """
    mask = np.asarray(binary_mask) > 0
    area = float(mask.sum()) * px_size ** 2

    coords = df[["x [nm]", "y [nm]"]].to_numpy()

    # Map each point to its mask pixel (same convention as extract_radial_density_profile)
    x_idx = np.clip(np.floor((coords[:, 0] - x_min) / px_size).astype(int), 0, mask.shape[0] - 1)
    y_idx = np.clip(np.floor((coords[:, 1] - y_min) / px_size).astype(int), 0, mask.shape[1] - 1)

    # K is estimated inside the mask only. A point outside it (e.g. a simulated domain spilling over the
    # nuclear edge) has almost none of its search disk inside, so the edge correction below would weight
    # its neighbour counts by up to 1 / 1e-3 and put spikes into g(r).
    inside = mask[x_idx, y_idx]
    coords, x_idx, y_idx = coords[inside], x_idx[inside], y_idx[inside]
    n = len(coords)
    lam = n / area

    radii = np.arange(dr, r_max + dr, dr)

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
        inside = fftconvolve(mask_f, disk, mode="same")  # disk pixels inside the mask
        # Normalise by the pixelated disk's own area, not pi r^2: at radii of a few pixels the two
        # differ a lot (5 px vs pi px^2 at r = px_size) and the ratio jumps from one radius to the
        # next, which made interior points' frac oscillate around 1 and put spikes and zeros into
        # g(r) whenever dr < px_size
        frac = np.clip(inside[x_idx, y_idx] / disk.sum(), 1e-3, 1.0)

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


# --- Domain size and spacing -----------------------------------------------------

def calculate_geometric_gt_target(df, min_points=10, percentile=75.0):
    """
    Computes a robust characteristic domain diameter from ground-truth
    cluster_label groups, for size-matching calibration (see domain_size_cost /
    optimize_dbscan_size_matching). Ported from this project's original
    calibration approach, which predates the ARI-based search.

    For each true domain (excluding noise and domains with fewer than
    min_points), computes 2x the `percentile`-th percentile of point-to-
    centroid distance in 2D (xy only) — a robust stand-in for diameter that
    rejects outlier points without needing a convex hull, and stays stable on
    the sparse, jagged point clouds seen at low labelling-efficiency density
    fractions. domain_size_cost measures predicted DBSCAN clusters with the
    exact same statistic, so the two are directly comparable.

    Args:
        df:          DataFrame with [x [nm], y [nm], cluster_label] (ground
                     truth; noise rows with cluster_label == -1 are excluded)
        min_points:  minimum points a true domain needs to contribute a size;
                     domains with fewer points are skipped, not penalised —
                     unlike domain_size_cost's treatment of small PREDICTED
                     clusters (see that function's docstring for why the two
                     are handled asymmetrically)
        percentile:  percentile of point-to-centroid distance defining the
                     domain radius (75 matches the original calibration); must
                     match the percentile passed to domain_size_cost

    Returns:
        float: mean domain diameter (nm) across all qualifying true domains
    """
    domains_only = df[df["cluster_label"] != -1]
    sizes = []
    for _, group in domains_only.groupby("cluster_label"):
        if len(group) < min_points:
            continue
        coords = group[["x [nm]", "y [nm]"]].to_numpy()
        centroid = coords.mean(axis=0)
        distances = np.linalg.norm(coords - centroid, axis=1)
        sizes.append(2.0 * np.percentile(distances, percentile))

    if not sizes:
        raise ValueError("No true domain had >= min_points points; cannot compute a target diameter.")
    return float(np.mean(sizes))


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
