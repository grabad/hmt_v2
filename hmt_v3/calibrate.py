"""
DBSCAN parameter calibration against simulated ground truth (manuscript Fig. 4).

A fixed DBSCAN (eps, min_samples) over- or under-merges nanodomains as the
labelling density changes, which biases measured domain sizes. The calibration
simulates a nucleus matching a real cell (simulate.generate_nucleus), thins it
to a ladder of densities (preprocess.sample_lower_densities), finds the best
(eps, min_samples) at each density, and fits

    eps(density)         = eps_coef * density ** eps_exp
    min_samples(density) = min_samples_slope * density + min_samples_intercept

which cluster.cluster_dbscan_density_corrected then evaluates per cell.

    size matching (used by the drivers)
        optimize_dbscan_size_matching, minimising domain_size_cost: mean
        predicted domain diameter vs. the ground-truth diameter
    ground-truth diagnostics
        diagnose_domain_packing (do simulated domains physically overlap?),
        sweep_domain_scale (pick generate_nucleus's domain_scale)
    ARI / F1 objective (earlier approach, kept for comparison)
        optimize_epsilon, optimize_dbscan_params, epsilon_cost_cluster,
        domain_detection_f1
"""
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics.cluster import contingency_matrix

from .measure import calculate_geometric_gt_target
from .simulate import generate_nucleus


# --- Size-matching calibration (used by the drivers) -----------------------------

def domain_size_cost(eps, coords, gt_target, min_samples=8, min_points=10,
                     percentile=75.0, continent_multiple=4.0, use_z=False):
    """
    Cost function for size-matching calibration: |mean predicted-cluster
    diameter - gt_target|, in nm. Ported from this project's original
    calibration approach (predates the ARI-based search in optimize_epsilon /
    optimize_dbscan_params).

    Unlike epsilon_cost_cluster (1 - ARI), this never compares point-level
    labels — it only asks whether DBSCAN's clusters are, ON AVERAGE, the right
    physical size. That is a much smoother, less brittle search signal: a
    handful of genuinely hard-to-recover domains barely moves a mean diameter,
    where ARI's point-pair sum and domain_detection_f1's exact IoU threshold
    can collapse from the same handful of errors. The tradeoff is that it is a
    coarser criterion — it cannot distinguish "right average size, wrong
    domain-to-domain correspondence" from a genuinely correct recovery. Treat
    it as a search objective, not a substitute for validating a converged
    calibration with domain_detection_f1 / ARI.

    Cluster diameter uses the same 2x-percentile-radius statistic as
    calculate_geometric_gt_target, so the two are directly comparable.
    Predicted clusters with fewer than min_points are still counted, but with
    a small fixed placeholder size rather than skipped (matching the original
    implementation) — this penalises DBSCAN fragmenting domains into tiny,
    noise-like litter by pulling the mean size down, rather than ignoring the
    fragments outright. Clusters whose diameter exceeds
    continent_multiple * gt_target are flagged as merged domains ("continents")
    and given a large fixed penalty scaled to gt_target — this is what makes
    the search actively avoid over-merging, rather than just averaging over it
    (a few oversized "continents" would otherwise barely nudge a large mean).

    Args:
        eps:                DBSCAN neighbourhood radius in nm
        coords:             (N, 3) array of [x, y, z] in nm
        gt_target:          target mean domain diameter in nm, from
                            calculate_geometric_gt_target
        min_samples:        DBSCAN min_samples parameter
        min_points:         minimum points a predicted cluster needs before
                            its actual size (rather than the placeholder) is
                            measured; must match calculate_geometric_gt_target
        percentile:         percentile of point-to-centroid distance defining
                            cluster radius; must match calculate_geometric_gt_target
        continent_multiple: a cluster with diameter > continent_multiple *
                            gt_target is treated as multiple merged domains
        use_z:              if True, cluster in 3D (XYZ); the diameter itself
                            is still always measured in 2D (xy), matching the
                            original implementation

    Returns:
        float >= 0: |mean predicted diameter - gt_target| in nm (lower is
        better); a large fixed value if DBSCAN found no valid clusters
    """
    fit_coords = coords[:, :3] if use_z else coords[:, :2]
    labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(fit_coords)

    valid = labels != -1
    if not np.any(valid):
        return 1e5

    xy = coords[valid, :2]
    valid_labels = labels[valid]
    continent_cutoff = continent_multiple * gt_target
    continent_penalty = continent_cutoff * 5.0

    sizes = []
    for cluster_id in np.unique(valid_labels):
        cluster_xy = xy[valid_labels == cluster_id]
        if len(cluster_xy) < min_points:
            sizes.append(20.0)  # matches the original implementation's fixed placeholder
            continue
        centroid = cluster_xy.mean(axis=0)
        distances = np.linalg.norm(cluster_xy - centroid, axis=1)
        diameter = 2.0 * np.percentile(distances, percentile)
        sizes.append(continent_penalty if diameter > continent_cutoff else diameter)

    if not sizes:
        return 1e5

    return abs(float(np.mean(sizes)) - gt_target)


def optimize_dbscan_size_matching(fraction_dfs, calibration_area_nm2, gt_target=None,
                                  min_samples_range=range(3, 21), continent_multiple=4.0,
                                  min_points=10, percentile=75.0, cost_threshold=50.0,
                                  use_z=False, show_plot=True, verbose=True):
    """
    Finds the DBSCAN (epsilon, min_samples) pair whose predicted clusters best
    match a target domain SIZE, for each density fraction, then fits power-law
    calibration curves for both parameters vs. density.

    This is a port of this project's original calibration approach
    (legacy/hmt_functions.py: epsilon_cost / optimize_epsilon_hybrid), extended to search
    min_samples jointly with epsilon (the original only optimised epsilon at a
    fixed min_samples). It minimises domain_size_cost: the absolute difference
    between DBSCAN's mean predicted-cluster diameter and a ground-truth target
    diameter (calculate_geometric_gt_target) — it never compares point-level
    cluster labels at all, purely whether DBSCAN's clusters are, on average,
    the right physical size. See domain_size_cost's docstring for the
    reasoning and its known limitation (right average size doesn't guarantee
    correct domain-to-domain correspondence — it can be fooled when true
    domains overlap heavily, since merging them then doesn't inflate the
    merged cluster's footprint much).

    gt_target is computed once — from fraction_dfs[0] via
    calculate_geometric_gt_target — and reused for every fraction, since
    domain diameter is a physical property of the simulation and shouldn't be
    re-derived from increasingly sparse subsampled data (sample_lower_densities
    already orders fractions from full density down, so fraction_dfs[0] is the
    full/un-subsampled fraction). Pass gt_target explicitly to override.

    Args:
        fraction_dfs:          list of DataFrames, one per density fraction,
                               each with [x [nm], y [nm], z [nm], cluster_label]
        calibration_area_nm2:  nucleus area in nm² used to compute density
        gt_target:             target mean domain diameter in nm; if None,
                               computed from fraction_dfs[0]
        min_samples_range:     iterable of candidate min_samples values to search
        continent_multiple:    forwarded to domain_size_cost
        min_points:            forwarded to domain_size_cost / calculate_geometric_gt_target
        percentile:            forwarded to domain_size_cost / calculate_geometric_gt_target
        cost_threshold:        a fraction is kept in the eps/min_samples calibration
                               fit only if its best size-match error is below
                               this many nm
        use_z:                 if True, cluster in 3D (XYZ); if False, 2D (XY only)
        show_plot:             if True, plots eps/min_samples calibration curves
                               plus the size-match-error diagnostic
        verbose:               if True (default), prints per-min_samples progress
                               within each fraction and a final summary

    Returns:
        dict with keys:
            'gt_target': the target diameter used (nm)
            'eps_coef', 'eps_exp', 'eps_r_squared':
                eps = eps_coef * density^eps_exp, and its log-log fit R²
            'min_samples_slope', 'min_samples_intercept', 'min_samples_r_squared':
                min_samples = slope * density + intercept, and its linear fit R²
            'densities', 'best_epsilons', 'best_min_samples', 'best_size_errors':
                values for fractions that passed cost_threshold (used in the fits)
            'all_densities', 'all_epsilons', 'all_min_samples', 'all_size_errors':
                the same, for EVERY fraction that produced a result, regardless
                of cost_threshold
    """
    if gt_target is None:
        gt_target = calculate_geometric_gt_target(fraction_dfs[0], min_points=min_points,
                                                   percentile=percentile)

    min_samples_list = list(min_samples_range)
    densities = []
    best_epsilons = []
    best_min_samples_list = []
    best_size_errors = []

    all_densities = []
    all_epsilons = []
    all_min_samples = []
    all_size_errors = []

    dims = "3D (XYZ)" if use_z else "2D (XY)"
    n_combos = len(fraction_dfs) * len(min_samples_list)
    if verbose:
        print(f"Starting joint (eps, min_samples) SIZE-MATCHING optimisation over {len(fraction_dfs)} "
              f"fractions [{dims}], target diameter {gt_target:.1f} nm, min_samples in "
              f"[{min_samples_list[0]}, {min_samples_list[-1]}] ({n_combos} combinations total)...")

    combo_count = 0
    for i, df in enumerate(fraction_dfs):
        coords = df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()

        if len(coords) < min(min_samples_list):
            if verbose:
                print(f"Fraction {i+1}/{len(fraction_dfs)}: Skipped (too few localisations)")
            combo_count += len(min_samples_list)
            continue

        density = len(coords) / calibration_area_nm2
        if verbose:
            print(f"Fraction {i+1}/{len(fraction_dfs)}: {len(coords):,} locs, "
                  f"density {density:.6e} locs/nm²")

        frac_best_eps = None
        frac_best_min_samples = None
        frac_best_error = np.inf
        frac_start = time.perf_counter()

        for min_samples in min_samples_list:
            combo_count += 1

            if len(coords) < min_samples:
                if verbose:
                    print(f"    [{combo_count}/{n_combos}] min_samples={min_samples:2d}: "
                          f"skipped (fewer locs than min_samples)")
                continue

            # --- Stage 1: coarse grid search over epsilon ---
            coarse_grid = np.arange(10, 225, 15.0)
            best_coarse_eps = coarse_grid[0]
            best_coarse_cost = np.inf

            for test_eps in coarse_grid:
                cost = domain_size_cost(test_eps, coords, gt_target, min_samples,
                                        min_points=min_points, percentile=percentile,
                                        continent_multiple=continent_multiple, use_z=use_z)
                if cost < best_coarse_cost:
                    best_coarse_cost = cost
                    best_coarse_eps = test_eps

            # --- Stage 2: bounded scalar minimisation within ±15 nm of coarse winner ---
            search_min = max(10.0, best_coarse_eps - 15.0)
            search_max = min(225.0, best_coarse_eps + 15.0)

            res = minimize_scalar(
                domain_size_cost,
                bounds=(search_min, search_max),
                args=(coords, gt_target, min_samples, min_points, percentile, continent_multiple, use_z),
                method='bounded',
                options={'xatol': 0.2},
            )

            if not res.success:
                if verbose:
                    print(f"    [{combo_count}/{n_combos}] min_samples={min_samples:2d}: optimisation failed")
                continue

            error = res.fun
            is_new_best = error < frac_best_error
            if verbose:
                marker = "  <- new best for this fraction" if is_new_best else ""
                print(f"    [{combo_count}/{n_combos}] min_samples={min_samples:2d}: "
                      f"eps={res.x:6.2f} nm  size error={error:.1f} nm{marker}")
            if is_new_best:
                frac_best_error = error
                frac_best_eps = res.x
                frac_best_min_samples = min_samples

        elapsed = time.perf_counter() - frac_start

        if frac_best_eps is None:
            if verbose:
                print(f"Fraction {i+1}/{len(fraction_dfs)}: Density {density:.6e} | "
                      f"FAILED (no successful optimisation)  [{elapsed:.1f}s]")
            continue

        all_densities.append(density)
        all_epsilons.append(frac_best_eps)
        all_min_samples.append(frac_best_min_samples)
        all_size_errors.append(frac_best_error)

        passed = frac_best_error < cost_threshold
        if passed:
            densities.append(density)
            best_epsilons.append(frac_best_eps)
            best_min_samples_list.append(frac_best_min_samples)
            best_size_errors.append(frac_best_error)

        if verbose:
            status = "" if passed else f"  | above cost_threshold={cost_threshold:.0f}nm, excluded from calibration fit"
            print(f"Fraction {i+1}/{len(fraction_dfs)}: Density {density:.6e} -> "
                  f"Eps {frac_best_eps:.2f} nm, min_samples {frac_best_min_samples}  "
                  f"(size error {frac_best_error:.1f} nm){status}  [{elapsed:.1f}s]")

    result = {
        'gt_target': gt_target,
        'eps_coef': 0.0, 'eps_exp': 0.0, 'eps_r_squared': 0.0,
        'min_samples_slope': 0.0, 'min_samples_intercept': 0.0, 'min_samples_r_squared': 0.0,
        'densities': densities, 'best_epsilons': best_epsilons,
        'best_min_samples': best_min_samples_list, 'best_size_errors': best_size_errors,
        'all_densities': all_densities, 'all_epsilons': all_epsilons,
        'all_min_samples': all_min_samples, 'all_size_errors': all_size_errors,
    }

    if len(densities) < 2:
        if verbose:
            print(f"\nNot enough fractions passed cost_threshold={cost_threshold:.0f}nm for a calibration fit "
                  f"({len(densities)}/{len(all_densities)} passed) — no eps/min_samples calibration curve.")
        return result

    # --- Fit eps as a power law in log-log space, min_samples as a linear function of density ---
    log_dens = np.log10(densities)
    dens_arr = np.array(densities, dtype=float)

    log_eps = np.log10(best_epsilons)
    eps_par = np.polyfit(log_dens, log_eps, 1)
    eps_exp, eps_coef = eps_par[0], 10 ** eps_par[1]
    eps_r_squared = np.corrcoef(log_dens, log_eps)[0, 1] ** 2

    ms_arr = np.array(best_min_samples_list, dtype=float)
    ms_par = np.polyfit(dens_arr, ms_arr, 1)
    ms_slope, ms_intercept = ms_par[0], ms_par[1]
    ms_r_squared = np.corrcoef(dens_arr, ms_arr)[0, 1] ** 2

    result.update({
        'eps_coef': eps_coef, 'eps_exp': eps_exp, 'eps_r_squared': eps_r_squared,
        'min_samples_slope': ms_slope, 'min_samples_intercept': ms_intercept,
        'min_samples_r_squared': ms_r_squared,
    })

    if verbose:
        print(f"\nCalibration complete (target diameter {gt_target:.1f} nm):")
        print(f"  Eps         = {eps_coef:.2f} * Density^{eps_exp:.2f}  (R² = {eps_r_squared:.4f})")
        print(f"  min_samples = {ms_slope:.4g} * Density + {ms_intercept:.2f}  (R² = {ms_r_squared:.4f})")

    if show_plot:
        dens_smooth = np.linspace(dens_arr.min(), dens_arr.max(), 200)

        fig, axes = plt.subplots(2, 1, figsize=(7, 8))

        axes[0].scatter(dens_arr, best_epsilons, color='blue', label='Optimised epsilons')
        axes[0].plot(dens_smooth, eps_coef * dens_smooth ** eps_exp, color='red', linestyle='--',
                    label=f'Fit: eps = {eps_coef:.2f} · dens^{eps_exp:.2f}')
        axes[0].set_xlabel('Density [locs/nm²]')
        axes[0].set_ylabel('Epsilon [nm]')
        axes[0].set_title('Epsilon Calibration')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        axes[1].scatter(dens_arr, best_min_samples_list, color='blue', label='Optimised min_samples')
        axes[1].plot(dens_smooth, ms_slope * dens_smooth + ms_intercept, color='red', linestyle='--',
                    label=f'Fit: min_samples = {ms_slope:.4g} · dens + {ms_intercept:.2f}')
        axes[1].set_xlabel('Density [locs/nm²]')
        axes[1].set_ylabel('min_samples')
        axes[1].set_title('min_samples Calibration')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.suptitle('DBSCAN Calibration (size-matched, joint eps + min_samples)')
        plt.tight_layout()
        plt.show()

    return result


# --- Ground-truth diagnostics ----------------------------------------------------

def diagnose_domain_packing(df, domain_scale=None, use_z=True):
    """
    Checks whether true domains (ground-truth cluster_label) are packed close
    enough to geometrically overlap, independent of any clustering algorithm.

    For each true domain (cluster_label != -1), computes its centroid and a
    characteristic radius (the RMS distance of its own points from that
    centroid, unless domain_scale overrides it — see Args). It then finds each
    domain's nearest neighbouring domain by centroid distance and compares that
    distance against the sum of the two domains' radii.

    If a domain's nearest-neighbour distance is LESS than the sum of the two
    radii, their point clouds are expected to interpenetrate — no DBSCAN
    (eps, min_samples) choice, and no choice of use_z, can separate two domains
    whose localizations are genuinely interspersed in space, because that
    information isn't in the data. A high frac_overlapping here means an
    ARI/F1 ceiling is coming from how densely generate_nucleus's Thomas-process
    seed placement packs domains (fit to match the real inter-domain g(r)), not
    from the clustering search — run this before spending more time tuning
    DBSCAN on a fraction that shows a high overlap fraction.

    Args:
        df:           DataFrame with [x [nm], y [nm], z [nm], cluster_label]
                      (ground truth; noise rows with cluster_label == -1 are
                      excluded)
        domain_scale: if given, use this fixed radius (nm) for every domain
                      instead of each domain's own point scatter — pass the
                      r_domain that generate_nucleus reports so small/sparse
                      domains (few points, especially after subsampling to a
                      low density fraction) don't get an unreliable radius
                      estimated from just a handful of points
        use_z:        if True (default), compute centroids/radii/distances in
                      3D (XYZ); if False, 2D (XY only) — match whatever you're
                      about to cluster with

    Returns:
        dict with keys:
            'n_domains':        number of true domains found
            'centroids':        (n_domains, 2 or 3) array of domain centroids
            'radii':            (n_domains,) array of per-domain radii (nm)
            'nn_distance':      (n_domains,) nearest-neighbour centroid
                                 distance from each domain to its closest
                                 other domain (nm)
            'overlap_margin':   (n_domains,) nn_distance - (own radius +
                                 neighbour's radius); negative means that pair
                                 geometrically overlaps
            'frac_overlapping': fraction of domains whose nearest neighbour
                                 overlaps them
    """
    cols = ["x [nm]", "y [nm]", "z [nm]"] if use_z else ["x [nm]", "y [nm]"]
    domains_only = df[df["cluster_label"] != -1]

    centroids = []
    radii = []
    for _, sub in domains_only.groupby("cluster_label")[cols]:
        pts = sub.to_numpy()
        centroid = pts.mean(axis=0)
        centroids.append(centroid)
        if domain_scale is not None:
            radii.append(float(domain_scale))
        else:
            radii.append(float(np.sqrt(np.mean(np.sum((pts - centroid) ** 2, axis=1)))))

    centroids = np.asarray(centroids)
    radii = np.asarray(radii)
    n_domains = len(centroids)

    if n_domains < 2:
        return {
            'n_domains': n_domains, 'centroids': centroids, 'radii': radii,
            'nn_distance': np.array([]), 'overlap_margin': np.array([]),
            'frac_overlapping': 0.0,
        }

    tree = cKDTree(centroids)
    dist, idx = tree.query(centroids, k=2)  # column 0 is the self-match
    nn_distance = dist[:, 1]
    overlap_margin = nn_distance - (radii + radii[idx[:, 1]])
    frac_overlapping = float(np.mean(overlap_margin < 0))

    return {
        'n_domains': n_domains,
        'centroids': centroids,
        'radii': radii,
        'nn_distance': nn_distance,
        'overlap_margin': overlap_margin,
        'frac_overlapping': frac_overlapping,
    }


def sweep_domain_scale(df, contour_bands, x_min, y_min, domain_scales, crop_bounds,
                       noise_fraction=None, min_samples_range=range(5, 55, 10),
                       use_z=True, ari_threshold=0.3, iou_threshold=0.5,
                       generate_kwargs=None, rng_seed=0, show_plot=True, verbose=True):
    """
    Sweeps generate_nucleus's domain_scale to find a value where simulated
    domains are BOTH geometrically separable (diagnose_domain_packing) and
    density-resolvable by DBSCAN (domain_detection_f1) — rather than tuning
    either one in isolation, which is what led to the domain_scale=100 (90%
    geometric overlap) vs domain_scale=300 (0% domain recovery from density
    dilution/over-merging) whiplash: a value that fixes one can break the
    other, so both need to be checked together.

    For each candidate domain_scale: reruns generate_nucleus with a fixed rng
    seed (so only domain_scale/noise_fraction differ across candidates, not
    the random draw), crops the result to crop_bounds (matching your
    calibration crop), measures geometric overlap with diagnose_domain_packing,
    and runs ONE quick optimize_dbscan_params call on that single fraction
    (not a full density sweep) to get its ARI/F1. Checking only the crop's own
    density (no subsampling) is enough: if a domain_scale can't be separated by
    DBSCAN at full density, it won't be recoverable at any lower density either
    — subsampling only makes an already-broken layout sparser, not more
    separable (this is exactly why F1 came back 0 at every density fraction
    for domain_scale=300 rather than just the low-density ones).

    Args:
        df:                 masked single-channel DataFrame (e.g. me3_df), as
                            passed to generate_nucleus
        contour_bands, x_min, y_min: from create_radial_contours / mask_origin,
                            as passed to generate_nucleus
        domain_scales:      list of candidate domain_scale values (nm) to try
        crop_bounds:        (xmin, xmax, ymin, ymax) in nm — generate_nucleus's
                            noisy output is cropped to this window before
                            measuring packing/F1; use the same crop as your
                            calibration
        noise_fraction:     fixed noise_fraction passed to every generate_nucleus
                            call, or None (default) to let each domain_scale
                            re-derive its own g(r)-based estimate. None is
                            recommended: the g(r) noise estimate is itself
                            sensitive to domain_scale (see generate_nucleus),
                            so holding noise_fraction fixed while only varying
                            domain_scale silently conflates the two
        min_samples_range:  candidate min_samples values for the quick DBSCAN
                            check — kept coarse by default since this sweep is
                            meant to narrow down a region, not be the final
                            calibration; pass a finer range once you've picked
                            a domain_scale
        use_z:              cluster in 3D (XYZ) if True, else 2D (XY only)
        ari_threshold:      forwarded to optimize_dbscan_params (only affects
                            its internal power-law fit, which this function
                            doesn't use — irrelevant here, kept for parity)
        iou_threshold:      forwarded to domain_detection_f1
        generate_kwargs:    extra kwargs forwarded to generate_nucleus (e.g.
                            match_spacing, noise_mode); do not include
                            domain_scale, noise_fraction, rng, or verbose —
                            those are set by this function
        rng_seed:           fixed seed reused for every candidate domain_scale,
                            so differences across candidates reflect
                            domain_scale/noise_fraction, not a different
                            random draw
        show_plot:          if True, plots geometric overlap fraction and
                            ARI/F1 side by side vs domain_scale
        verbose:            if True, prints one summary line per candidate

    Returns:
        pd.DataFrame, one row per domain_scale, with columns:
            domain_scale, noise_fraction_used, n_domains_in_region,
            frac_overlapping, median_overlap_margin, eps, min_samples,
            ari, f1, precision, recall
    """
    generate_kwargs = generate_kwargs or {}
    xmin, xmax, ymin, ymax = crop_bounds
    calibration_area_nm2 = (xmax - xmin) * (ymax - ymin)

    records = []
    if verbose:
        print(f"Sweeping domain_scale over {list(domain_scales)}...")

    for domain_scale in domain_scales:
        rng = np.random.default_rng(rng_seed)
        np.random.seed(rng_seed)  # spawn_nanodomains / add_noise_locs draw from the global state, not rng
        _, noisy_df, info = generate_nucleus(
            df, contour_bands, x_min, y_min,
            domain_scale=domain_scale, noise_fraction=noise_fraction,
            rng=rng, verbose=False, **generate_kwargs,
        )

        cropped = noisy_df[noisy_df["x [nm]"].between(xmin, xmax) &
                           noisy_df["y [nm]"].between(ymin, ymax)].copy()

        packing = diagnose_domain_packing(cropped, domain_scale=domain_scale, use_z=use_z)

        dbscan_result = optimize_dbscan_params(
            [cropped], calibration_area_nm2, min_samples_range=min_samples_range,
            use_z=use_z, ari_threshold=ari_threshold, iou_threshold=iou_threshold,
            show_plot=False, verbose=False,
        )

        has_result = len(dbscan_result['all_aris']) > 0
        eps = dbscan_result['all_epsilons'][0] if has_result else np.nan
        n_min_samples = dbscan_result['all_min_samples'][0] if has_result else np.nan
        ari = dbscan_result['all_aris'][0] if has_result else np.nan
        f1 = dbscan_result['all_f1s'][0] if has_result else np.nan
        precision = dbscan_result['all_precisions'][0] if has_result else np.nan
        recall = dbscan_result['all_recalls'][0] if has_result else np.nan

        records.append({
            'domain_scale': domain_scale,
            'noise_fraction_used': info['f_noise'],
            'n_domains_in_region': packing['n_domains'],
            'frac_overlapping': packing['frac_overlapping'],
            'median_overlap_margin': (float(np.median(packing['overlap_margin']))
                                      if packing['n_domains'] >= 2 else np.nan),
            'eps': eps, 'min_samples': n_min_samples,
            'ari': ari, 'f1': f1, 'precision': precision, 'recall': recall,
        })

        if verbose:
            eps_str = f"{eps:.1f} nm, min_samples={n_min_samples:.0f}" if has_result else "no result"
            ari_str = f"ARI={ari:.3f}, F1={f1:.3f} [P={precision:.2f} R={recall:.2f}]" if has_result else ""
            print(f"domain_scale={domain_scale:6.0f} nm | noise_fraction={info['f_noise']:.2f} | "
                  f"n_domains={packing['n_domains']:4d} | overlap={packing['frac_overlapping']:.0%} | "
                  f"{eps_str} -> {ari_str}")

    results_df = pd.DataFrame.from_records(records)

    if show_plot and len(results_df) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

        axes[0].plot(results_df['domain_scale'], results_df['frac_overlapping'], 'o-', color='steelblue')
        axes[0].set_xlabel('domain_scale [nm]')
        axes[0].set_ylabel('Fraction of domains overlapping\ntheir nearest neighbour')
        axes[0].set_title('Geometric separability\n(diagnose_domain_packing)')
        axes[0].set_ylim(-0.05, 1.05)
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(results_df['domain_scale'], results_df['f1'], 'o-', color='darkorange', label='F1')
        axes[1].plot(results_df['domain_scale'], results_df['ari'], 's--', color='steelblue', alpha=0.7, label='ARI')
        axes[1].set_xlabel('domain_scale [nm]')
        axes[1].set_ylabel('Score')
        axes[1].set_title('DBSCAN recoverability\n(single fraction, no subsampling)')
        axes[1].set_ylim(-0.05, 1.05)
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.suptitle('domain_scale sweep: separability vs. recoverability — look for both curves healthy at once')
        plt.tight_layout()
        plt.show()

    return results_df


# --- ARI / F1 calibration (earlier objective, kept for comparison) ---------------

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


def domain_detection_f1(true_labels, pred_labels, iou_threshold=0.5):
    """
    Domain-level detection F1: matches each true domain to its best-overlapping
    predicted DBSCAN cluster and scores whether that overlap is good enough to
    count the domain as "found", rather than scoring every point pair the way
    ARI does.

    Why this complements ARI (see epsilon_cost_cluster): ARI is a pairwise
    agreement score over ALL points, noise included. With the 15-50% background
    noise typical of generate_nucleus, noise-noise pairs can dominate the pair
    count, and a handful of adjacent domains merging can crater ARI even when
    most individual domains are cleanly recovered — ARI has no notion of "this
    specific domain was found." This metric asks that question directly: of the
    N true domains, how many have a predicted cluster covering most of their
    points and little else (Intersection-over-Union >= iou_threshold — the
    standard object-detection matching convention)? That is insensitive to the
    total noise count and to how many OTHER domains exist, so it does not
    inherit ARI's sensitivity to having many small clusters.

    Noise (label == -1) is excluded on both sides: a true domain can never be
    "matched" to predicted noise, and predicted noise is never scored as a
    spurious cluster (both would otherwise trivially inflate/deflate precision).

    Args:
        true_labels:   (N,) ground-truth cluster_label array (-1 = noise)
        pred_labels:   (N,) DBSCAN-predicted label array (-1 = noise)
        iou_threshold: minimum Intersection-over-Union for a match to count as
                       "found" (0.5 is the standard object-detection convention)

    Returns:
        dict with keys:
            'precision', 'recall', 'f1': domain-detection scores
            'tp', 'fp', 'fn':            matched / spurious / missed domain counts
            'n_true_domains', 'n_pred_clusters': counts excluding noise
    """
    true_labels = np.asarray(true_labels)
    pred_labels = np.asarray(pred_labels)

    true_unique = np.unique(true_labels)
    pred_unique = np.unique(pred_labels)

    C = contingency_matrix(true_labels, pred_labels, sparse=True).tocsr()
    true_sizes = np.asarray(C.sum(axis=1)).ravel()
    pred_sizes = np.asarray(C.sum(axis=0)).ravel()

    pred_noise_idx = np.where(pred_unique == -1)[0]
    pred_noise_col = int(pred_noise_idx[0]) if len(pred_noise_idx) else None

    matched_pred_cols = set()
    tp = 0

    for i, t in enumerate(true_unique):
        if t == -1:
            continue
        row = C.getrow(i)
        if row.nnz == 0:
            continue
        cols, overlaps = row.indices, row.data
        if pred_noise_col is not None:
            keep = cols != pred_noise_col
            cols, overlaps = cols[keep], overlaps[keep]
        if len(cols) == 0:
            continue

        best = np.argmax(overlaps)
        j, overlap = cols[best], overlaps[best]
        union = true_sizes[i] + pred_sizes[j] - overlap
        iou = overlap / union if union > 0 else 0.0
        if iou >= iou_threshold:
            tp += 1
            matched_pred_cols.add(j)

    n_true_domains = int(np.sum(true_unique != -1))
    n_pred_clusters = int(np.sum(pred_unique != -1))
    fn = n_true_domains - tp
    fp = n_pred_clusters - len(matched_pred_cols)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        'precision': precision, 'recall': recall, 'f1': f1,
        'tp': tp, 'fp': fp, 'fn': fn,
        'n_true_domains': n_true_domains, 'n_pred_clusters': n_pred_clusters,
    }


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

        if res.success and (1.0 - res.fun) > 0.1:
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


def _plot_ari_vs_f1(ax_ari, ax_f1, densities, aris, f1s, precisions, recalls, iou_threshold):
    """
    Draws the ARI-vs-domain-detection-F1 diagnostic used by optimize_dbscan_params
    onto a pair of existing axes, so it can be shown either standalone (when too
    few fractions pass ari_threshold for a power-law fit) or as a row within the
    full calibration figure.
    """
    order = np.argsort(densities)
    dens_arr = np.asarray(densities)[order]

    ax_ari.plot(dens_arr, np.asarray(aris)[order], 'o-', color='steelblue')
    ax_ari.set_xlabel('Density [locs/nm²]')
    ax_ari.set_ylabel('ARI')
    ax_ari.set_title('ARI at winning (eps, min_samples)')
    ax_ari.set_ylim(-0.05, 1.05)
    ax_ari.grid(True, alpha=0.3)

    ax_f1.plot(dens_arr, np.asarray(f1s)[order], 'o-', color='darkorange', label='F1')
    ax_f1.plot(dens_arr, np.asarray(precisions)[order], 's--', color='seagreen', alpha=0.6, label='Precision')
    ax_f1.plot(dens_arr, np.asarray(recalls)[order], '^--', color='indianred', alpha=0.6, label='Recall')
    ax_f1.set_xlabel('Density [locs/nm²]')
    ax_f1.set_ylabel('Domain detection score')
    ax_f1.set_title(f'Domain-detection F1 (IoU ≥ {iou_threshold}) at same params')
    ax_f1.set_ylim(-0.05, 1.05)
    ax_f1.legend()
    ax_f1.grid(True, alpha=0.3)


def optimize_dbscan_params(fraction_dfs, calibration_area_nm2, min_samples_range=range(3, 21),
                           use_z=False, ari_threshold=0.3, iou_threshold=0.5,
                           show_plot=True, verbose=True):
    """
    Finds the DBSCAN (epsilon, min_samples) pair that maximises ARI against
    ground-truth cluster_label for each density fraction, then fits power-law
    calibration curves for BOTH parameters as functions of localization density.

    Extends optimize_epsilon() by searching min_samples jointly with epsilon
    instead of holding it fixed. min_samples interacts with epsilon — the
    noise/domain separation a given epsilon achieves shifts as min_samples
    changes — so optimizing epsilon alone at a fixed min_samples can miss the
    true joint optimum. For each candidate min_samples in min_samples_range,
    this runs the same two-stage search as optimize_epsilon (coarse grid, then
    bounded refinement around the coarse winner) and keeps whichever
    (eps, min_samples) pair achieves the highest ARI for that fraction.

    The search itself is still driven by ARI (unchanged objective). At each
    fraction's winning (eps, min_samples), DBSCAN is re-run once more and scored
    with domain_detection_f1 as well, purely for comparison — see that
    function's docstring for why ARI and domain-detection F1 can disagree (ARI
    is pairwise and noise-sensitive; F1 asks whether each individual domain was
    recovered). This F1 comparison is computed and reported for EVERY fraction
    that produced a result, independent of ari_threshold — if every fraction's
    ARI happens to fall below threshold (e.g. ARI is uniformly low), you still
    get the F1/precision/recall comparison instead of nothing, which is the
    whole point of checking whether ARI is even the right benchmark here.

    Each DataFrame in fraction_dfs must have columns
    [x [nm], y [nm], z [nm], cluster_label], where cluster_label is the
    ground-truth assigned by spawn_nanodomains / generate_nucleus (use -1 for
    noise). ari_threshold only controls which fractions feed the eps/min_samples
    power-law FIT (same convention as optimize_epsilon) — it does not gate
    whether a fraction gets reported or plotted at all; see 'all_*' below.

    Args:
        fraction_dfs:          list of DataFrames, one per density fraction
        calibration_area_nm2:  nucleus area in nm² used to compute density
        min_samples_range:     iterable of candidate min_samples values to search
        use_z:                 if True, cluster in 3D (XYZ); if False, 2D (XY only)
        ari_threshold:         minimum ARI for a fraction to be kept in the
                                eps/min_samples power-law fit (does not affect
                                whether F1 is computed/reported for that fraction)
        iou_threshold:         IoU threshold passed to domain_detection_f1 for the
                                comparison score (does not affect the ARI-driven search)
        show_plot:             if True, plots the ARI-vs-F1 diagnostic; also plots
                                the eps/min_samples calibration curves if at least
                                2 fractions passed ari_threshold
        verbose:               if True (default), prints per-min_samples progress within
                                each fraction (candidate eps/ARI, running best, elapsed
                                time) in addition to the per-fraction and final summary
                                lines. The full grid can take a while (two-stage epsilon
                                search x every min_samples x every fraction), so this is
                                on by default to make long runs easy to follow. Set to
                                False for a silent run.

    Returns:
        dict with keys:
            'eps_coef', 'eps_exp':                 eps = eps_coef * density^eps_exp
            'min_samples_coef', 'min_samples_exp':  min_samples = coef * density^exp
            'eps_r_squared', 'min_samples_r_squared': R² of each log-log fit
            'densities', 'best_epsilons', 'best_min_samples',
            'best_aris', 'best_f1s', 'best_precisions', 'best_recalls':
                values for fractions that passed ari_threshold (used in the fits)
            'all_densities', 'all_epsilons', 'all_min_samples',
            'all_aris', 'all_f1s', 'all_precisions', 'all_recalls':
                the same, for EVERY fraction that produced a result, regardless
                of ari_threshold — use these to inspect the ARI-vs-F1 comparison
                when few or no fractions passed the threshold
    """
    min_samples_list = list(min_samples_range)
    densities = []
    best_epsilons = []
    best_min_samples_list = []
    best_aris = []
    best_f1s = []
    best_precisions = []
    best_recalls = []

    # Unthresholded: every fraction that produced ANY result, regardless of
    # ari_threshold, so the ARI-vs-F1 comparison survives even when every
    # fraction's ARI falls below threshold.
    all_densities = []
    all_epsilons = []
    all_min_samples = []
    all_aris = []
    all_f1s = []
    all_precisions = []
    all_recalls = []

    dims = "3D (XYZ)" if use_z else "2D (XY)"
    n_combos = len(fraction_dfs) * len(min_samples_list)
    if verbose:
        print(f"Starting joint (eps, min_samples) ARI optimisation over {len(fraction_dfs)} fractions "
              f"[{dims}], min_samples in [{min_samples_list[0]}, {min_samples_list[-1]}] "
              f"({n_combos} (fraction, min_samples) combinations total)...")

    combo_count = 0
    for i, df in enumerate(fraction_dfs):
        coords = df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
        true_labels = df["cluster_label"].to_numpy()

        if len(coords) < min(min_samples_list):
            if verbose:
                print(f"Fraction {i+1}/{len(fraction_dfs)}: Skipped (too few localisations)")
            combo_count += len(min_samples_list)
            continue

        density = len(coords) / calibration_area_nm2
        if verbose:
            print(f"Fraction {i+1}/{len(fraction_dfs)}: {len(coords):,} locs, "
                  f"density {density:.6e} locs/nm²")

        frac_best_eps = None
        frac_best_min_samples = None
        frac_best_ari = -1.0
        frac_start = time.perf_counter()

        for min_samples in min_samples_list:
            combo_count += 1

            if len(coords) < min_samples:
                if verbose:
                    print(f"    [{combo_count}/{n_combos}] min_samples={min_samples:2d}: "
                          f"skipped (fewer locs than min_samples)")
                continue

            # --- Stage 1: coarse grid search over epsilon ---
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

            if not res.success:
                if verbose:
                    print(f"    [{combo_count}/{n_combos}] min_samples={min_samples:2d}: optimisation failed")
                continue

            ari = 1.0 - res.fun
            is_new_best = ari > frac_best_ari
            if verbose:
                marker = "  <- new best for this fraction" if is_new_best else ""
                print(f"    [{combo_count}/{n_combos}] min_samples={min_samples:2d}: "
                      f"eps={res.x:6.2f} nm  ARI={ari:.4f}{marker}")
            if is_new_best:
                frac_best_ari = ari
                frac_best_eps = res.x
                frac_best_min_samples = min_samples

        elapsed = time.perf_counter() - frac_start

        if frac_best_eps is None:
            # every min_samples candidate failed to optimise at all
            if verbose:
                print(f"Fraction {i+1}/{len(fraction_dfs)}: Density {density:.6e} | "
                      f"FAILED (no successful optimisation)  [{elapsed:.1f}s]")
            continue

        # One extra DBSCAN fit at the ARI-winning params, purely to compute the
        # domain-detection F1 comparison score (see domain_detection_f1). Computed
        # unconditionally — not gated by ari_threshold — so a uniformly-low-ARI
        # run still yields a comparison instead of nothing.
        fit_coords = coords[:, :3] if use_z else coords[:, :2]
        pred_labels = DBSCAN(eps=frac_best_eps, min_samples=frac_best_min_samples,
                             n_jobs=-1).fit_predict(fit_coords)
        f1_result = domain_detection_f1(true_labels, pred_labels, iou_threshold=iou_threshold)

        all_densities.append(density)
        all_epsilons.append(frac_best_eps)
        all_min_samples.append(frac_best_min_samples)
        all_aris.append(frac_best_ari)
        all_f1s.append(f1_result['f1'])
        all_precisions.append(f1_result['precision'])
        all_recalls.append(f1_result['recall'])

        passed = frac_best_ari > ari_threshold
        if passed:
            densities.append(density)
            best_epsilons.append(frac_best_eps)
            best_min_samples_list.append(frac_best_min_samples)
            best_aris.append(frac_best_ari)
            best_f1s.append(f1_result['f1'])
            best_precisions.append(f1_result['precision'])
            best_recalls.append(f1_result['recall'])

        if verbose:
            status = "" if passed else f"  | below ari_threshold={ari_threshold}, excluded from power-law fit"
            print(f"Fraction {i+1}/{len(fraction_dfs)}: Density {density:.6e} -> "
                  f"Eps {frac_best_eps:.2f} nm, min_samples {frac_best_min_samples}  "
                  f"(ARI {frac_best_ari:.4f}, domain F1 {f1_result['f1']:.4f} "
                  f"[P={f1_result['precision']:.2f} R={f1_result['recall']:.2f}, "
                  f"{f1_result['tp']} tp / {f1_result['fp']} fp / {f1_result['fn']} fn])"
                  f"{status}  [{elapsed:.1f}s]")

    result = {
        'eps_coef': 0.0, 'eps_exp': 0.0, 'eps_r_squared': 0.0,
        'min_samples_coef': 0.0, 'min_samples_exp': 0.0, 'min_samples_r_squared': 0.0,
        'densities': densities, 'best_epsilons': best_epsilons,
        'best_min_samples': best_min_samples_list, 'best_aris': best_aris,
        'best_f1s': best_f1s, 'best_precisions': best_precisions, 'best_recalls': best_recalls,
        'all_densities': all_densities, 'all_epsilons': all_epsilons, 'all_min_samples': all_min_samples,
        'all_aris': all_aris, 'all_f1s': all_f1s, 'all_precisions': all_precisions, 'all_recalls': all_recalls,
    }

    if len(densities) < 2:
        if verbose:
            print(f"\nNot enough fractions passed ari_threshold={ari_threshold} for a power-law fit "
                  f"({len(densities)}/{len(all_densities)} passed) — no eps/min_samples calibration curve.")
        if show_plot and len(all_densities) > 0:
            fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
            _plot_ari_vs_f1(axes[0], axes[1], all_densities, all_aris, all_f1s,
                            all_precisions, all_recalls, iou_threshold)
            plt.suptitle('ARI vs domain-detection F1 (all fractions — too few passed ari_threshold to fit eps/min_samples)')
            plt.tight_layout()
            plt.show()
        return result

    # --- Fit power laws in log-log space ---
    log_dens = np.log10(densities)

    log_eps = np.log10(best_epsilons)
    eps_par = np.polyfit(log_dens, log_eps, 1)
    eps_exp, eps_coef = eps_par[0], 10 ** eps_par[1]
    eps_r_squared = np.corrcoef(log_dens, log_eps)[0, 1] ** 2

    log_min_samples = np.log10(best_min_samples_list)
    ms_par = np.polyfit(log_dens, log_min_samples, 1)
    ms_exp, ms_coef = ms_par[0], 10 ** ms_par[1]
    ms_r_squared = np.corrcoef(log_dens, log_min_samples)[0, 1] ** 2

    result.update({
        'eps_coef': eps_coef, 'eps_exp': eps_exp, 'eps_r_squared': eps_r_squared,
        'min_samples_coef': ms_coef, 'min_samples_exp': ms_exp, 'min_samples_r_squared': ms_r_squared,
    })

    if verbose:
        print(f"\nCalibration complete:")
        print(f"  Eps         = {eps_coef:.2f} * Density^{eps_exp:.2f}  (R² = {eps_r_squared:.4f})")
        print(f"  min_samples = {ms_coef:.2f} * Density^{ms_exp:.2f}  (R² = {ms_r_squared:.4f})")

    if show_plot:
        dens_arr = np.array(densities)
        dens_smooth = np.linspace(dens_arr.min(), dens_arr.max(), 200)

        fig, axes = plt.subplots(3, 2, figsize=(12, 12))

        axes[0, 0].scatter(log_dens, log_eps, color='blue', label='Optimised epsilons')
        axes[0, 0].plot(log_dens, eps_par[0] * log_dens + eps_par[1], color='red', linestyle='--',
                        label=f'Fit: eps = {eps_coef:.2f} · dens^{eps_exp:.2f}')
        axes[0, 0].set_xlabel('Log10(Density [locs/nm²])')
        axes[0, 0].set_ylabel('Log10(Epsilon [nm])')
        axes[0, 0].set_title('Epsilon Calibration — Log-Log Scale')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        axes[0, 1].scatter(dens_arr, best_epsilons, color='blue', label='Optimised epsilons')
        axes[0, 1].plot(dens_smooth, eps_coef * dens_smooth ** eps_exp, color='red', linestyle='--',
                        label=f'Fit: eps = {eps_coef:.2f} · dens^{eps_exp:.2f}')
        axes[0, 1].set_xlabel('Density [locs/nm²]')
        axes[0, 1].set_ylabel('Epsilon [nm]')
        axes[0, 1].set_title('Epsilon Calibration — Standard Scale')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)

        axes[1, 0].scatter(log_dens, log_min_samples, color='blue', label='Optimised min_samples')
        axes[1, 0].plot(log_dens, ms_par[0] * log_dens + ms_par[1], color='red', linestyle='--',
                        label=f'Fit: min_samples = {ms_coef:.2f} · dens^{ms_exp:.2f}')
        axes[1, 0].set_xlabel('Log10(Density [locs/nm²])')
        axes[1, 0].set_ylabel('Log10(min_samples)')
        axes[1, 0].set_title('min_samples Calibration — Log-Log Scale')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)

        axes[1, 1].scatter(dens_arr, best_min_samples_list, color='blue', label='Optimised min_samples')
        axes[1, 1].plot(dens_smooth, ms_coef * dens_smooth ** ms_exp, color='red', linestyle='--',
                        label=f'Fit: min_samples = {ms_coef:.2f} · dens^{ms_exp:.2f}')
        axes[1, 1].set_xlabel('Density [locs/nm²]')
        axes[1, 1].set_ylabel('min_samples')
        axes[1, 1].set_title('min_samples Calibration — Standard Scale')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)

        # Diagnostic: does ARI actually track domain-detection F1? Uses ALL
        # fractions that produced a result (not just the ari_threshold-passing
        # subset used for the fits above), so it stays informative even when
        # most fractions failed to pass threshold. If the two diverge (e.g. F1
        # stays high while ARI collapses), that points at ARI's noise-pair-count
        # sensitivity rather than a real drop in cluster quality.
        _plot_ari_vs_f1(axes[2, 0], axes[2, 1], all_densities, all_aris, all_f1s,
                        all_precisions, all_recalls, iou_threshold)

        plt.suptitle(f'DBSCAN Calibration (ARI-optimised, joint eps + min_samples)')
        plt.tight_layout()
        plt.show()

    return result
