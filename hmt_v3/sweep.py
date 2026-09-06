"""
Parameter-sweep, density-subsampling and perturbation experiments for the toy
centroid model in ``hmt_v3.postprocess``.

Everything here is built on the primitives already in ``postprocess``:

* generation   -- ``simulate_toy_centroids`` (symmetric integration model),
                  wrapped here as ``simulate_pair_centroids`` (parametrised by a
                  TOTAL count + an me3/ac balance instead of two raw counts), plus
                  a new ``simulate_one_sided_centroids`` (asymmetric: one channel
                  is placed independently, the other clusters onto it).
* quantification -- every centroid metric in ``postprocess`` collapsed to a flat
                  dict of scalars by ``summarize_metrics`` so that one (me3, ac)
                  point set -> one DataFrame row.

Four experiment drivers, all returning tidy (long-form) DataFrames:

* ``oat_parameter_sweep``          -- one-at-a-time sweep of the five basic
                                     parameters (field size, colocalization
                                     radius = generation pair jitter, integration
                                     level, number of localizations, me3/ac
                                     balance), each varied while all others are
                                     held at a baseline.
* ``one_sided_clustering_sweep``   -- same idea for the asymmetric model: sweep
                                     the fraction of the follower channel that
                                     clusters onto the independently-placed
                                     anchor channel.
* ``density_subsample_experiment`` -- ONE dense realization, randomly subsampled
                                     to a range of retention fractions (many
                                     draws each) so the underlying structure is
                                     fixed and only sampling density changes.
* ``perturbation_experiment``      -- ONE dataset, many copies with every point
                                     independently jittered by a small Gaussian,
                                     to measure metric stability under positional
                                     noise.

``summarize_stability`` reduces any of those to per-group mean / std / CV.
"""

import numpy as np
import pandas as pd

from . import postprocess as pp


# --- shared helpers ---------------------------------------------------------

def _xy(df):
    """Accept a centroid DataFrame (``x [nm]`` / ``y [nm]`` columns) or a raw
    (N, 2) array and always return a float (N, 2) array."""
    if isinstance(df, np.ndarray):
        return np.asarray(df, dtype=float)
    return df[["x [nm]", "y [nm]"]].to_numpy(dtype=float)


#: Every scalar column ``summarize_metrics`` produces, in output order. Used by
#: the drivers to build rows and by ``summarize_stability`` to know which columns
#: are metrics rather than parameters/bookkeeping.
METRIC_COLUMNS = [
    "n_me3", "n_ac",
    "coloc_frac_me3", "coloc_frac_ac",
    "gcross_at_ref", "gcross_peak", "gcross_auc",
    "sncr_me3_norm_at_ref", "sncr_ac_norm_at_ref",
    "ripley_Lab_at_ref", "ripley_Lba_at_ref", "ripley_Lab_max", "ripley_Lba_max",
    "delaunay_heterotypic_fraction", "delaunay_assortativity",
    "fr_observed_fraction", "fr_null_fraction", "fr_z", "fr_p",
    "mst_heterotypic_overall", "mst_heterotypic_minus_chance",
    "mesh_n_crossings", "mesh_me3_crossing_fraction", "mesh_ac_crossing_fraction",
    "mesh_z", "mesh_p",
]


# --- generation -----------------------------------------------------------

def simulate_pair_centroids(n_total=300, me3_fraction=0.5, integration_level=0.5,
                            field_size_nm=10000.0, pair_jitter_nm=60.0, rng=None):
    """
    ``postprocess.simulate_toy_centroids`` re-parametrised by a TOTAL centroid
    count and an me3/ac balance, so "number of localizations" and "balance" are
    two independent knobs rather than two raw per-channel counts.

    Args:
        n_total:           total number of centroids across both channels
        me3_fraction:      fraction of ``n_total`` assigned to me3
                           (0.5 -> balanced, ->1 me3-dominated, ->0 ac-dominated)
        integration_level: passed straight through -- fraction of the smaller
                           channel that is spatially paired to the other
                           (symmetric model: both members of a pair are jittered
                           around a shared parent)
        field_size_nm:     side length of the square field (nm)
        pair_jitter_nm:    stdev of the offset between the two members of a pair
                           -- the generation-time "colocalization radius"
        rng:               ``np.random.default_rng()`` instance; None -> fresh

    Returns:
        me3_df, ac_df -- DataFrames [x [nm], y [nm], pair_id], exactly as
        ``simulate_toy_centroids``.
    """
    rng = rng or np.random.default_rng()
    n_total = int(n_total)
    n_me3 = int(round(float(np.clip(me3_fraction, 0.0, 1.0)) * n_total))
    n_me3 = max(0, min(n_total, n_me3))
    n_ac = n_total - n_me3
    return pp.simulate_toy_centroids(
        n_domains=n_me3, n_domains_ac=n_ac, integration_level=integration_level,
        field_size_nm=field_size_nm, pair_jitter_nm=pair_jitter_nm, rng=rng)


def simulate_one_sided_centroids(n_total=300, me3_fraction=0.5, cluster_fraction=0.5,
                                 field_size_nm=10000.0, cluster_jitter_nm=60.0,
                                 follower="me3", rng=None):
    """
    Asymmetric (one-sided) clustering. The ANCHOR channel is a plain independent
    CSR point process -- it never "knows" about the other channel. The FOLLOWER
    channel places a ``cluster_fraction`` of its points near randomly chosen
    anchor points (isotropic Gaussian offset, ``cluster_jitter_nm``) and the rest
    as independent CSR.

    Contrast with ``simulate_pair_centroids`` / ``simulate_toy_centroids``, where
    integration is symmetric: BOTH members of a pair are jittered around a shared
    hidden parent, so neither channel is independent of the other. Here exactly
    one channel keeps a clean CSR distribution.

    Args:
        n_total, me3_fraction: total count and me3/ac split (as
                               ``simulate_pair_centroids``)
        cluster_fraction:      fraction of the follower channel that clusters onto
                               the anchor (0 -> both channels independent CSR,
                               1 -> every follower point sits on an anchor point)
        field_size_nm:         square field side length (nm)
        cluster_jitter_nm:     stdev of the follower-to-anchor offset (nm)
        follower:              "me3" (ac is the independent anchor, me3 clusters
                               onto it) or "ac" (me3 is the anchor)
        rng:                   ``np.random.default_rng()`` instance; None -> fresh

    Returns:
        me3_df, ac_df -- DataFrames [x [nm], y [nm], pair_id]. On a clustered
        follower point ``pair_id`` is the index of its anchor point (into the
        anchor channel's own DataFrame); it is -1 on independent follower points
        and on every anchor point.
    """
    rng = rng or np.random.default_rng()
    n_total = int(n_total)
    n_me3 = max(0, min(n_total, int(round(float(np.clip(me3_fraction, 0.0, 1.0)) * n_total))))
    n_ac = n_total - n_me3
    cf = float(np.clip(cluster_fraction, 0.0, 1.0))

    if follower == "me3":
        n_follow, n_anchor = n_me3, n_ac
    elif follower == "ac":
        n_follow, n_anchor = n_ac, n_me3
    else:
        raise ValueError("follower must be 'me3' or 'ac'")

    anchor_xy = rng.uniform(0, field_size_nm, size=(n_anchor, 2))

    n_clustered = int(round(cf * n_follow)) if n_anchor > 0 else 0
    n_free = n_follow - n_clustered

    if n_clustered > 0:
        parent_idx = rng.integers(0, n_anchor, size=n_clustered)
        follow_clustered = anchor_xy[parent_idx] + rng.normal(0.0, cluster_jitter_nm, size=(n_clustered, 2))
        clustered_pair_id = parent_idx.astype(float)
    else:
        follow_clustered = np.empty((0, 2))
        clustered_pair_id = np.empty(0)

    follow_free = rng.uniform(0, field_size_nm, size=(n_free, 2))
    follow_xy = np.vstack([follow_clustered, follow_free])
    follow_pair_id = np.concatenate([clustered_pair_id, np.full(n_free, -1.0)])

    anchor_df = pd.DataFrame(
        np.column_stack([anchor_xy, np.full(n_anchor, -1.0)]),
        columns=["x [nm]", "y [nm]", "pair_id"])
    follow_df = pd.DataFrame(
        np.column_stack([follow_xy, follow_pair_id]),
        columns=["x [nm]", "y [nm]", "pair_id"])

    return (follow_df, anchor_df) if follower == "me3" else (anchor_df, follow_df)


# --- unified scalar quantification ---------------------------------------

def summarize_metrics(me3, ac, field_size_nm, *, coloc_radius_nm=150.0,
                      reference_radius_nm=250.0, r_max=800.0, dr=25.0,
                      fr_permutations=200, mesh_permutations=0, rng=None):
    """
    Run every centroid metric in ``hmt_v3.postprocess`` on one (me3, ac) point
    set and collapse each to scalars -- so a whole sweep becomes one DataFrame.

    Any metric that is undefined for the given input (too few points inside the
    border-excluded region, a degenerate triangulation, an empty channel, ...)
    yields ``np.nan`` for its column rather than raising: callers always get a
    full-width row.

    Scalar columns (see ``METRIC_COLUMNS``):
        n_me3, n_ac                          point counts
        coloc_frac_me3 / _ac                 colocalization_fraction, me3->ac / ac->me3
        gcross_at_ref                        cross_pair_correlation g_ab at reference_radius_nm
        gcross_peak                          max of g_ab over r <= r_max
        gcross_auc                           sum((g_ab - 1) * dr) over r <= r_max  (nm; >0 attraction)
        sncr_me3_norm_at_ref / _ac_          density-normalized self/nonself ratio at reference_radius_nm
        ripley_Lab_at_ref / _Lba_            (L(r) - r) at reference_radius_nm, me3-centered / ac-centered
        ripley_Lab_max / _Lba_max            max of (L(r) - r) over r <= r_max
        delaunay_heterotypic_fraction        raw heterotypic edge fraction (pooled Delaunay)
        delaunay_assortativity               Newman's channel assortativity (chance-corrected; <0 mixed)
        fr_observed_fraction / fr_null_fraction  Friedman-Rafsky cross-edge fraction, observed / permutation null
        fr_z / fr_p                          Friedman-Rafsky z-score / two-sided permutation p
        mst_heterotypic_overall              final cumulative heterotypic-merge fraction (MST merge curve)
        mst_heterotypic_minus_chance         same minus the exchangeability chance fraction
        mesh_n_crossings                     separate me3/ac Delaunay meshes: number of edge crossings
        mesh_me3_crossing_fraction / _ac_    fraction of that channel's edges crossing the other mesh
        mesh_z / mesh_p                      NaN unless mesh_permutations > 0

    Args:
        me3, ac:              centroid DataFrames or (N, 2) arrays
        field_size_nm:        square field side length (nm)
        coloc_radius_nm:      radius for colocalization_fraction (analysis side --
                              kept fixed across sweeps by convention)
        reference_radius_nm:  r at which the r-resolved metrics are sampled for
                              their "_at_ref" scalar (snapped to the nearest bin)
        r_max, dr:            radius sweep extent / step for every r-based metric
        fr_permutations:      label permutations for friedman_rafsky_test
                              (< 1 skips it -> fr_* columns NaN)
        mesh_permutations:    permutations for the mesh_overlap null (0 skips it,
                              off by default -- each one rebuilds two
                              triangulations, so it is the pricey knob)
        rng:                  np.random.default_rng() for the permutation nulls

    Returns:
        dict {column: scalar} with exactly the keys in ``METRIC_COLUMNS``.
    """
    rng = rng or np.random.default_rng()
    me3_xy = _xy(me3)
    ac_xy = _xy(ac)
    out = {k: np.nan for k in METRIC_COLUMNS}
    out["n_me3"] = len(me3_xy)
    out["n_ac"] = len(ac_xy)

    def at_ref(r, y):
        r = np.asarray(r, dtype=float)
        y = np.asarray(y, dtype=float)
        if r.size == 0 or not np.any(np.isfinite(y)):
            return np.nan
        return float(y[int(np.argmin(np.abs(r - reference_radius_nm)))])

    try:
        out["coloc_frac_me3"] = pp.colocalization_fraction(me3_xy, ac_xy, coloc_radius_nm)
        out["coloc_frac_ac"] = pp.colocalization_fraction(ac_xy, me3_xy, coloc_radius_nm)
    except Exception:
        pass

    try:
        r, g = pp.cross_pair_correlation(me3_xy, ac_xy, field_size_nm, r_max=r_max, dr=dr)
        out["gcross_at_ref"] = at_ref(r, g)
        if np.any(np.isfinite(g)):
            out["gcross_peak"] = float(np.nanmax(g))
            out["gcross_auc"] = float(np.nansum(np.asarray(g) - 1.0) * dr)
    except Exception:
        pass

    try:
        s = pp.self_nonself_contact_ratio(me3_xy, ac_xy, field_size_nm, r_max=r_max, dr=dr)
        out["sncr_me3_norm_at_ref"] = at_ref(s["r"], s["me3_ratio_norm"])
        out["sncr_ac_norm_at_ref"] = at_ref(s["r"], s["ac_ratio_norm"])
    except Exception:
        pass

    try:
        k = pp.ripleys_k_cross(me3_xy, ac_xy, field_size_nm, r_max=r_max, dr=dr)
        lab = np.asarray(k["L_ab"]) - np.asarray(k["r"])
        lba = np.asarray(k["L_ba"]) - np.asarray(k["r"])
        out["ripley_Lab_at_ref"] = at_ref(k["r"], lab)
        out["ripley_Lba_at_ref"] = at_ref(k["r"], lba)
        if np.any(np.isfinite(lab)):
            out["ripley_Lab_max"] = float(np.nanmax(lab))
        if np.any(np.isfinite(lba)):
            out["ripley_Lba_max"] = float(np.nanmax(lba))
    except Exception:
        pass

    try:
        d = pp.delaunay_channel_mixing(me3_xy, ac_xy)
        out["delaunay_heterotypic_fraction"] = d["heterotypic_fraction"]
        out["delaunay_assortativity"] = d["assortativity"]
    except Exception:
        pass

    if fr_permutations and int(fr_permutations) >= 1:
        try:
            f = pp.friedman_rafsky_test(me3_xy, ac_xy, n_permutations=int(fr_permutations), rng=rng)
            out["fr_observed_fraction"] = f["observed_fraction"]
            out["fr_null_fraction"] = f["null_mean_fraction"]
            out["fr_z"] = f["z_score"]
            out["fr_p"] = f["p_value"]
        except Exception:
            pass

    try:
        m = pp.mst_merge_curve(me3_xy, ac_xy, r_max=r_max)
        chf = np.asarray(m["cum_heterotypic_frac"], dtype=float)
        if chf.size:
            out["mst_heterotypic_overall"] = float(chf[-1])
            out["mst_heterotypic_minus_chance"] = float(chf[-1] - m["chance_fraction"])
    except Exception:
        pass

    try:
        mo = pp.mesh_overlap(me3_xy, ac_xy, n_permutations=int(mesh_permutations), rng=rng)
        out["mesh_n_crossings"] = mo["n_crossings"]
        out["mesh_me3_crossing_fraction"] = mo["me3_crossing_fraction"]
        out["mesh_ac_crossing_fraction"] = mo["ac_crossing_fraction"]
        if "z_score" in mo:
            out["mesh_z"] = mo["z_score"]
            out["mesh_p"] = mo["p_value"]
    except Exception:
        pass

    return out


# --- experiment 1: one-at-a-time parameter sweep -------------------------

#: The point every OAT panel is anchored to. One entry per "basic parameter".
BASELINE = {
    "field_size_nm": 10000.0,   # square field side length (nm)
    "pair_jitter_nm": 60.0,     # generation-time colocalization radius (nm)
    "integration_level": 0.5,   # fraction of the smaller channel spatially paired
    "n_total": 300,             # number of localizations (centroids), both channels
    "me3_fraction": 0.5,        # me3/ac balance (0.5 = balanced)
}

#: Values each parameter is swept through while the other four stay at BASELINE.
DEFAULT_SWEEP_GRID = {
    "field_size_nm":     [5000.0, 7500.0, 10000.0, 12500.0, 15000.0],
    "pair_jitter_nm":    [20.0, 40.0, 60.0, 100.0, 150.0],
    "integration_level": [0.0, 0.25, 0.5, 0.75, 1.0],
    "n_total":           [100, 200, 300, 500, 800],
    "me3_fraction":      [0.3, 0.4, 0.5, 0.6, 0.7],
}


def _params_to_centroids(params, rng):
    return simulate_pair_centroids(
        n_total=params["n_total"], me3_fraction=params["me3_fraction"],
        integration_level=params["integration_level"],
        field_size_nm=params["field_size_nm"], pair_jitter_nm=params["pair_jitter_nm"],
        rng=rng)


def oat_parameter_sweep(baseline=None, grid=None, n_replicates=8, seed=0,
                        coloc_radius_nm=150.0, reference_radius_nm=250.0,
                        r_max=800.0, dr=25.0, fr_permutations=200,
                        mesh_permutations=0, progress=True):
    """
    One-at-a-time sweep of the five basic parameters (see ``BASELINE`` /
    ``DEFAULT_SWEEP_GRID``). For each parameter, every value in its grid is used
    while the other four are held at ``baseline``; ``n_replicates`` independent
    realizations are drawn per grid point from a single seeded RNG (so the whole
    sweep is reproducible and each draw is still fresh).

    The r-based metrics use a fixed absolute ``r_max`` / ``reference_radius_nm``
    across every panel so their scalars stay comparable; in the ``field_size_nm``
    panel the smallest fields therefore keep a smaller fraction of their points
    as border-safe references -- a mild, documented confound of that panel only.

    Returns:
        tidy DataFrame, one row per (swept_param, swept_value, replicate), with
        columns: ``swept_param``, ``swept_value``, ``replicate``, the five
        parameter columns (the swept one overridden), then every column in
        ``METRIC_COLUMNS``.
    """
    baseline = dict(BASELINE if baseline is None else baseline)
    grid = DEFAULT_SWEEP_GRID if grid is None else grid
    rng = np.random.default_rng(seed)
    rows = []
    for pname, values in grid.items():
        for val in values:
            params = dict(baseline)
            params[pname] = val
            for rep in range(n_replicates):
                me3_df, ac_df = _params_to_centroids(params, rng)
                m = summarize_metrics(
                    me3_df, ac_df, field_size_nm=params["field_size_nm"],
                    coloc_radius_nm=coloc_radius_nm, reference_radius_nm=reference_radius_nm,
                    r_max=r_max, dr=dr, fr_permutations=fr_permutations,
                    mesh_permutations=mesh_permutations, rng=rng)
                rows.append({"swept_param": pname, "swept_value": val,
                             "replicate": rep, **params, **m})
        if progress:
            print(f"  oat_parameter_sweep: done {pname}")
    return pd.DataFrame(rows)


# --- experiment 2: one-sided (asymmetric) clustering sweep --------------

def one_sided_clustering_sweep(cluster_fractions=(0.0, 0.25, 0.5, 0.75, 1.0),
                               followers=("me3", "ac"), n_replicates=8, seed=0,
                               baseline=None, coloc_radius_nm=150.0,
                               reference_radius_nm=250.0, r_max=800.0, dr=25.0,
                               fr_permutations=200, mesh_permutations=0):
    """
    Sweep the one-sided model: the follower channel clusters onto an
    independently-placed anchor channel (``simulate_one_sided_centroids``).
    ``cluster_fraction`` is swept for each choice of ``follower``; everything
    else stays at ``baseline`` (``pair_jitter_nm`` is reused as the follower
    cluster jitter). Compare directly against
    ``oat_parameter_sweep`` on ``integration_level`` to see how each metric
    responds to symmetric vs asymmetric integration at the same nominal level.

    Returns:
        tidy DataFrame, one row per (follower, cluster_fraction, replicate),
        with a ``follower`` column, ``cluster_fraction``, ``replicate``, the
        baseline parameter columns, then ``METRIC_COLUMNS``.
    """
    baseline = dict(BASELINE if baseline is None else baseline)
    rng = np.random.default_rng(seed)
    rows = []
    for follower in followers:
        for cf in cluster_fractions:
            for rep in range(n_replicates):
                me3_df, ac_df = simulate_one_sided_centroids(
                    n_total=baseline["n_total"], me3_fraction=baseline["me3_fraction"],
                    cluster_fraction=cf, field_size_nm=baseline["field_size_nm"],
                    cluster_jitter_nm=baseline["pair_jitter_nm"], follower=follower, rng=rng)
                m = summarize_metrics(
                    me3_df, ac_df, field_size_nm=baseline["field_size_nm"],
                    coloc_radius_nm=coloc_radius_nm, reference_radius_nm=reference_radius_nm,
                    r_max=r_max, dr=dr, fr_permutations=fr_permutations,
                    mesh_permutations=mesh_permutations, rng=rng)
                rows.append({"follower": follower, "cluster_fraction": cf,
                             "replicate": rep, **baseline, **m})
    return pd.DataFrame(rows)


# --- experiment 3: density subsampling ---------------------------------

def density_subsample_experiment(me3_df, ac_df, field_size_nm,
                                 fractions=(1.0, 0.8, 0.6, 0.4, 0.25, 0.1, 0.05),
                                 n_draws=12, seed=0, coloc_radius_nm=150.0,
                                 reference_radius_nm=250.0, r_max=800.0, dr=25.0,
                                 fr_permutations=100, mesh_permutations=0):
    """
    Fixed structure, varying sampling density. Takes ONE realized (dense) pair of
    centroid sets and, for each retention ``fraction``, draws ``n_draws``
    independent random subsets without replacement (each channel thinned by the
    same fraction), recomputing ``summarize_metrics`` on each.

    Every draw at a given fraction is a subset of the SAME points, so the spread
    across draws is pure sampling / estimator noise, and any drift of the
    per-fraction mean as the fraction falls is genuine density-dependence of the
    metric. ``fraction == 1.0`` is the full set and is evaluated once.

    Returns:
        tidy DataFrame, one row per (fraction, draw), columns: ``fraction``,
        ``draw``, then ``METRIC_COLUMNS`` (``n_me3`` / ``n_ac`` record the kept
        counts).
    """
    rng = np.random.default_rng(seed)
    me3_xy = _xy(me3_df)
    ac_xy = _xy(ac_df)
    n_me3, n_ac = len(me3_xy), len(ac_xy)
    rows = []
    for frac in fractions:
        frac = float(frac)
        n_draws_eff = 1 if frac >= 1.0 else n_draws
        for draw in range(n_draws_eff):
            if frac >= 1.0:
                sub_me3, sub_ac = me3_xy, ac_xy
            else:
                k_me3 = max(1, int(round(frac * n_me3)))
                k_ac = max(1, int(round(frac * n_ac)))
                sub_me3 = me3_xy[rng.choice(n_me3, size=k_me3, replace=False)]
                sub_ac = ac_xy[rng.choice(n_ac, size=k_ac, replace=False)]
            m = summarize_metrics(
                sub_me3, sub_ac, field_size_nm=field_size_nm,
                coloc_radius_nm=coloc_radius_nm, reference_radius_nm=reference_radius_nm,
                r_max=r_max, dr=dr, fr_permutations=fr_permutations,
                mesh_permutations=mesh_permutations, rng=rng)
            rows.append({"fraction": frac, "draw": draw, **m})
    return pd.DataFrame(rows)


# --- experiment 4: positional perturbation ---------------------------

def perturbation_experiment(me3_df, ac_df, field_size_nm,
                            sigmas_nm=(5.0, 10.0, 20.0, 40.0, 80.0),
                            n_repeats=30, seed=0, clip_to_field=True,
                            coloc_radius_nm=150.0, reference_radius_nm=250.0,
                            r_max=800.0, dr=25.0, fr_permutations=100,
                            mesh_permutations=0):
    """
    Metric stability under positional noise. For each sigma in ``sigmas_nm``,
    makes ``n_repeats`` perturbed copies of the input dataset -- every point
    independently displaced by an isotropic 2D Gaussian of that stdev (optionally
    clipped back into ``[0, field_size_nm]``) -- and recomputes
    ``summarize_metrics`` on each copy. A ``sigma_nm == 0`` group (the untouched
    dataset, one row) is prepended as a noise-free reference.

    Downstream: group by ``sigma_nm`` and read off mean / std / CV per metric
    (``summarize_stability``) to rank the metrics by robustness.

    Returns:
        tidy DataFrame, one row per (sigma_nm, repeat), columns: ``sigma_nm``,
        ``repeat``, then ``METRIC_COLUMNS``.
    """
    rng = np.random.default_rng(seed)
    me3_xy = _xy(me3_df)
    ac_xy = _xy(ac_df)
    sigmas = [0.0] + [float(s) for s in sigmas_nm]
    rows = []
    for sigma in sigmas:
        reps = 1 if sigma == 0.0 else n_repeats
        for rep in range(reps):
            if sigma == 0.0:
                p_me3, p_ac = me3_xy, ac_xy
            else:
                p_me3 = me3_xy + rng.normal(0.0, sigma, size=me3_xy.shape)
                p_ac = ac_xy + rng.normal(0.0, sigma, size=ac_xy.shape)
                if clip_to_field:
                    p_me3 = np.clip(p_me3, 0.0, field_size_nm)
                    p_ac = np.clip(p_ac, 0.0, field_size_nm)
            m = summarize_metrics(
                p_me3, p_ac, field_size_nm=field_size_nm,
                coloc_radius_nm=coloc_radius_nm, reference_radius_nm=reference_radius_nm,
                r_max=r_max, dr=dr, fr_permutations=fr_permutations,
                mesh_permutations=mesh_permutations, rng=rng)
            rows.append({"sigma_nm": sigma, "repeat": rep, **m})
    return pd.DataFrame(rows)


# --- stability reduction ---------------------------------------------

def summarize_stability(df, group_cols, metric_cols=None):
    """
    Per-group mean / std / coefficient-of-variation for every metric column.

    Args:
        df:          any tidy DataFrame from the drivers above
        group_cols:  column(s) defining a group (e.g. ``["sigma_nm"]`` for a
                     perturbation run, ``["fraction"]`` for a density run,
                     ``["swept_param", "swept_value"]`` for an OAT run)
        metric_cols: metric columns to reduce; default is every
                     ``METRIC_COLUMNS`` entry present in ``df``

    Returns:
        wide DataFrame indexed by ``group_cols`` (reset to columns) with
        ``<metric>_mean`` / ``<metric>_std`` / ``<metric>_cv`` for each metric.
        ``cv`` is ``std / |mean|`` (NaN where mean is 0 or the group has one row).
    """
    if metric_cols is None:
        metric_cols = [c for c in METRIC_COLUMNS if c in df.columns]
    group_cols = list(group_cols)
    g = df.groupby(group_cols, dropna=False)[metric_cols]
    agg = g.agg(["mean", "std"])
    agg.columns = [f"{metric}_{stat}" for metric, stat in agg.columns]
    for m in metric_cols:
        mean = agg[f"{m}_mean"]
        std = agg[f"{m}_std"]
        agg[f"{m}_cv"] = std / mean.abs().replace(0.0, np.nan)
    return agg.reset_index()
