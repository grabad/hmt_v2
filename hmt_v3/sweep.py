"""
Toy two-mark centroid model for validating the spatial metrics in
``hmt_v3.postprocess``: synthetic me3 / ac nanodomain centroids with a known
amount of spatial integration, parameter-sweep, density-subsampling and
perturbation experiments, and plots of every metric across an integration sweep.

* generation   -- ``simulate_toy_centroids`` (symmetric integration model) and
                  ``generate_integration_sweep`` (one draw per integration level),
                  wrapped as ``simulate_pair_centroids`` (parametrised by a
                  TOTAL count + an me3/ac balance instead of two raw counts), plus
                  ``simulate_one_sided_centroids`` (asymmetric: one channel
                  is placed independently, the other clusters onto it).
* quantification -- every centroid metric in ``postprocess`` collapsed to a flat
                  dict of scalars by ``summarize_metrics`` so that one (me3, ac)
                  point set -> one DataFrame row.
* plotting     -- ``plot_sncr_sweep``, ``plot_ripley_k_sweep``,
                  ``plot_delaunay_sweep``, ``plot_mst_sweep``,
                  ``plot_merge_curve_sweep``, ``plot_mesh_overlap_sweep``: one
                  metric per figure, one curve / panel per integration level of a
                  ``generate_integration_sweep`` result.

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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection

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

def simulate_toy_centroids(n_domains=150, integration_level=0.5, field_size_nm=10000.0,
                            pair_jitter_nm=60.0, n_domains_ac=None, rng=None):
    """
    Two synthetic sets of nanodomain centroids (me3, ac) inside a square field,
    with a single knob controlling how spatially "integrated" (co-localized)
    the two channels are.

    integration_level in [0, 1]:
        0 -> me3 and ac centroids are two fully independent CSR (uniform
             random) point processes -- no relationship between channels.
        1 -> every centroid in the smaller-count channel is paired to one in
             the other, offset by a small Gaussian jitter (pair_jitter_nm) --
             the channels sit on top of each other as much as their relative
             abundance allows.
        0 < x < 1 -> that fraction of the smaller channel's domains are
             paired; the rest are independent, unpaired CSR centroids in
             each channel.

    Args:
        n_domains:        number of me3 domains
        integration_level: fraction of domains that are spatially paired
        field_size_nm:     side length of the square field (nm)
        pair_jitter_nm:    stdev of the 2D offset between a paired me3/ac centroid
        n_domains_ac:      number of ac domains; None (default) makes it equal
                           to n_domains. Set this different from n_domains to
                           test an abundance imbalance between channels --
                           colocalization_fraction / cross_pair_correlation /
                           self_nonself_contact_ratio's *_norm curves are
                           designed to stay comparable across that imbalance,
                           while the raw self_nonself_contact_ratio curves are
                           NOT (see that function's docstring) -- this is the
                           case that makes the raw-vs-normalized distinction
                           actually visible, unlike the balanced default.
        rng:               np.random.default_rng() instance; None makes a fresh one

    Returns:
        me3_df, ac_df: DataFrames [x [nm], y [nm], pair_id]
            pair_id >= 0 links a paired me3/ac centroid across channels;
            pair_id == -1 marks an independent (unpaired) centroid.
    """
    rng = rng or np.random.default_rng()
    integration_level = float(np.clip(integration_level, 0.0, 1.0))
    n_domains_ac = n_domains if n_domains_ac is None else n_domains_ac

    n_paired = int(round(integration_level * min(n_domains, n_domains_ac)))
    n_me3_indep = n_domains - n_paired
    n_ac_indep = n_domains_ac - n_paired

    parent_xy = rng.uniform(0, field_size_nm, size=(n_paired, 2))
    me3_paired = parent_xy + rng.normal(0, pair_jitter_nm, size=(n_paired, 2))
    ac_paired = parent_xy + rng.normal(0, pair_jitter_nm, size=(n_paired, 2))
    pair_ids = np.arange(n_paired)

    me3_indep = rng.uniform(0, field_size_nm, size=(n_me3_indep, 2))
    ac_indep = rng.uniform(0, field_size_nm, size=(n_ac_indep, 2))

    me3_xy = np.vstack([me3_paired, me3_indep])
    ac_xy = np.vstack([ac_paired, ac_indep])
    me3_pair_id = np.concatenate([pair_ids, np.full(n_me3_indep, -1)])
    ac_pair_id = np.concatenate([pair_ids, np.full(n_ac_indep, -1)])

    me3_df = pd.DataFrame(np.column_stack([me3_xy, me3_pair_id]), columns=["x [nm]", "y [nm]", "pair_id"])
    ac_df = pd.DataFrame(np.column_stack([ac_xy, ac_pair_id]), columns=["x [nm]", "y [nm]", "pair_id"])
    return me3_df, ac_df

def generate_integration_sweep(integration_levels=(0.0, 0.25, 0.5, 0.75, 1.0), n_domains=150,
                               n_domains_ac=None, field_size_nm=10000.0, pair_jitter_nm=60.0, rng=None):
    """
    Generates me3/ac toy centroids at each integration_level in one call, so
    every downstream metric (colocalization_fraction, cross_pair_correlation,
    self_nonself_contact_ratio, or a custom one) is computed on and plotted
    against the SAME generated centroids instead of each caller re-simulating
    its own draw. The balanced case (equal me3/ac counts) and the imbalanced
    case (n_domains_ac != n_domains) both go through this one function --
    just pass n_domains_ac to get the imbalanced distribution.

    Args:
        integration_levels: iterable of integration_level values to sweep
        n_domains, n_domains_ac, field_size_nm, pair_jitter_nm: forwarded to
            simulate_toy_centroids at each level (see its docstring)
        rng: np.random.default_rng() instance, reused (not reset) across
             levels so the whole sweep is reproducible from one seed while
             each level still draws fresh points; None makes a fresh one

    Returns:
        list of dicts, one per level, each with keys:
            'integration_level', 'me3_seeds', 'ac_seeds'
    """
    rng = rng or np.random.default_rng()
    sweep = []
    for level in integration_levels:
        me3_seeds, ac_seeds = simulate_toy_centroids(
            n_domains=n_domains, n_domains_ac=n_domains_ac, integration_level=level,
            field_size_nm=field_size_nm, pair_jitter_nm=pair_jitter_nm, rng=rng)
        sweep.append({"integration_level": level, "me3_seeds": me3_seeds, "ac_seeds": ac_seeds})
    return sweep


def simulate_pair_centroids(n_total=300, me3_fraction=0.5, integration_level=0.5,
                            field_size_nm=10000.0, pair_jitter_nm=60.0, rng=None):
    """
    ``simulate_toy_centroids`` re-parametrised by a TOTAL centroid
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
    return simulate_toy_centroids(
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
    """simulate_pair_centroids called with a BASELINE-style parameter dict."""
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

    Args:
        baseline:        parameter dict with the ``BASELINE`` keys; None -> ``BASELINE``
        grid:            {parameter: values to sweep}; None -> ``DEFAULT_SWEEP_GRID``
        n_replicates:    independent realizations per grid point
        seed:            seed of the single RNG used for the whole sweep
        coloc_radius_nm, reference_radius_nm, r_max, dr, fr_permutations,
        mesh_permutations: forwarded to ``summarize_metrics``
        progress:        print a line as each parameter finishes

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

    Args:
        cluster_fractions: follower fractions that cluster onto the anchor
        followers:         which channel(s) follow: any of "me3", "ac"
        n_replicates:      independent realizations per (follower, fraction)
        seed:              seed of the single RNG used for the whole sweep
        baseline:          parameter dict with the ``BASELINE`` keys; None -> ``BASELINE``
        coloc_radius_nm, reference_radius_nm, r_max, dr, fr_permutations,
        mesh_permutations: forwarded to ``summarize_metrics``

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

    Args:
        me3_df, ac_df:  the dense centroid sets (DataFrames or (N, 2) arrays)
        field_size_nm:  square field side length (nm)
        fractions:      retention fractions to thin both channels to
        n_draws:        random subsets per fraction below 1.0
        seed:           seed of the RNG used for subsampling and permutations
        coloc_radius_nm, reference_radius_nm, r_max, dr, fr_permutations,
        mesh_permutations: forwarded to ``summarize_metrics``

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

    Args:
        me3_df, ac_df:  the dataset to perturb (DataFrames or (N, 2) arrays)
        field_size_nm:  square field side length (nm)
        sigmas_nm:      Gaussian displacement stdevs (nm) to test
        n_repeats:      perturbed copies per sigma
        seed:           seed of the RNG used for jitter and permutations
        clip_to_field:  clip displaced points back into [0, field_size_nm]
        coloc_radius_nm, reference_radius_nm, r_max, dr, fr_permutations,
        mesh_permutations: forwarded to ``summarize_metrics``

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


# --- plotting: metric curves across an integration sweep ------------

def plot_sncr_sweep(sweep, field_size_nm, r_max=800.0, dr=25.0, title=None):
    """
    Computes self_nonself_contact_ratio for every level in `sweep` (as returned
    by generate_integration_sweep) and plots the raw and density-normalized
    me3-centered/ac-centered curves in a 2x2 grid, one color per integration_level.

    Args:
        sweep:      list of dicts from generate_integration_sweep, each with
                    'integration_level', 'me3_seeds', 'ac_seeds'
        field_size_nm, r_max, dr: forwarded to self_nonself_contact_ratio
        title:      figure suptitle; a sensible default is used if None

    Returns:
        (fig, axes): the created matplotlib Figure and 2x2 Axes array
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    cmap = plt.get_cmap("viridis")
    me3_expected = ac_expected = None

    for row in sweep:
        level = row["integration_level"]
        ratios = pp.self_nonself_contact_ratio(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            field_size_nm=field_size_nm, r_max=r_max, dr=dr)
        me3_expected, ac_expected = ratios["me3_expected_ratio"], ratios["ac_expected_ratio"]

        color = cmap(level)
        axes[0, 0].plot(ratios["r"], ratios["me3_ratio_raw"], marker="o", markersize=3, color=color, label=f"level={level}")
        axes[0, 1].plot(ratios["r"], ratios["ac_ratio_raw"], marker="o", markersize=3, color=color, label=f"level={level}")
        axes[1, 0].plot(ratios["r"], ratios["me3_ratio_norm"], marker="o", markersize=3, color=color, label=f"level={level}")
        axes[1, 1].plot(ratios["r"], ratios["ac_ratio_norm"], marker="o", markersize=3, color=color, label=f"level={level}")

    axes[0, 0].axhline(me3_expected, color="gray", linestyle="--", label=f"CSR (chance, raw={me3_expected:.2f})")
    axes[0, 1].axhline(ac_expected, color="gray", linestyle="--", label=f"CSR (chance, raw={ac_expected:.2f})")
    axes[1, 0].axhline(1.0, color="gray", linestyle="--", label="CSR (chance)")
    axes[1, 1].axhline(1.0, color="gray", linestyle="--", label="CSR (chance)")

    axes[0, 0].set_title("me3-centered (raw)")
    axes[0, 1].set_title("ac-centered (raw)")
    axes[1, 0].set_title("me3-centered (normalized)")
    axes[1, 1].set_title("ac-centered (normalized)")

    axes[0, 0].set_ylabel("self/nonself count ratio")
    axes[1, 0].set_ylabel("normalized self/nonself ratio")
    for ax in axes[1, :]:
        ax.set_xlabel("r (nm)")
    for ax in axes.ravel():
        ax.legend(fontsize=7)

    plt.suptitle(title or "Self-nonself contact ratio")
    plt.tight_layout()
    plt.show()
    return fig, axes

def plot_ripley_k_sweep(sweep, field_size_nm, r_max=800.0, dr=25.0, title=None):
    """
    Computes ripleys_k_cross for every level in `sweep` (as returned by
    generate_integration_sweep) and plots the detrended L_ab(r) - r /
    L_ba(r) - r curves, one color per integration_level.

    Args:
        sweep:      list of dicts from generate_integration_sweep, each with
                    'integration_level', 'me3_seeds', 'ac_seeds'
        field_size_nm, r_max, dr: forwarded to ripleys_k_cross
        title:      figure suptitle; a sensible default is used if None

    Returns:
        (fig, axes): the created matplotlib Figure and 1x2 Axes array
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    cmap = plt.get_cmap("viridis")

    for row in sweep:
        level = row["integration_level"]
        k = pp.ripleys_k_cross(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            field_size_nm=field_size_nm, r_max=r_max, dr=dr)

        color = cmap(level)
        axes[0].plot(k["r"], k["L_ab"] - k["r"], marker="o", markersize=3, color=color, label=f"level={level}")
        axes[1].plot(k["r"], k["L_ba"] - k["r"], marker="o", markersize=3, color=color, label=f"level={level}")

    axes[0].axhline(0.0, color="gray", linestyle="--", label="CSR (no association)")
    axes[1].axhline(0.0, color="gray", linestyle="--", label="CSR (no association)")

    axes[0].set_title("me3-centered")
    axes[1].set_title("ac-centered")
    axes[0].set_ylabel("L(r) - r (nm)")
    for ax in axes:
        ax.set_xlabel("r (nm)")
        ax.legend(fontsize=7)

    plt.suptitle(title or "Ripley's K/L")
    plt.tight_layout()
    plt.show()
    return fig, axes

def _draw_graph(ax, points, edges, labels, heterotypic_color="darkorange", homotypic_color="lightgray",
                heterotypic_width=1.2, homotypic_width=0.6):
    """
    Draws a point+edge graph (Delaunay triangulation or MST) on `ax`: nodes
    colored by channel (green=me3, red=ac), edges colored by whether they
    connect the same channel (homotypic, faint gray) or different channels
    (heterotypic, highlighted). Shared by plot_delaunay_sweep and
    plot_mst_sweep since both draw the same kind of point+edge graph.
    """
    heterotypic = labels[edges[:, 0]] != labels[edges[:, 1]]
    segments = points[edges]  # (E, 2, 2): E edges, 2 endpoints, [x, y] each

    ax.add_collection(LineCollection(segments[~heterotypic], colors=homotypic_color,
                                     linewidths=homotypic_width, zorder=1))
    ax.add_collection(LineCollection(segments[heterotypic], colors=heterotypic_color,
                                     linewidths=heterotypic_width, zorder=2))
    ax.scatter(points[labels == 0, 0], points[labels == 0, 1], s=15, color="green", zorder=3, label="me3")
    ax.scatter(points[labels == 1, 0], points[labels == 1, 1], s=15, color="red", zorder=3, label="ac")
    ax.set_aspect("equal")

def plot_delaunay_sweep(sweep, title=None):
    """
    For every level in `sweep` (from generate_integration_sweep), draws the
    Delaunay triangulation over the pooled me3+ac centroids side by side
    across integration_level (heterotypic edges highlighted, homotypic edges
    faint), then plots heterotypic_fraction (raw) and assortativity
    (chance-corrected) vs. integration_level underneath -- see
    delaunay_channel_mixing for what each measures.

    Args:
        sweep: list of dicts from generate_integration_sweep
        title: suptitle for the graph row; a sensible default is used if None

    Returns:
        list of dicts from delaunay_channel_mixing, one per level (sweep order)
    """
    results = [
        pp.delaunay_channel_mixing(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy())
        for row in sweep
    ]

    fig, axes = plt.subplots(1, len(sweep), figsize=(4 * len(sweep), 4), sharey=True)
    axes = np.atleast_1d(axes)
    for row, result, ax in zip(sweep, results, axes):
        _draw_graph(ax, result["points"], result["edges"], result["labels"])
        ax.set_title(f"level={row['integration_level']}\nheterotypic={result['heterotypic_fraction']:.2f}")

    axes[0].legend(loc="upper right", fontsize=9)
    plt.suptitle(title or "Delaunay triangulation, channel mixing")
    plt.tight_layout()
    plt.show()

    levels = [row["integration_level"] for row in sweep]
    fig, axes2 = plt.subplots(1, 2, figsize=(11, 4.5))

    axes2[0].plot(levels, [r["heterotypic_fraction"] for r in results], "o-", color="darkorange")
    axes2[0].set_title("heterotypic edge fraction (raw)")
    axes2[0].set_xlabel("integration_level")
    axes2[0].set_ylim(-0.05, 1.05)

    axes2[1].plot(levels, [r["assortativity"] for r in results], "o-", color="steelblue")
    axes2[1].axhline(0.0, color="gray", linestyle="--", label="no mixing preference")
    axes2[1].set_title("assortativity (chance-corrected)")
    axes2[1].set_xlabel("integration_level")
    axes2[1].set_ylim(-1.05, 1.05)
    axes2[1].legend(fontsize=8)

    plt.tight_layout()
    plt.show()
    return results

def plot_mst_sweep(sweep, n_permutations=500, rng=None, title=None):
    """
    For every level in `sweep` (from generate_integration_sweep), draws the
    minimum spanning tree over the pooled me3+ac centroids side by side
    across integration_level (cross-type edges highlighted), then plots the
    Friedman-Rafsky cross-edge fraction (observed vs. permutation null) and
    z-score vs. integration_level underneath -- see friedman_rafsky_test for
    what each measures.

    Args:
        sweep:          list of dicts from generate_integration_sweep
        n_permutations: forwarded to friedman_rafsky_test
        rng:            np.random.default_rng() instance, reused across
                        levels; None makes a fresh one
        title:          suptitle for the tree row; a sensible default is
                        used if None

    Returns:
        list of dicts from friedman_rafsky_test, one per level (sweep order)
    """
    rng = rng or np.random.default_rng()
    results = [
        pp.friedman_rafsky_test(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            n_permutations=n_permutations, rng=rng)
        for row in sweep
    ]

    fig, axes = plt.subplots(1, len(sweep), figsize=(4 * len(sweep), 4), sharey=True)
    axes = np.atleast_1d(axes)
    for row, result, ax in zip(sweep, results, axes):
        _draw_graph(ax, result["points"], result["edges"], result["labels"],
                    heterotypic_width=1.5, homotypic_width=0.8)
        ax.set_title(f"level={row['integration_level']}\nz={result['z_score']:.1f}")

    axes[0].legend(loc="upper right", fontsize=9)
    plt.suptitle(title or "Minimum spanning tree, channel mixing")
    plt.tight_layout()
    plt.show()

    levels = [row["integration_level"] for row in sweep]
    fig, axes2 = plt.subplots(1, 2, figsize=(11, 4.5))

    axes2[0].plot(levels, [r["observed_fraction"] for r in results], "o-", color="darkorange", label="observed")
    axes2[0].plot(levels, [r["null_mean_fraction"] for r in results], "o--", color="gray", label="permutation null")
    axes2[0].set_title("cross-type edge fraction")
    axes2[0].set_xlabel("integration_level")
    axes2[0].set_ylim(-0.05, 1.05)
    axes2[0].legend(fontsize=8)

    axes2[1].plot(levels, [r["z_score"] for r in results], "o-", color="steelblue")
    axes2[1].axhline(0.0, color="gray", linestyle="--")
    axes2[1].axhline(1.96, color="gray", linestyle=":", linewidth=1, label="|z| = 1.96")
    axes2[1].axhline(-1.96, color="gray", linestyle=":", linewidth=1)
    axes2[1].set_title("Friedman-Rafsky z-score")
    axes2[1].set_xlabel("integration_level")
    axes2[1].legend(fontsize=8)

    plt.tight_layout()
    plt.show()
    return results

def plot_merge_curve_sweep(sweep, r_max=800.0, title=None):
    """
    Computes mst_merge_curve for every level in `sweep` (from
    generate_integration_sweep) and plots the cumulative heterotypic-merge
    fraction vs. merge radius, one color per integration_level -- how EARLY
    (at what length scale) cross-channel merging starts dominating, as a
    complement to plot_mst_sweep's single aggregate number per level.

    Args:
        sweep: list of dicts from generate_integration_sweep
        r_max: forwarded to mst_merge_curve
        title: figure title; a sensible default is used if None

    Returns:
        list of dicts from mst_merge_curve, one per level (sweep order)
    """
    results = []
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("viridis")
    chance = None

    for row in sweep:
        level = row["integration_level"]
        curve = pp.mst_merge_curve(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            r_max=r_max)
        results.append(curve)
        chance = curve["chance_fraction"]

        ax.plot(curve["radius"], curve["cum_heterotypic_frac"], color=cmap(level), label=f"level={level}")

    ax.axhline(chance, color="gray", linestyle="--", label=f"chance ({chance:.2f})")
    ax.set_xlabel("merge radius (nm)")
    ax.set_ylabel("cumulative heterotypic-merge fraction")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=8)
    ax.set_title(title or "MST merge curve (H0 persistence)")
    plt.tight_layout()
    plt.show()
    return results

def plot_mesh_overlap_sweep(sweep, n_permutations=100, rng=None, title=None):
    """
    For every level in `sweep` (from generate_integration_sweep), draws the
    SEPARATE me3 and ac Delaunay meshes overlaid (me3 edges green, ac edges
    red, actual crossing points marked with an orange x) side by side across
    integration_level, then plots the crossing fraction (raw, both
    directions) and -- if n_permutations > 0 -- the permutation z-score vs.
    integration_level underneath. See mesh_overlap for what each measures.

    Args:
        sweep:          list of dicts from generate_integration_sweep
        n_permutations: forwarded to mesh_overlap; each permutation rebuilds
                        BOTH triangulations from scratch, so this is pricier
                        per-call than friedman_rafsky_test's permutations --
                        kept modest by default. 0 skips the null entirely
                        (crossing counts/fractions only, no z-score panel)
        rng:            np.random.default_rng() instance, reused across
                        levels; None makes a fresh one
        title:          suptitle for the mesh row; a sensible default is
                        used if None

    Returns:
        list of dicts from mesh_overlap, one per level (sweep order)
    """
    rng = rng or np.random.default_rng()
    results = [
        pp.mesh_overlap(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            n_permutations=n_permutations, rng=rng)
        for row in sweep
    ]
    has_null = "z_score" in results[0]

    fig, axes = plt.subplots(1, len(sweep), figsize=(4 * len(sweep), 4), sharey=True)
    axes = np.atleast_1d(axes)
    for row, result, ax in zip(sweep, results, axes):
        me3_pts, ac_pts = result["me3_points"], result["ac_points"]
        me3_segs = me3_pts[result["me3_edges"]]
        ac_segs = ac_pts[result["ac_edges"]]

        ax.add_collection(LineCollection(me3_segs, colors="green", linewidths=0.7, alpha=0.6, zorder=1))
        ax.add_collection(LineCollection(ac_segs, colors="red", linewidths=0.7, alpha=0.6, zorder=1))
        ax.scatter(me3_pts[:, 0], me3_pts[:, 1], s=12, color="green", zorder=3, label="me3")
        ax.scatter(ac_pts[:, 0], ac_pts[:, 1], s=12, color="red", zorder=3, label="ac")
        if len(result["crossing_points"]) > 0:
            ax.scatter(result["crossing_points"][:, 0], result["crossing_points"][:, 1],
                      s=4, color="darkorange", marker="x", linewidths=0.7, zorder=4, label="crossing")

        z_label = f"  z={result['z_score']:.1f}" if has_null else ""
        ax.set_title(f"level={row['integration_level']}\ncrossings={result['n_crossings']}{z_label}")
        ax.set_aspect("equal")

    mesh_title = title or "Delaunay mesh overlap, me3 vs. ac"
    axes[0].legend(loc="upper right", fontsize=8)
    plt.suptitle(mesh_title)
    plt.tight_layout()
    plt.show()

    levels = [row["integration_level"] for row in sweep]
    fig, axes2 = plt.subplots(1, 3 if has_null else 1, figsize=(5.5 * (3 if has_null else 1), 4.5))
    axes2 = np.atleast_1d(axes2)

    axes2[0].plot(levels, [r["me3_crossing_fraction"] for r in results], "o-", color="green", label="me3 edges")
    axes2[0].plot(levels, [r["ac_crossing_fraction"] for r in results], "o-", color="red", label="ac edges")
    axes2[0].set_title("fraction of edges crossing the other mesh")
    axes2[0].set_xlabel("integration_level")
    axes2[0].set_ylim(-0.05, 1.05)
    axes2[0].legend(fontsize=8)

    if has_null:
        axes2[1].plot(levels, [r["n_crossings"] for r in results], "o-", color="darkorange", label="observed")
        axes2[1].plot(levels, [r["null_mean_crossings"] for r in results], "o--", color="gray", label="permutation null")
        axes2[1].set_title("total edge crossings")
        axes2[1].set_xlabel("integration_level")
        axes2[1].legend(fontsize=8)

        axes2[2].plot(levels, [r["z_score"] for r in results], "o-", color="steelblue")
        axes2[2].axhline(0.0, color="gray", linestyle="--")
        axes2[2].axhline(1.96, color="gray", linestyle=":", linewidth=1, label="|z| = 1.96")
        axes2[2].axhline(-1.96, color="gray", linestyle=":", linewidth=1)
        axes2[2].set_title("mesh-overlap z-score")
        axes2[2].set_xlabel("integration_level")
        axes2[2].legend(fontsize=8)

    plt.suptitle(f"{mesh_title} -- crossing statistics")
    plt.tight_layout()
    plt.show()
    return results
