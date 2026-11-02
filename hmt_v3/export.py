"""
Export simulated localizations in a ThunderSTORM-importable CSV (for rendering
in ImageJ / Fiji), with per-localization uncertainty bootstrapped from a real
channel so rendered simulations look like the real data.
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


# --- ThunderSTORM export ---------------------------------------------------------

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

    return out
