"""
Primary/secondary antibody + fluorophore placement on an SR-EV monomer cloud.

Each level (primary, secondary, fluorophore) is built the same way: an
anchor point (the parent's position) plus a random offset (a uniform
direction on the sphere times a linker length), proposed for a whole wave
of candidates at once and checked for steric clashes in batched
`scipy.spatial.cKDTree` queries rather than a per-pair Python loop. Both
"how many secondaries per primary" and "how many fluorophores per
secondary" are treated as REQUESTED counts — steric rejection is what
determines the REALIZED count, so labeling density and physical crowding
interact the way they would in a real sample instead of a hand-picked
constant always being satisfied.
"""
from dataclasses import dataclass, field
from typing import Union, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

CountSpec = Union[int, Tuple[str, float], Tuple[str, int, float]]


@dataclass
class LabelingConfig:
    """
    All tunable knobs for antibody/fluorophore placement, meant to be
    created once and threaded unchanged through place_primaries /
    place_secondaries / place_fluorophores / label_monomers — see
    pipeline.SrEvSimConfig, which nests this alongside the other
    stage configs so no parameter can be declared in one place and read
    from a different, silently-diverged default somewhere else (the bug
    found in the legacy MATLAB pipeline's `anti_length`).

    n_secondary_per_primary / n_dye_per_secondary accept:
        int                    -> a fixed requested count
        ("poisson", lam)       -> Poisson(lam), clipped to max_count_per_parent
        ("binomial", n, p)     -> Binomial(n, p), clipped to max_count_per_parent
    """
    labeling_efficiency: float = 0.25
    primary_linker_length_nm: float = 10.0
    primary_linker_jitter_nm: float = 2.0
    n_secondary_per_primary: CountSpec = 1
    secondary_linker_length_nm: float = 10.0
    secondary_linker_jitter_nm: float = 2.0
    n_dye_per_secondary: CountSpec = field(default_factory=lambda: ("poisson", 3.0))
    dye_offset_length_nm: float = 3.0
    dye_offset_jitter_nm: float = 1.0
    exclusion_radius_nm: float = 5.0
    max_resample_rounds: int = 6
    max_count_per_parent: int = 8


def _sample_uniform_sphere(n, rng):
    """
    n uniform-random unit vectors on the sphere, via normalize(N(0, I_3))
    (Marsaglia's method) — NOT independent theta/phi draws. Independent
    draws over theta in [0, 2pi) and phi in some [a, b] (as the legacy
    MATLAB sim_anti.m did, with mismatched, asymmetric ranges for the two
    antibody segments) oversample the poles; this method is exactly
    uniform on the sphere regardless of range choices.
    """
    v = rng.normal(size=(n, 3))
    norms = np.linalg.norm(v, axis=1, keepdims=True)
    norms[norms == 0] = 1.0  # measure-zero exact-origin draw; avoid a NaN direction
    return v / norms


def _sample_linker_length(n, mean_nm, jitter_nm, rng):
    """Truncated-normal linker length, floored well above zero so a wild low tail can't collapse or invert the offset."""
    lengths = rng.normal(mean_nm, jitter_nm, size=n)
    return np.clip(lengths, 0.1, None)


def _sample_counts(spec, n, rng, max_count):
    """Draws n requested child-counts from a CountSpec (see LabelingConfig), clipped to max_count so a runaway tail can't force unbounded placement attempts."""
    if isinstance(spec, (int, np.integer)):
        return np.full(n, int(spec), dtype=np.int64)
    kind = spec[0]
    if kind == "poisson":
        _, lam = spec
        counts = rng.poisson(lam, size=n)
    elif kind == "binomial":
        _, n_trials, p = spec
        counts = rng.binomial(n_trials, p, size=n)
    else:
        raise ValueError(f"unknown count spec kind: {kind!r}")
    return np.clip(counts, 0, max_count)


def _place_with_steric_rejection(anchor_xyz, parent_id, sample_length_fn, exclusion_radius,
                                  reference_trees, rng, max_rounds=6):
    """
    Proposes one candidate offset point per row of anchor_xyz (the caller
    has already repeated each parent's anchor position once per
    REQUESTED child), and accepts/rejects the whole wave against:
      - each static `reference_trees` entry (already-placed structure
        this wave must not clash with), via `tree.query(candidates, k=1)`
      - every other candidate in the same wave, via
        `cKDTree(candidates).query_pairs(exclusion_radius)` and a greedy
        keep-first thinning pass

    Only the rejected subset is redrawn (fresh direction + length) each
    round, so the pending set shrinks geometrically and total cost stays
    close to O(N log N) rather than the legacy O(N^2) pairwise loop.
    Candidates that are still rejected after `max_rounds` are dropped —
    this is what lets a crowded region realize fewer children than
    requested instead of forcing every request to succeed regardless of
    physical plausibility.

    Args:
        anchor_xyz:      (M, 3) anchor position per requested child
        parent_id:       (M,) id of the requesting parent per row
        sample_length_fn: callable(n, rng) -> (n,) linker lengths, called
                          each round on the current pending-subset size
        exclusion_radius: applied uniformly to both reference-tree and
                          intra-wave checks
        reference_trees: list of prebuilt cKDTree (or None) instances;
                          entries with `tree.n == 0` are skipped
        rng:             np.random.Generator
        max_rounds:      retry budget on the shrinking pending subset

    Returns:
        (accepted_xyz, accepted_parent_id, accepted_anchor_xyz, realized_n):
            accepted_xyz:        (K, 3) accepted candidate positions, K <= M
            accepted_parent_id:  (K,) parent id per accepted row
            accepted_anchor_xyz: (K, 3) the anchor each accepted row came
                                 from (gives linker_length_nm via a norm,
                                 and IS the epitope/parent position callers
                                 need, with no separate lookup required)
            realized_n:          dict {parent_id: accepted-child count},
                                 present for every parent_id in the input
                                 (0 if none of its requests survived)
    """
    n = len(anchor_xyz)
    pending = np.arange(n)
    accepted_xyz, accepted_parent, accepted_anchor = [], [], []
    # Points this same call already accepted in an EARLIER round. Without also
    # checking new candidates against this, two candidates from the same
    # placement level that happen to succeed in different retry rounds are
    # never compared to each other (only against reference_trees and their
    # own round's wave), letting rounds silently clash with one another.
    accepted_so_far = np.empty((0, 3))

    for _round in range(max_rounds):
        if len(pending) == 0:
            break

        cur_anchor = anchor_xyz[pending]
        directions = _sample_uniform_sphere(len(pending), rng)
        lengths = sample_length_fn(len(pending), rng)
        candidates = cur_anchor + directions * lengths[:, None]

        ok = np.ones(len(pending), dtype=bool)
        for tree in reference_trees:
            if tree is not None and tree.n > 0:
                dist, _ = tree.query(candidates, k=1)
                ok &= dist >= exclusion_radius

        if len(accepted_so_far) > 0:
            prev_tree = cKDTree(accepted_so_far)
            dist_prev, _ = prev_tree.query(candidates, k=1)
            ok &= dist_prev >= exclusion_radius

        ok_idx = np.nonzero(ok)[0]
        if len(ok_idx) > 1:
            wave_tree = cKDTree(candidates[ok_idx])
            pairs = wave_tree.query_pairs(exclusion_radius, output_type="ndarray")
            if len(pairs):
                rejected_local = np.zeros(len(ok_idx), dtype=bool)
                for a, b in pairs:
                    if not rejected_local[a]:
                        rejected_local[b] = True
                ok[ok_idx[rejected_local]] = False

        if np.any(ok):
            accepted_xyz.append(candidates[ok])
            accepted_parent.append(parent_id[pending[ok]])
            accepted_anchor.append(cur_anchor[ok])
            accepted_so_far = np.vstack([accepted_so_far, candidates[ok]])

        pending = pending[~ok]

    accepted_xyz = np.concatenate(accepted_xyz, axis=0) if accepted_xyz else np.empty((0, 3))
    accepted_parent = (np.concatenate(accepted_parent, axis=0) if accepted_parent
                        else np.empty((0,), dtype=parent_id.dtype))
    accepted_anchor = np.concatenate(accepted_anchor, axis=0) if accepted_anchor else np.empty((0, 3))

    realized_counts = pd.Series(accepted_parent).value_counts()
    realized_n = {pid: int(realized_counts.get(pid, 0)) for pid in np.unique(parent_id)}

    return accepted_xyz, accepted_parent, accepted_anchor, realized_n


def place_primaries(monomers_df, config, rng=None):
    """
    Labels a `labeling_efficiency` fraction of monomers with a primary
    antibody, offset from its epitope by a random direction times
    `primary_linker_length_nm` (+ jitter), rejected against the full
    monomer backbone and against every other primary in the same wave.

    Assumes monomers_df["monomer_id"] is a dense 0..n-1 range matching
    row order (true for io.parse_config_txt's output), so a
    monomer_id can be used directly as a positional index.

    Args:
        monomers_df: DataFrame from io.parse_config_txt, needs
                     monomer_id, x [nm], y [nm], z [nm]
        config:      LabelingConfig
        rng:         np.random.Generator (default_rng() if None)

    Returns:
        DataFrame: primary_id, monomer_id, x/y/z [nm], epitope_x/y/z [nm]
                   (denormalized copy of the labeled monomer's position),
                   linker_length_nm (realized, i.e. the true distance to
                   the epitope, not the pre-jitter target)
    """
    rng = rng if rng is not None else np.random.default_rng()

    xyz = monomers_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy(dtype=float)
    monomer_ids = monomers_df["monomer_id"].to_numpy()

    is_labeled = rng.random(len(xyz)) < config.labeling_efficiency
    anchor_xyz = xyz[is_labeled]
    parent_id = monomer_ids[is_labeled]

    backbone_tree = cKDTree(xyz)
    length_fn = lambda n, rng: _sample_linker_length(
        n, config.primary_linker_length_nm, config.primary_linker_jitter_nm, rng
    )

    accepted_xyz, accepted_parent, accepted_anchor, _ = _place_with_steric_rejection(
        anchor_xyz, parent_id, length_fn, config.exclusion_radius_nm,
        reference_trees=[backbone_tree], rng=rng, max_rounds=config.max_resample_rounds,
    )

    return pd.DataFrame({
        "primary_id": np.arange(len(accepted_xyz), dtype=np.int64),
        "monomer_id": accepted_parent,
        "x [nm]": accepted_xyz[:, 0],
        "y [nm]": accepted_xyz[:, 1],
        "z [nm]": accepted_xyz[:, 2],
        "epitope_x [nm]": accepted_anchor[:, 0],
        "epitope_y [nm]": accepted_anchor[:, 1],
        "epitope_z [nm]": accepted_anchor[:, 2],
        "linker_length_nm": np.linalg.norm(accepted_xyz - accepted_anchor, axis=1),
    })


def place_secondaries(primaries_df, monomers_df, config, rng=None):
    """
    Requests `n_secondary_per_primary` secondaries per primary, offset by
    `secondary_linker_length_nm` (+ jitter), rejected against the monomer
    backbone, every already-placed primary, and every other secondary in
    the same wave (which is what naturally caps how many secondaries a
    single crowded primary ends up realizing).

    Args:
        primaries_df: output of place_primaries
        monomers_df:  output of io.parse_config_txt (backbone reference)
        config:       LabelingConfig
        rng:          np.random.Generator (default_rng() if None)

    Returns:
        DataFrame: secondary_id, primary_id, x/y/z [nm], linker_length_nm,
                   requested_n (this primary's requested count, repeated
                   per accepted secondary row), realized_n (this primary's
                   final accepted count, repeated per row) — requested_n
                   > realized_n for a row is direct evidence steric
                   rejection actually capped that primary
    """
    rng = rng if rng is not None else np.random.default_rng()

    monomer_xyz = monomers_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy(dtype=float)
    primary_xyz = primaries_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy(dtype=float)
    primary_ids = primaries_df["primary_id"].to_numpy()

    requested = _sample_counts(config.n_secondary_per_primary, len(primary_ids), rng,
                                config.max_count_per_parent)
    anchor_xyz = np.repeat(primary_xyz, requested, axis=0)
    parent_id = np.repeat(primary_ids, requested)

    backbone_tree = cKDTree(monomer_xyz)
    primaries_tree = cKDTree(primary_xyz) if len(primary_xyz) else None
    length_fn = lambda n, rng: _sample_linker_length(
        n, config.secondary_linker_length_nm, config.secondary_linker_jitter_nm, rng
    )

    accepted_xyz, accepted_parent, accepted_anchor, realized_n = _place_with_steric_rejection(
        anchor_xyz, parent_id, length_fn, config.exclusion_radius_nm,
        reference_trees=[backbone_tree, primaries_tree], rng=rng,
        max_rounds=config.max_resample_rounds,
    )

    requested_by_primary = pd.Series(requested, index=primary_ids)
    realized_by_primary = pd.Series(realized_n)

    return pd.DataFrame({
        "secondary_id": np.arange(len(accepted_xyz), dtype=np.int64),
        "primary_id": accepted_parent,
        "x [nm]": accepted_xyz[:, 0],
        "y [nm]": accepted_xyz[:, 1],
        "z [nm]": accepted_xyz[:, 2],
        "linker_length_nm": np.linalg.norm(accepted_xyz - accepted_anchor, axis=1),
        "requested_n": requested_by_primary.reindex(accepted_parent).to_numpy(),
        "realized_n": realized_by_primary.reindex(accepted_parent).to_numpy(),
    })


def place_fluorophores(secondaries_df, primaries_df, monomers_df, config, rng=None):
    """
    Requests `n_dye_per_secondary` fluorophores per secondary (the
    degree-of-labeling knob), offset by `dye_offset_length_nm` (+
    jitter), rejected against the monomer backbone, every primary, every
    secondary, and every other fluorophore in the same wave.

    Args:
        secondaries_df: output of place_secondaries
        primaries_df:   output of place_primaries (backbone reference)
        monomers_df:    output of io.parse_config_txt (backbone reference)
        config:         LabelingConfig
        rng:            np.random.Generator (default_rng() if None)

    Returns:
        DataFrame: fluorophore_id, secondary_id, x/y/z [nm] (the emitter's
                   TRUE static position — the primary scoring target for
                   everything downstream), requested_n/realized_n (same
                   per-parent QC convention as place_secondaries)
    """
    rng = rng if rng is not None else np.random.default_rng()

    monomer_xyz = monomers_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy(dtype=float)
    primary_xyz = primaries_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy(dtype=float)
    secondary_xyz = secondaries_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy(dtype=float)
    secondary_ids = secondaries_df["secondary_id"].to_numpy()

    requested = _sample_counts(config.n_dye_per_secondary, len(secondary_ids), rng,
                                config.max_count_per_parent)
    anchor_xyz = np.repeat(secondary_xyz, requested, axis=0)
    parent_id = np.repeat(secondary_ids, requested)

    backbone_tree = cKDTree(monomer_xyz)
    primaries_tree = cKDTree(primary_xyz) if len(primary_xyz) else None
    secondaries_tree = cKDTree(secondary_xyz) if len(secondary_xyz) else None
    length_fn = lambda n, rng: _sample_linker_length(
        n, config.dye_offset_length_nm, config.dye_offset_jitter_nm, rng
    )

    accepted_xyz, accepted_parent, accepted_anchor, realized_n = _place_with_steric_rejection(
        anchor_xyz, parent_id, length_fn, config.exclusion_radius_nm,
        reference_trees=[backbone_tree, primaries_tree, secondaries_tree], rng=rng,
        max_rounds=config.max_resample_rounds,
    )

    requested_by_secondary = pd.Series(requested, index=secondary_ids)
    realized_by_secondary = pd.Series(realized_n)

    return pd.DataFrame({
        "fluorophore_id": np.arange(len(accepted_xyz), dtype=np.int64),
        "secondary_id": accepted_parent,
        "x [nm]": accepted_xyz[:, 0],
        "y [nm]": accepted_xyz[:, 1],
        "z [nm]": accepted_xyz[:, 2],
        "requested_n": requested_by_secondary.reindex(accepted_parent).to_numpy(),
        "realized_n": realized_by_secondary.reindex(accepted_parent).to_numpy(),
    })


def label_monomers(monomers_df, config, rng=None):
    """
    Runs the full primary -> secondary -> fluorophore placement chain.

    Args:
        monomers_df: output of io.parse_config_txt
        config:      LabelingConfig
        rng:         np.random.Generator (default_rng() if None; pass an
                     explicit seeded Generator for a reproducible run)

    Returns:
        (primaries_df, secondaries_df, fluorophores_df)
    """
    rng = rng if rng is not None else np.random.default_rng()
    primaries_df = place_primaries(monomers_df, config, rng)
    secondaries_df = place_secondaries(primaries_df, monomers_df, config, rng)
    fluorophores_df = place_fluorophores(secondaries_df, primaries_df, monomers_df, config, rng)
    return primaries_df, secondaries_df, fluorophores_df
