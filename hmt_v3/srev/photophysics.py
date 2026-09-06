"""
Fluorophore blinking, as a discrete frame-synchronous multi-state Markov
chain run vectorized across every emitter at once.

State per emitter: Off (dark, not yet or no longer emitting) -> On
(fluorescing) -> Off | Bleached (irreversible). This intentionally omits a
separate short-lived Triplet/dark sub-state (a documented extension point,
not implemented here) to keep the parameter surface to exactly the three
rate constants this module actually needs. A persistent per-emitter state
variable evolving across frames is what gives real temporal correlation
within one molecule's blink train "for free" -- the legacy MATLAB pipeline
instead drew blink-count, frame, and intensity from three independent
empirical marginals, which has no such correlation.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

_OFF, _ON, _BLEACHED = 0, 1, 2


@dataclass
class PhotophysicsConfig:
    """
    Rate constants are in ms^-1 (matching the never-wired-up k2/k3/k4
    convention in the legacy MATLAB master2.m), converted to a per-frame
    transition probability via `1 - exp(-k * frame_time_ms)`.

    mean_photons_per_frame / photon_rate_lognormal_sigma set a per-emitter
    brightness drawn ONCE per fluorophore (lognormal, so brightness is
    always positive and right-skewed like real dye populations) -- this is
    a photon RATE, not a noised count. Poisson shot noise is applied
    exactly once downstream, at the camera stage (render), so it is
    never introduced here.
    """
    k_on_per_ms: float = 1e-5
    k_off_per_ms: float = 0.02
    k_bleach_per_ms: float = 5e-4
    frame_time_ms: float = 10.0
    mean_photons_per_frame: float = 300.0
    photon_rate_lognormal_sigma: float = 0.3
    n_frames: int = 10000


def simulate_blinking(fluorophores_df, config, rng=None):
    """
    Steps the Markov chain across config.n_frames for every fluorophore in
    fluorophores_df, vectorized per frame (one boolean comparison per
    transition type over all emitters at once) -- never a per-emitter
    Python loop, and never a materialized (n_frames, n_emitters) array.
    Only emitters actually ending a blink this frame are appended to a
    compact growing event log, so memory scales with total blink count,
    not frames x emitters.

    On -> {Off, Bleached} is modeled as two competing exponential-rate
    processes within one frame: the combined probability of leaving On is
    `1 - exp(-(k_off + k_bleach) * dt)`, and conditional on leaving, the
    probability the exit was a bleach is `k_bleach / (k_off + k_bleach)` --
    the standard way to combine two competing discrete-time transitions
    without treating them as independent (which would double-count exits).

    Args:
        fluorophores_df: output of labeling.place_fluorophores,
                          needs fluorophore_id
        config:           PhotophysicsConfig
        rng:              np.random.Generator (default_rng() if None)

    Returns:
        DataFrame: blink_id, fluorophore_id, frame_start, frame_end,
                   n_frames_on, total_photons (expected, not yet
                   Poisson-noised), end_reason ("off" | "bleach" |
                   "truncated" -- still On when the movie ended)
    """
    rng = rng if rng is not None else np.random.default_rng()
    n = len(fluorophores_df)
    fluorophore_ids = fluorophores_df["fluorophore_id"].to_numpy()

    mean_photons = rng.lognormal(
        mean=np.log(config.mean_photons_per_frame),
        sigma=config.photon_rate_lognormal_sigma,
        size=n,
    )

    p_on = 1.0 - np.exp(-config.k_on_per_ms * config.frame_time_ms)
    k_leave = config.k_off_per_ms + config.k_bleach_per_ms
    p_leave_on = 1.0 - np.exp(-k_leave * config.frame_time_ms) if k_leave > 0 else 0.0
    p_bleach_given_leave = (config.k_bleach_per_ms / k_leave) if k_leave > 0 else 0.0

    state = np.full(n, _OFF, dtype=np.int8)
    blink_start = np.full(n, -1, dtype=np.int64)
    photons_accum = np.zeros(n, dtype=np.float64)

    ev_fluor, ev_start, ev_end, ev_photons, ev_reason = [], [], [], [], []

    def _close(idx_arr, frame, reason):
        if len(idx_arr) == 0:
            return
        ev_fluor.extend(fluorophore_ids[idx_arr].tolist())
        ev_start.extend(blink_start[idx_arr].tolist())
        ev_end.extend([frame] * len(idx_arr))
        ev_photons.extend(photons_accum[idx_arr].tolist())
        ev_reason.extend([reason] * len(idx_arr))

    for frame in range(config.n_frames):
        activate = state == _OFF
        if p_on > 0 and np.any(activate):
            activate &= rng.random(n) < p_on
            if np.any(activate):
                state[activate] = _ON
                blink_start[activate] = frame

        on_mask = state == _ON
        photons_accum[on_mask] += mean_photons[on_mask]

        if p_leave_on > 0 and np.any(on_mask):
            leaving = on_mask & (rng.random(n) < p_leave_on)
            if np.any(leaving):
                is_bleach = leaving & (rng.random(n) < p_bleach_given_leave)
                is_off = leaving & ~is_bleach

                _close(np.nonzero(is_bleach)[0], frame, "bleach")
                _close(np.nonzero(is_off)[0], frame, "off")

                state[is_bleach] = _BLEACHED
                state[is_off] = _OFF
                blink_start[leaving] = -1
                photons_accum[leaving] = 0.0

    _close(np.nonzero(state == _ON)[0], config.n_frames - 1, "truncated")

    ev_start_arr = np.array(ev_start, dtype=np.int64)
    ev_end_arr = np.array(ev_end, dtype=np.int64)

    return pd.DataFrame({
        "blink_id": np.arange(len(ev_fluor), dtype=np.int64),
        "fluorophore_id": np.array(ev_fluor, dtype=fluorophore_ids.dtype),
        "frame_start": ev_start_arr,
        "frame_end": ev_end_arr,
        "n_frames_on": ev_end_arr - ev_start_arr + 1,
        "total_photons": np.array(ev_photons, dtype=np.float64),
        "end_reason": ev_reason,
    })
