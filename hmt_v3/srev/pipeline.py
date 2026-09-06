"""
Top-level SR-EV orchestrator: chains (optional compartment filter) ->
labeling -> photophysics -> rendering -> ground-truth bookkeeping into one
call, from an already-parsed monomers_df through to a written TIFF stack
and a saved ground-truth run directory.

SrEvSimConfig nests every stage's own config dataclass (LabelingConfig,
PhotophysicsConfig, CameraConfig, RenderConfig) so a parameter is declared
once and read from that single instance everywhere it's used -- the
structural fix for the legacy MATLAB pipeline's anti_length bug, where a
parameter was declared in the driver script but the function that should
have read it silently used its own hardcoded value instead.

Scope: one emission channel per run_pipeline call, matching
export_to_thunderstorm's own documented convention of running once per
mark/channel (hmt_v3/simulate.py). A two-channel acquisition (e.g. one dye
on heterochromatin, another on euchromatin, mirroring hmt_v1's dual-color
setup) is two separate run_pipeline calls with different `compartment_filter`
values -- not one call forking internally -- so each channel's labeling
density, dye brightness, and camera noise stay fully independent, exactly
as two real antibody/dye channels would be.

IMPORTANT -- compartment classification vs. spatial cropping order: CN
(coordination number, compartments.compute_coordination_number) counts real
neighbors within a fixed radius, so computing it on an already-cropped patch
undercounts monomers near the crop boundary and biases them toward
euchromatin. SR-EV configs can have 10^5-10^6 monomers, so cropping down to
a small, fast-to-render patch before labeling is typical -- if you also want
compartment_filter to be correct, classify compartments on the FULL,
uncropped monomers_df first (compartments.classify_compartments), THEN crop.
select_labeling_targets/simulate_channel use an existing "compartment"
column as-is rather than recomputing it, so this ordering is respected
automatically once you've done it once.
"""
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from . import compartments as srev_compartments
from . import ground_truth as srev_ground_truth
from . import labeling as srev_labeling
from . import photophysics as srev_photophysics
from . import render as srev_render


@dataclass
class SrEvSimConfig:
    """
    One run's full configuration. n_frames lives ONLY on `photophysics`
    (PhotophysicsConfig.n_frames) and is reused by run_pipeline for both
    the blinking simulation length and the number of frames written to
    the TIFF stack -- a single source of truth, not two independently
    settable frame counts that could silently diverge.

    compartment_filter: None (label every monomer), "heterochromatin", or
    "euchromatin" -- see module docstring for the classify-before-crop
    ordering requirement if you want this to be accurate on a cropped
    monomers_df.
    """
    labeling: srev_labeling.LabelingConfig = field(default_factory=srev_labeling.LabelingConfig)
    photophysics: srev_photophysics.PhotophysicsConfig = field(default_factory=srev_photophysics.PhotophysicsConfig)
    camera: srev_render.CameraConfig = field(default_factory=srev_render.CameraConfig)
    render: srev_render.RenderConfig = field(default_factory=srev_render.RenderConfig)

    compartment_filter: str = None
    coordination_radius_nm: float = 11.5
    heterochromatin_cn_threshold: int = 6
    frame_margin_nm: float = 300.0
    seed: int = None


def select_labeling_targets(monomers_df, config):
    """
    Applies config.compartment_filter, if set. If monomers_df already has
    a "compartment" column (from a prior compartments.classify_compartments
    call -- the recommended path, see module docstring), filters on it
    directly; otherwise computes it on whatever monomers_df was passed in
    (which undercounts CN near any existing crop boundary -- fine if
    monomers_df is still the full, uncropped cloud).

    Args:
        monomers_df: output of io.parse_config_txt, optionally already
                     passed through compartments.classify_compartments
        config:      SrEvSimConfig

    Returns:
        monomers_df, filtered to the requested compartment (or unchanged
        if compartment_filter is None)
    """
    if config.compartment_filter is None:
        return monomers_df

    if config.compartment_filter not in ("heterochromatin", "euchromatin"):
        raise ValueError(
            "compartment_filter must be None, 'heterochromatin', or 'euchromatin', "
            f"got {config.compartment_filter!r}"
        )

    if "compartment" not in monomers_df.columns:
        monomers_df = srev_compartments.classify_compartments(
            monomers_df, config.coordination_radius_nm, config.heterochromatin_cn_threshold
        )

    return monomers_df[monomers_df["compartment"] == config.compartment_filter]


def simulate_channel(monomers_df, config, rng=None):
    """
    Runs one full channel up through photophysics: compartment filter ->
    labeling -> blinking. Does not render or write anything -- see
    run_pipeline for that.

    Args:
        monomers_df: output of io.parse_config_txt (optionally
                     compartment-labeled already, see module docstring)
        config:      SrEvSimConfig
        rng:         np.random.Generator (default_rng(config.seed) if None)

    Returns:
        dict with keys "monomers" (the filtered, re-indexed monomers
        actually used), "primaries", "secondaries", "fluorophores",
        "blink_events"
    """
    rng = rng if rng is not None else np.random.default_rng(config.seed)

    targets = select_labeling_targets(monomers_df, config)
    if len(targets) == 0:
        raise ValueError(
            f"no monomers left after compartment_filter={config.compartment_filter!r} -- "
            "check coordination_radius_nm/heterochromatin_cn_threshold, or that "
            "monomers_df actually spans both compartments"
        )

    # labeling.place_primaries assumes a dense 0..n-1 monomer_id matching row order
    # (true for io.parse_config_txt's own output); a compartment filter breaks that
    # invariant on the original monomers_df's ids, so it's rebuilt here.
    targets = targets.reset_index(drop=True).copy()
    targets["monomer_id"] = np.arange(len(targets), dtype=np.int64)

    primaries_df, secondaries_df, fluorophores_df = srev_labeling.label_monomers(
        targets, config.labeling, rng=rng
    )
    blink_events_df = srev_photophysics.simulate_blinking(
        fluorophores_df, config.photophysics, rng=rng
    )

    return {
        "monomers": targets,
        "primaries": primaries_df,
        "secondaries": secondaries_df,
        "fluorophores": fluorophores_df,
        "blink_events": blink_events_df,
    }


def run_pipeline(monomers_df, config, output_dir, source_path=None, input_meta=None, rng=None):
    """
    Full single-channel pipeline: simulate_channel's tables -> shift into
    frame-relative coordinates -> render + noise -> write TIFF -> verify
    referential integrity -> save ground truth.

    Args:
        monomers_df: output of io.parse_config_txt (optionally
                     compartment-labeled already, see module docstring)
        config:      SrEvSimConfig
        output_dir:  directory for the TIFF stack + ground-truth tables
                     (created if missing)
        source_path: the config-N.txt path this run came from, if any --
                     recorded in run_config.json purely for provenance,
                     not read
        input_meta:  io.parse_config_txt's returned metadata dict, if any
                     -- recorded in run_config.json alongside source_path
        rng:         np.random.Generator (default_rng(config.seed) if None)

    Returns:
        dict: simulate_channel's tables (post coordinate-shift), plus
        "frame_shape" (the auto-sized (height, width) actually used) and
        "written" (save_run's return value: paths to every written file,
        plus "tiff")
    """
    rng = rng if rng is not None else np.random.default_rng(config.seed)
    output_dir = Path(output_dir)

    tables = simulate_channel(monomers_df, config, rng=rng)

    fluor = tables["fluorophores"]
    frame_shape, x_offset_nm, y_offset_nm = srev_render.fit_frame_to_points(
        fluor["x [nm]"].to_numpy(), fluor["y [nm]"].to_numpy(),
        config.camera.pixel_size_nm, config.frame_margin_nm,
    )
    for name in ("monomers", "primaries", "secondaries", "fluorophores"):
        shifted = tables[name].copy()
        shifted["x [nm]"] += x_offset_nm
        shifted["y [nm]"] += y_offset_nm
        tables[name] = shifted

    camera_config = replace(config.camera, frame_shape=frame_shape)

    output_dir.mkdir(parents=True, exist_ok=True)
    tiff_path = output_dir / "stack.tif"
    srev_render.render_stack_to_tiff(
        tables["fluorophores"], tables["blink_events"], config.photophysics.n_frames,
        camera_config, config.render, str(tiff_path), rng=rng,
    )

    srev_ground_truth.check_referential_integrity(
        tables["monomers"], tables["primaries"], tables["secondaries"],
        tables["fluorophores"], tables["blink_events"],
    )
    written = srev_ground_truth.save_run(
        output_dir, tables["monomers"], tables["primaries"], tables["secondaries"],
        tables["fluorophores"], tables["blink_events"],
        config=config,
        extra_metadata={
            "source_config": str(source_path) if source_path else None,
            "input_meta": input_meta,
            "coordinate_offset_nm": [x_offset_nm, y_offset_nm],
            "frame_shape": list(frame_shape),
        },
    )
    written["tiff"] = str(tiff_path)

    return {**tables, "frame_shape": frame_shape, "written": written}
