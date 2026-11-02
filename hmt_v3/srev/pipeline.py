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

Two entry points:

run_pipeline -- ONE emission channel per call, matching
export_to_thunderstorm's own documented convention of running once per
mark/channel (hmt_v3/export.py): one compartment_filter, one dye, one
conventional (zeroth-order only) `stack.tif`. Two such calls give two
separate movies that never interact.

run_multichannel_pipeline -- the dual-color sSMLM acquisition hmt_v1's
master2.m actually simulated: every channel (e.g. AF647 on
heterochromatin, CF680 on euchromatin) is labeled and blinked
INDEPENDENTLY, with its own LabelingConfig / PhotophysicsConfig / seed,
exactly as in run_pipeline -- and then all channels' emitters are merged
into ONE emitter set and rendered together into a zeroth-order
`spatial.tif` plus a first-order `spectral.tif` (srev.spectral). Both dyes
blink in the same frames of the same stacks; which is which is only
recoverable from the spectral image, so a downstream RainbowSTORM run is
what separates them again, and the merged ground truth carries each
emitter's `dye` and `emission_centroid_nm` to score that separation
against. Camera/PSF/spectral configs are per-INSTRUMENT (shared across
channels); labeling/photophysics are per-CHANNEL.

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
import json
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np
import pandas as pd

from . import compartments as srev_compartments
from . import ground_truth as srev_ground_truth
from . import labeling as srev_labeling
from . import photophysics as srev_photophysics
from . import render as srev_render
from . import spectral as srev_spectral


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

    dye: the spectral.DyeSpectrum conjugated to this channel's secondary
    antibody (e.g. spectral.AF647 / spectral.CF680). Required by
    run_multichannel_pipeline, where it is what makes the channel
    distinguishable in the spectral stack; optional for run_pipeline,
    which only records it in the fluorophore table if set.

    channel_name: label used for this channel's output subdirectory and
    its "channel" ground-truth column in run_multichannel_pipeline;
    defaults to compartment_filter, then the dye name.
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
    dye: srev_spectral.DyeSpectrum = None
    channel_name: str = None

    def resolved_channel_name(self):
        """channel_name, else compartment_filter, else the dye's name, else "channel"."""
        if self.channel_name:
            return self.channel_name
        if self.compartment_filter:
            return self.compartment_filter
        if self.dye is not None:
            return self.dye.name
        return "channel"


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

    If config.dye is set, fluorophores gains two columns: "dye" (its name)
    and "spectral_shift_nm" (this molecule's spectral-heterogeneity shift,
    Normal(0, dye.heterogeneity_nm), drawn here so it is part of the
    channel's own seeded randomness).

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

    if config.dye is not None:
        fluorophores_df["dye"] = config.dye.name
        fluorophores_df["spectral_shift_nm"] = (
            rng.normal(0.0, config.dye.heterogeneity_nm, size=len(fluorophores_df))
            if config.dye.heterogeneity_nm > 0 else np.zeros(len(fluorophores_df))
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


_SHIFTED_TABLES = ("monomers", "primaries", "secondaries", "fluorophores")


def _shift_tables(tables, x_offset_nm, y_offset_nm):
    """Copies of the coordinate-bearing tables with (x_offset_nm, y_offset_nm) added -- see render.fit_frame_to_points."""
    out = dict(tables)
    for name in _SHIFTED_TABLES:
        if name in out:
            shifted = out[name].copy()
            shifted["x [nm]"] += x_offset_nm
            shifted["y [nm]"] += y_offset_nm
            out[name] = shifted
    return out


def merge_channels(channel_results, channel_configs, spectral_config):
    """
    Combines the per-channel simulate_channel outputs into ONE emitter set
    for joint rendering: a merged fluorophores table and a merged
    blink_events table whose fluorophore_id / blink_id are globally unique
    (each channel's ids are offset by the running count), with the
    channel-level identity every row would otherwise lose:

        channel                 the channel's resolved_channel_name()
        dye                     the dye's name (already set per channel)
        channel_fluorophore_id  the id in that channel's own saved tables
        spectral_shift_nm       the molecule's heterogeneity shift
        emission_centroid_nm    spectral.emission_centroid_nm(dye) +
                                spectral_shift_nm, over
                                spectral_config.wavelength_window_nm --
                                the per-molecule ground truth a
                                RainbowSTORM centroid is scored against

    The per-channel tables themselves are left untouched (their own ids
    stay dense 0..n-1 and self-consistent for ground_truth.save_run).

    Args:
        channel_results: list of simulate_channel outputs, one per channel
        channel_configs: the matching SrEvSimConfig list (dye per channel)
        spectral_config: spectral.SpectralConfig (for the centroid window)

    Returns:
        (fluorophores_df, blink_events_df): the merged tables, both with a
        "channel" column
    """
    fluor_parts, blink_parts = [], []
    fluor_offset, blink_offset = 0, 0
    for result, config in zip(channel_results, channel_configs):
        name = config.resolved_channel_name()
        fl = result["fluorophores"].copy()
        bl = result["blink_events"].copy()
        centroid = srev_spectral.emission_centroid_nm(config.dye, spectral_config)

        fl["channel"] = name
        fl["channel_fluorophore_id"] = fl["fluorophore_id"]
        fl["fluorophore_id"] = fl["fluorophore_id"] + fluor_offset
        fl["emission_centroid_nm"] = centroid + fl["spectral_shift_nm"]

        bl["channel"] = name
        bl["fluorophore_id"] = bl["fluorophore_id"] + fluor_offset
        bl["blink_id"] = bl["blink_id"] + blink_offset

        fluor_parts.append(fl)
        blink_parts.append(bl)
        fluor_offset += len(fl)
        blink_offset += len(bl)

    fluorophores_df = pd.concat(fluor_parts, ignore_index=True)
    blink_events_df = pd.concat(blink_parts, ignore_index=True)
    return fluorophores_df, blink_events_df


def _check_shared_instrument(channel_configs):
    """
    The camera, spatial PSF and frame count are hardware / acquisition
    properties of ONE recording, so every channel rendered into the same
    stacks must agree on them; raises naming the first field that differs.
    """
    ref = channel_configs[0]
    for i, cfg in enumerate(channel_configs[1:], start=1):
        for label, a, b in (
            ("camera", asdict(ref.camera), asdict(cfg.camera)),
            ("render", asdict(ref.render), asdict(cfg.render)),
            ("photophysics.n_frames", ref.photophysics.n_frames, cfg.photophysics.n_frames),
            ("frame_margin_nm", ref.frame_margin_nm, cfg.frame_margin_nm),
        ):
            if a != b:
                raise ValueError(
                    f"channel {i} ({cfg.resolved_channel_name()!r}) has a different {label} than "
                    f"channel 0 ({ref.resolved_channel_name()!r}) -- all channels of one acquisition "
                    "share the instrument, so these must match"
                )


def run_multichannel_pipeline(monomers_df, channel_configs, spectral_config, output_dir,
                              source_path=None, input_meta=None, seed=None):
    """
    Dual-(or N-)color sSMLM acquisition: simulate every channel
    independently, merge their emitters, and render ONE spatial + ONE
    spectral stack in which all dyes blink together (see module
    docstring). The RainbowSTORM-facing outputs in `output_dir`:

        spatial.tif      zeroth-order movie (both dyes; astigmatic PSF, so
                         ThunderSTORM's 3D fit still applies)
        spectral.tif     first-order movie, frame-aligned with spatial.tif
        run_config.json  every channel's config, the SpectralConfig, and
                         under "spectral_calibration" the wavelength ->
                         x-offset relation (parameters + a lookup table) to
                         enter into RainbowSTORM's calibration
        fluorophores.parquet, blink_events.parquet,
        resolved_localizations_gt.parquet
                         MERGED ground truth (merge_channels), with
                         channel / dye / emission_centroid_nm per row --
                         resolved_localizations_gt is the table to match a
                         RainbowSTORM localization list against (one row
                         per blink, true x/y/z, true spectral centroid)
        channels/<name>/ each channel's full ground_truth.save_run output
                         (monomers -> primaries -> secondaries ->
                         fluorophores -> blink_events, ids as in that
                         channel alone), for structure-level analysis

    Every table is shifted into the shared frame's coordinate system
    (render.fit_frame_to_points convention; the offset is recorded in
    run_config.json) -- the same frame for both stacks, sized by
    spectral.fit_frame_to_spatial_and_spectral so dispersed spectra fit.

    Args:
        monomers_df:     output of io.parse_config_txt (compartment-labeled
                         already if channels use compartment_filter -- see
                         module docstring)
        channel_configs: list of SrEvSimConfig, one per channel, each with
                         a `dye` set and distinct resolved_channel_name();
                         camera / render / n_frames must match across them
        spectral_config: spectral.SpectralConfig -- the shared spectral path
        output_dir:      directory for everything above (created if missing)
        source_path:     provenance only, as in run_pipeline
        input_meta:      provenance only, as in run_pipeline
        seed:            RNG seed for the RENDERING noise (each channel's
                         labeling/photophysics uses its own config.seed)

    Returns:
        dict: "channels" ({name: that channel's shifted simulate_channel
        tables}), "fluorophores" / "blink_events" (merged, shifted),
        "frame_shape", "kernels" (per-dye dispersed kernels used),
        "written" (paths to every file written)
    """
    if len(channel_configs) == 0:
        raise ValueError("run_multichannel_pipeline needs at least one channel config")
    names = [cfg.resolved_channel_name() for cfg in channel_configs]
    if len(set(names)) != len(names):
        raise ValueError(f"channel names must be distinct, got {names} -- set channel_name on each config")
    for cfg, name in zip(channel_configs, names):
        if cfg.dye is None:
            raise ValueError(f"channel {name!r} has no dye -- set SrEvSimConfig.dye (e.g. spectral.AF647)")
    _check_shared_instrument(channel_configs)

    output_dir = Path(output_dir)
    rng = np.random.default_rng(seed)
    instrument = channel_configs[0]
    n_frames = instrument.photophysics.n_frames

    results = [
        simulate_channel(monomers_df, cfg, rng=np.random.default_rng(cfg.seed))
        for cfg in channel_configs
    ]
    fluorophores_df, blink_events_df = merge_channels(results, channel_configs, spectral_config)

    dyes = list({cfg.dye.name: cfg.dye for cfg in channel_configs}.values())
    kernels = srev_spectral.build_kernels(dyes, spectral_config)
    frame_shape, x_offset_nm, y_offset_nm = srev_spectral.fit_frame_to_spatial_and_spectral(
        fluorophores_df["x [nm]"].to_numpy(), fluorophores_df["y [nm]"].to_numpy(),
        kernels, instrument.camera, spectral_config, instrument.render, instrument.frame_margin_nm,
    )
    fluorophores_df = _shift_tables({"fluorophores": fluorophores_df}, x_offset_nm, y_offset_nm)["fluorophores"]
    results = [_shift_tables(r, x_offset_nm, y_offset_nm) for r in results]
    camera_config = replace(instrument.camera, frame_shape=frame_shape)

    output_dir.mkdir(parents=True, exist_ok=True)
    spatial_path = output_dir / "spatial.tif"
    spectral_path = output_dir / "spectral.tif"
    srev_spectral.render_dual_stack_to_tiff(
        fluorophores_df, blink_events_df, n_frames, camera_config, instrument.render,
        spectral_config, dyes, str(spatial_path), str(spectral_path), rng=rng,
    )

    written = {"spatial_tiff": str(spatial_path), "spectral_tiff": str(spectral_path), "channels": {}}
    for result, cfg, name in zip(results, channel_configs, names):
        srev_ground_truth.check_referential_integrity(
            result["monomers"], result["primaries"], result["secondaries"],
            result["fluorophores"], result["blink_events"],
        )
        written["channels"][name] = srev_ground_truth.save_run(
            output_dir / "channels" / name,
            result["monomers"], result["primaries"], result["secondaries"],
            result["fluorophores"], result["blink_events"],
            config=cfg,
            extra_metadata={
                "channel": name,
                "part_of_multichannel_run": str(output_dir),
                "coordinate_offset_nm": [x_offset_nm, y_offset_nm],
                "frame_shape": list(frame_shape),
            },
        )

    assert blink_events_df["fluorophore_id"].isin(fluorophores_df["fluorophore_id"]).all(), (
        "merged blink_events reference a fluorophore_id missing from the merged fluorophores"
    )
    merged_tables = {
        "fluorophores": fluorophores_df,
        "blink_events": blink_events_df,
        "resolved_localizations_gt": srev_ground_truth.resolved_localizations_gt(
            fluorophores_df, blink_events_df
        ),
    }
    for name, df in merged_tables.items():
        path = output_dir / f"{name}.parquet"
        df.to_parquet(path, index=False)
        written[name] = str(path)

    calib_wavelengths = np.arange(spectral_config.wavelength_window_nm[0],
                                  spectral_config.wavelength_window_nm[1] + 1e-9, 10.0)
    run_config = {
        "channels": {name: asdict(cfg) for cfg, name in zip(channel_configs, names)},
        "spectral": asdict(spectral_config),
        "spectral_calibration": {
            "model": "offset_px = reference_offset_px + (wavelength_nm - reference_wavelength_nm) / dispersion_nm_per_px",
            "offset_is_relative_to": "the emitter's zeroth-order x position, in spatial.tif pixels; same y row",
            "dispersion_nm_per_px": spectral_config.dispersion_nm_per_px,
            "reference_wavelength_nm": spectral_config.reference_wavelength_nm,
            "reference_offset_px": spectral_config.reference_offset_px,
            "table": [
                {"wavelength_nm": float(w), "offset_px": float(srev_spectral.wavelength_to_offset_px(w, spectral_config))}
                for w in calib_wavelengths
            ],
        },
        "dyes": {
            dye.name: {**asdict(dye), "emission_centroid_nm": srev_spectral.emission_centroid_nm(dye, spectral_config)}
            for dye in dyes
        },
        "source_config": str(source_path) if source_path else None,
        "input_meta": input_meta,
        "coordinate_offset_nm": [x_offset_nm, y_offset_nm],
        "frame_shape": list(frame_shape),
        "n_frames": n_frames,
        "render_seed": seed,
    }
    config_path = output_dir / "run_config.json"
    with open(config_path, "w") as f:
        json.dump(run_config, f, indent=2, default=str)
    written["run_config"] = str(config_path)

    return {
        "channels": dict(zip(names, results)),
        "fluorophores": fluorophores_df,
        "blink_events": blink_events_df,
        "frame_shape": frame_shape,
        "kernels": kernels,
        "written": written,
    }
