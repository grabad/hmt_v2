"""
Ground-truth bookkeeping for one SR-EV simulation run: persisting the
monomer/primary/secondary/fluorophore/blink_event tables as a single run,
and deriving the localization-list projection a real ThunderSTORM run's
output later gets scored against.

Every table already carries its own foreign key -- primaries.monomer_id,
secondaries.primary_id, fluorophores.secondary_id,
blink_events.fluorophore_id -- so this module doesn't invent new keys, it
persists the run as one unit, checks those keys stay consistent, and
builds the one further-denormalized view scoring needs.
"""
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path

import pandas as pd

_TABLE_NAMES = ("monomers", "primaries", "secondaries", "fluorophores", "blink_events")


def resolved_localizations_gt(fluorophores_df, blink_events_df):
    """
    The direct ground-truth comparison target for a ThunderSTORM run on
    this simulation's TIFF stack: one row per blink event (a ThunderSTORM
    detection should correspond to roughly one row here, not one row per
    frame), with the emitter's TRUE static position joined in.

    Args:
        fluorophores_df: output of labeling.place_fluorophores
        blink_events_df: output of photophysics.simulate_blinking

    Returns:
        DataFrame: blink_id, fluorophore_id, frame_start, frame_end,
                   n_frames_on, total_photons, end_reason, x/y/z [nm]
                   (the true position, not what any fitter would recover)
    """
    return blink_events_df.merge(
        fluorophores_df[["fluorophore_id", "x [nm]", "y [nm]", "z [nm]"]],
        on="fluorophore_id", how="left",
    )


def save_run(output_dir, monomers_df, primaries_df, secondaries_df, fluorophores_df,
             blink_events_df, config=None, extra_metadata=None):
    """
    Persists one full SR-EV simulation run to `output_dir`: one Parquet
    file per core table plus resolved_localizations_gt.parquet, and a
    run_config.json capturing the run's configuration for reproducibility.

    Args:
        output_dir:                       directory to write into (created if missing)
        monomers_df ... blink_events_df:  the five core tables
        config:          the run's config object (e.g. a nested
                         SrEvSimConfig dataclass) -- serialized recursively
                         via dataclasses.asdict if it is a dataclass,
                         stored via json's default=str fallback otherwise
        extra_metadata:  additional dict merged into run_config.json
                         (e.g. input config-N.txt path, RNG seed, timestep)

    Returns:
        dict of {table_name: written path}, plus "run_config": its path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = {
        "monomers": monomers_df,
        "primaries": primaries_df,
        "secondaries": secondaries_df,
        "fluorophores": fluorophores_df,
        "blink_events": blink_events_df,
        "resolved_localizations_gt": resolved_localizations_gt(fluorophores_df, blink_events_df),
    }

    written = {}
    for name, df in tables.items():
        path = output_dir / f"{name}.parquet"
        df.to_parquet(path, index=False)
        written[name] = str(path)

    run_config = {"config": asdict(config) if is_dataclass(config) else config}
    if extra_metadata:
        run_config.update(extra_metadata)

    config_path = output_dir / "run_config.json"
    with open(config_path, "w") as f:
        json.dump(run_config, f, indent=2, default=str)
    written["run_config"] = str(config_path)

    return written


def load_run(input_dir):
    """
    Reads back a run saved by save_run.

    Args:
        input_dir: directory previously passed to save_run

    Returns:
        dict of {table_name: DataFrame} for whichever of monomers,
        primaries, secondaries, fluorophores, blink_events,
        resolved_localizations_gt are present, plus "run_config" (the
        parsed run_config.json dict) if present
    """
    input_dir = Path(input_dir)
    tables = {}
    for name in (*_TABLE_NAMES, "resolved_localizations_gt"):
        path = input_dir / f"{name}.parquet"
        if path.exists():
            tables[name] = pd.read_parquet(path)

    config_path = input_dir / "run_config.json"
    if config_path.exists():
        with open(config_path) as f:
            tables["run_config"] = json.load(f)

    return tables


def check_referential_integrity(monomers_df, primaries_df, secondaries_df,
                                 fluorophores_df, blink_events_df):
    """
    Verifies every foreign key in the placement/photophysics chain
    resolves to a real parent row. Raises with the specific orphan count
    on the first violated relationship, rather than passing silently or
    only reporting the first bad row.

    Returns:
        True (never False -- a failed check raises AssertionError)
    """
    checks = [
        ("primaries.monomer_id", primaries_df["monomer_id"], monomers_df["monomer_id"]),
        ("secondaries.primary_id", secondaries_df["primary_id"], primaries_df["primary_id"]),
        ("fluorophores.secondary_id", fluorophores_df["secondary_id"], secondaries_df["secondary_id"]),
        ("blink_events.fluorophore_id", blink_events_df["fluorophore_id"], fluorophores_df["fluorophore_id"]),
    ]
    for label, child_keys, parent_keys in checks:
        orphaned = ~child_keys.isin(parent_keys)
        n_orphaned = int(orphaned.sum())
        assert n_orphaned == 0, (
            f"{label}: {n_orphaned} of {len(child_keys)} rows reference a missing parent"
        )

    return True
