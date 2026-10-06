"""
Shared plumbing for the driver scripts: CLI flags, cell discovery, per-cell file
layout, and saving the figures that hmt_v3's plotting functions open.

Output layout (under config.OUTPUT_DIR):
    manifest.csv                     one row per cell: cell_id, condition, input paths
    cells/<cell_id>/                 everything about one nucleus
        me3_masked.parquet, ac_masked.parquet     masked localizations       (01)
        mask.npz                                  mask, contour bands, origin (01)
        preprocess.json, qc_*.png                                             (01)
        <ch>_clustered.parquet, <ch>_nanodomains.csv, cluster.json            (03)
        radial_profile.csv, sncr.csv, sncr_fit.csv, spatial.png               (04)
    calibration/<key>/               one calibration (per condition, cell, or pooled) (02)
        <ch>_calibration.json, <ch>_simulated.parquet, <ch>_size_vs_density.csv, *.png
    results/                         pooled tables, statistics, condition figures (05)
"""
import argparse
import contextlib
import json
import re
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: hmt_v3's plt.show() calls become no-ops, figures stay open to save
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C

warnings.filterwarnings("ignore", message="FigureCanvasAgg is non-interactive")
sys.path.insert(0, str(C.REPO_ROOT))
import hmt_v3 as hmt  # noqa: E402  (needs the repo root on sys.path first)


# --- CLI ---------------------------------------------------------------------

def parse_args(description):
    """Standard flags for every step. Applies --quick overrides to config in place."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--quick", action="store_true",
                        help="smoke-test settings (config.QUICK_OVERRIDES), written to a separate output dir")
    parser.add_argument("--cells", nargs="+", default=None, help="only process these cell_ids")
    parser.add_argument("--overwrite", action="store_true", help="redo cells whose outputs already exist")
    args = parser.parse_args()
    if args.quick:
        for key, value in C.QUICK_OVERRIDES.items():
            setattr(C, key, value)
        log(f"--quick: {', '.join(C.QUICK_OVERRIDES)} overridden; output -> {C.OUTPUT_DIR}")
    return args


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --- paths -------------------------------------------------------------------

def cells_root():
    return Path(C.OUTPUT_DIR) / "cells"


def cell_dir(cell_id):
    d = cells_root() / cell_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def calibration_dir(key):
    d = Path(C.OUTPUT_DIR) / "calibration" / key
    d.mkdir(parents=True, exist_ok=True)
    return d


def results_dir():
    d = Path(C.OUTPUT_DIR) / "results"
    d.mkdir(parents=True, exist_ok=True)
    return d


# --- cells -------------------------------------------------------------------

def assign_condition(cell_id):
    for name, pattern in C.CONDITION_PATTERNS:
        if re.search(pattern, cell_id, flags=re.IGNORECASE):
            return name
    raise ValueError(f"no CONDITION_PATTERNS entry matches {cell_id!r}; add a catch-all ('', r'')")


def discover_cells():
    """Pairs <cell_id>_me3.csv with <cell_id>_ac.csv in DATA_DIR."""
    data_dir = Path(C.DATA_DIR)
    if not data_dir.is_dir():
        raise FileNotFoundError(f"DATA_DIR does not exist: {data_dir} (set it in drivers/config.py)")
    rows = []
    for me3_path in sorted(data_dir.glob(f"*{C.ME3_SUFFIX}")):
        cell_id = me3_path.name[: -len(C.ME3_SUFFIX)]
        ac_path = data_dir / f"{cell_id}{C.AC_SUFFIX}"
        if not ac_path.exists():
            log(f"WARNING: {me3_path.name} has no matching {ac_path.name}; skipped")
            continue
        if cell_id in C.EXCLUDE_CELLS:
            continue
        rows.append({"cell_id": cell_id, "condition": assign_condition(cell_id),
                     "me3_path": str(me3_path), "ac_path": str(ac_path)})
    if not rows:
        raise FileNotFoundError(f"no *{C.ME3_SUFFIX} / *{C.AC_SUFFIX} pairs found in {data_dir}")
    return pd.DataFrame(rows)


def load_manifest(args=None):
    path = Path(C.OUTPUT_DIR) / "manifest.csv"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run 01_preprocess.py first")
    manifest = pd.read_csv(path)
    if args is not None and args.cells:
        missing = set(args.cells) - set(manifest["cell_id"])
        if missing:
            raise ValueError(f"--cells not in manifest: {sorted(missing)}")
        manifest = manifest[manifest["cell_id"].isin(args.cells)]
    return manifest.reset_index(drop=True)


def load_cell(cell_id):
    """Masked localizations + mask geometry written by 01_preprocess."""
    d = cells_root() / cell_id
    geom = np.load(d / "mask.npz")
    return {
        "me3": pd.read_parquet(d / "me3_masked.parquet"),
        "ac": pd.read_parquet(d / "ac_masked.parquet"),
        "binary_mask": geom["binary_mask"],
        "contour_bands": geom["contour_bands"],
        "x_min": float(geom["x_min"]),
        "y_min": float(geom["y_min"]),
        "bin_nm": float(geom["bin_nm"]),
        "area_nm2": float(geom["area_nm2"]),
    }


def calibration_key(cell_id, condition):
    if C.CALIBRATION_MODE == "per_condition":
        return condition
    if C.CALIBRATION_MODE == "per_cell":
        return cell_id
    if C.CALIBRATION_MODE == "pooled":
        return "pooled"
    raise ValueError(f"unknown CALIBRATION_MODE {C.CALIBRATION_MODE!r}")


# --- I/O helpers ---------------------------------------------------------------

def _json_default(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (Path, range)):
        return str(obj)
    raise TypeError(f"not JSON serializable: {type(obj)}")


def save_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=_json_default)


def load_json(path):
    with open(path) as f:
        return json.load(f)


@contextlib.contextmanager
def capture_figures(path_prefix):
    """
    Saves every matplotlib figure opened inside the block as <path_prefix>_<n>.png
    (or <path_prefix>.png if only one), then closes them. Lets the driver reuse
    hmt_v3's plotting functions, which build a figure and call plt.show().
    """
    before = set(plt.get_fignums())
    try:
        yield
    finally:
        new = [n for n in plt.get_fignums() if n not in before]
        for i, num in enumerate(new):
            fig = plt.figure(num)
            suffix = "" if len(new) == 1 else f"_{i + 1}"
            fig.savefig(f"{path_prefix}{suffix}.png", dpi=150, bbox_inches="tight")
            plt.close(fig)


def xyz_cols(use_z):
    return ["x [nm]", "y [nm]", "z [nm]"] if use_z else ["x [nm]", "y [nm]"]
