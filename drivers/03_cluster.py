"""
Step 03 -- density-corrected DBSCAN nanodomain clustering (Methods: "Nanodomain
formation and measurement with DBSCAN").

For each cell and mark:
  1. density = localizations / nucleus mask area (same definition as the
     calibration crop's), then eps and min_samples are read off that cell's
     calibration curves (cluster_dbscan_density_corrected)
  2. DBSCAN -> one cluster per nanodomain, noise = -1
  3. calculate_nanodomain_characteristics: per-nanodomain localization count,
     centroid, bounding-box size (mean of width and height -- the manuscript's
     "nanodomain size"), p75 diameter, hull area/volume, density

Writes cells/<cell_id>/<ch>_clustered.parquet, <ch>_nanodomains.csv and
cluster.json (the eps/min_samples used, and whether this cell's density fell
outside the calibrated range, i.e. the curve was extrapolated).
"""
import config as C
from common import (calibration_dir, calibration_key, cell_dir, hmt, load_cell, load_json,
                    load_manifest, log, parse_args, save_json)

KEEP_COLUMNS = ["id", "frame", "x [nm]", "y [nm]", "z [nm]", "uncertainty [nm]", "centroid [nm]", "cluster_label"]


def cluster_cell(row):
    cell_id, condition = row["cell_id"], row["condition"]
    key = calibration_key(cell_id, condition)
    cell = load_cell(cell_id)
    out = cell_dir(cell_id)
    summary = {"cell_id": cell_id, "condition": condition, "calibration": key}

    for channel in C.CHANNELS:
        cal_path = calibration_dir(key) / f"{channel}_calibration.json"
        if not cal_path.exists():
            raise FileNotFoundError(f"{cal_path} missing -- run 02_calibrate.py")
        calibration = load_json(cal_path)
        if calibration.get("use_z", C.USE_Z) != C.USE_Z:
            raise ValueError(f"{cal_path} was fit with use_z={calibration['use_z']} but config USE_Z={C.USE_Z}")

        clustered, params = hmt.cluster.cluster_dbscan_density_corrected(
            cell[channel], calibration, cell["area_nm2"], use_z=C.USE_Z,
            min_eps=C.MIN_EPS_NM, min_min_samples=C.MIN_MIN_SAMPLES)
        nanodomains = hmt.cluster.calculate_nanodomain_characteristics(clustered, percentile=C.SIZE_PERCENTILE)

        cols = [c for c in KEEP_COLUMNS if c in clustered.columns]
        clustered[cols].to_parquet(out / f"{channel}_clustered.parquet", index=False)
        nanodomains.drop(columns=["cluster_color", "centroid", "vertices"]).to_csv(
            out / f"{channel}_nanodomains.csv", index=False)

        lo, hi = min(calibration["densities"]), max(calibration["densities"])
        tol = C.EXTRAPOLATION_TOLERANCE
        extrapolated = not (lo * (1 - tol) <= params["density"] <= hi * (1 + tol))
        n_clustered = int((clustered["cluster_label"] != -1).sum())
        summary[channel] = {
            **params,
            "density_per_um2": params["density"] * 1e6,
            "calibrated_density_range_per_um2": [lo * 1e6, hi * 1e6],
            "extrapolated": extrapolated,
            "n_locs": len(clustered),
            "frac_locs_clustered": n_clustered / max(len(clustered), 1),
            "n_nanodomains": len(nanodomains),
            "mean_bb_size_nm": float(nanodomains["bb_size_2D"].mean()) if len(nanodomains) else None,
            "median_bb_size_nm": float(nanodomains["bb_size_2D"].median()) if len(nanodomains) else None,
            "mean_diameter_nm": float(nanodomains["geometric_diameter_2D"].mean()) if len(nanodomains) else None,
        }
        flag = "  <-- density outside calibrated range (extrapolated)" if extrapolated else ""
        log(f"{cell_id} [{channel}] density {params['density'] * 1e6:.0f}/um^2 -> eps {params['eps']:.1f} nm, "
            f"min_samples {params['min_samples']}: {len(nanodomains)} nanodomains, "
            f"mean size {summary[channel]['mean_bb_size_nm'] or float('nan'):.0f} nm, "
            f"{summary[channel]['frac_locs_clustered']:.0%} of locs clustered{flag}")

    save_json(summary, out / "cluster.json")


def main():
    args = parse_args(__doc__.splitlines()[1])
    manifest = load_manifest(args)
    for _, row in manifest.iterrows():
        if (cell_dir(row["cell_id"]) / "cluster.json").exists() and not args.overwrite:
            log(f"{row['cell_id']}: already clustered, skipping (--overwrite to redo)")
            continue
        cluster_cell(row)


if __name__ == "__main__":
    main()
