"""
Optional, before step 02 -- choose DOMAIN_SCALE_NM, the simulated domains' cutoff radius.

Step 02's calibration needs simulated ground-truth domains that are separate enough to have a
well-defined size, and that DBSCAN can still recover. For each candidate domain scale this rebuilds
the calibration source cell's simulated nucleus (generate_nucleus), crops it to the calibration crop
and reports, via calibrate.sweep_domain_scale:
    frac_overlapping   fraction of ground-truth domains overlapping their nearest neighbour  (want low)
    f1, ari            how well the best single DBSCAN run recovers those domains             (want high)
It suggests, per mark, the scale with the best F1 among those with at most MAX_OVERLAP overlap.
Set those in config.DOMAIN_SCALE_NM, then rerun 02_calibrate.py with --overwrite.

Writes calibration/<key>/<ch>_domain_scale_sweep.csv and <ch>_domain_scale_sweep.png.
Every candidate rebuilds a whole nucleus, so expect a few minutes per scale and mark.

    python drivers/sweep_domain_scale.py                          # reference condition's cell, both marks
    python drivers/sweep_domain_scale.py --channels me3 --scales 80 120 160
"""
import argparse
import importlib

import numpy as np

import config as C
from common import calibration_dir, calibration_key, capture_figures, hmt, load_cell, load_json, load_manifest, log

calibrate_step = importlib.import_module("02_calibrate")  # reuse its crop choice

DEFAULT_SCALES = [60, 100, 150, 200, 250, 300]
MAX_OVERLAP = 0.5                     # most separable acceptable: at most this fraction of domains overlapping
SWEEP_MIN_SAMPLES = range(5, 55, 10)  # coarse on purpose: this narrows down a scale, step 02 calibrates finely


def densest_cell(manifest, condition):
    """The cell step 02 calibrates a condition on (per_condition mode): its densest nucleus."""
    rows = manifest[manifest["condition"] == condition]
    if rows.empty:
        raise ValueError(f"no cells with condition {condition!r} in the manifest")
    density = [sum(load_json(f"{C.OUTPUT_DIR}/cells/{cell_id}/preprocess.json")[f"{ch}_density_per_um2"]
                   for ch in C.CHANNELS) for cell_id in rows["cell_id"]]
    return rows["cell_id"].iloc[int(np.argmax(density))]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--scales", nargs="+", type=float, default=DEFAULT_SCALES,
                        help="candidate domain scales in nm (default: %(default)s)")
    parser.add_argument("--channels", nargs="+", default=list(C.CHANNELS), choices=list(C.CHANNELS))
    parser.add_argument("--condition", default=C.CONDITION_ORDER[0],
                        help="sweep on this condition's calibration cell (default: %(default)s)")
    parser.add_argument("--cell", default=None, help="sweep on this cell_id instead")
    args = parser.parse_args()

    manifest = load_manifest()
    cell_id = args.cell or densest_cell(manifest, args.condition)
    condition = manifest.set_index("cell_id").loc[cell_id, "condition"]
    out = calibration_dir(calibration_key(cell_id, condition))
    cell = load_cell(cell_id)
    log(f"domain-scale sweep on {cell_id} [{condition}]: scales {args.scales} nm")

    suggested = {}
    for j, ch in enumerate(args.channels):
        # Same crop as the existing calibration if there is one, so the sweep checks what step 02 will use
        cal_path = out / f"{ch}_calibration.json"
        if cal_path.exists():
            bounds = tuple(load_json(cal_path)["crop_bounds_nm"])
        else:
            bounds, _ = calibrate_step.choose_crop(cell[ch], cell["binary_mask"], cell["x_min"], cell["y_min"],
                                                   cell["bin_nm"], C.CROP_SIZE_NM, C.CROP_MIN_INSIDE)
        log(f"  [{ch}] crop {bounds[1] - bounds[0]:.0f} nm square")
        with capture_figures(out / f"{ch}_domain_scale_sweep"):
            table = hmt.calibrate.sweep_domain_scale(
                cell[ch], cell["contour_bands"], cell["x_min"], cell["y_min"], args.scales, bounds,
                noise_fraction=C.NOISE_FRACTION[ch], min_samples_range=SWEEP_MIN_SAMPLES, use_z=C.USE_Z,
                generate_kwargs={"px_size": cell["bin_nm"], "seed_placement": C.SEED_PLACEMENT},
                rng_seed=C.RANDOM_SEED + j, show_plot=True, verbose=True)
        table.to_csv(out / f"{ch}_domain_scale_sweep.csv", index=False)

        ok = table[table["frac_overlapping"] <= MAX_OVERLAP].dropna(subset=["f1"])
        if len(ok):
            best = ok.loc[ok["f1"].idxmax()]
            suggested[ch] = float(best["domain_scale"])
            log(f"  [{ch}] suggested domain scale {suggested[ch]:.0f} nm: {best['frac_overlapping']:.0%} of domains "
                f"overlap, F1 {best['f1']:.2f}")
        else:
            log(f"  [{ch}] WARNING: every scale leaves more than {MAX_OVERLAP:.0%} of the simulated domains "
                f"overlapping (least: {table['frac_overlapping'].min():.0%}), so no scale gives domains with a "
                "well-defined size. See the csv before choosing.")

    if suggested:
        current = dict(C.DOMAIN_SCALE_NM)
        current.update(suggested)
        log(f"set in drivers/config.py:  DOMAIN_SCALE_NM = {current}")
    log(f"tables and plots -> {out}")


if __name__ == "__main__":
    main()
