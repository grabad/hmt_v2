"""
Step 01 -- nucleus masking and radial contours.

For each cell (a <cell_id>_me3.csv / <cell_id>_ac.csv pair in DATA_DIR):
  1. filter_axial: drop localizations piled up at the edges of the z range
  2. binarize_nucleus: smoothed 2D density -> filled mask of the single largest
     nucleus; localizations outside it are dropped
  3. mask_origin: the mask grid's (x_min, y_min) -- every later function that
     maps localizations onto the mask/contours must use this same origin
  4. create_radial_contours: 100 bands following the nuclear outline
     (0 = core ... 99 = periphery), for the radial profiles (Fig. 3b/c, 6c/f)

Writes cells/<cell_id>/{me3,ac}_masked.parquet, mask.npz, preprocess.json, QC
figures, and the manifest.csv every later step reads. Check the qc_mask*.png
files: a bad mask (nucleoli punched out, neighbouring cell included, rim cut
off) should be fixed with MASK_THRESH_OVERRIDES before running step 02.
"""
import numpy as np
import pandas as pd

import config as C
from common import (capture_figures, cell_dir, discover_cells, hmt, log, parse_args,
                    save_json)


def preprocess_cell(row):
    cell_id = row["cell_id"]
    out = cell_dir(cell_id)
    me3_raw = pd.read_csv(row["me3_path"])
    ac_raw = pd.read_csv(row["ac_path"])
    me3_filt = hmt.preprocess.filter_axial(me3_raw, padding=C.AXIAL_PADDING_NM)
    ac_filt = hmt.preprocess.filter_axial(ac_raw, padding=C.AXIAL_PADDING_NM)

    thresh = C.MASK_THRESH_OVERRIDES.get(cell_id, C.MASK_THRESH)
    if thresh == "auto":
        thresh = hmt.preprocess.auto_mask_threshold(
            me3_filt, ac_filt, bin_size=C.MASK_BIN_NM, sigma=C.MASK_SIGMA, factor=C.MASK_AUTO_FACTOR)

    with capture_figures(out / "qc_mask"):
        binary_mask, me3_masked, ac_masked = hmt.preprocess.binarize_nucleus(
            me3_filt, ac_filt, thresh=thresh, bin_size=C.MASK_BIN_NM, sigma=C.MASK_SIGMA,
            num_nuclei=1, show_plots=True)
    x_min, y_min = hmt.preprocess.mask_origin(me3_filt, ac_filt)

    with capture_figures(out / "qc_contours"):
        _, contour_bands = hmt.preprocess.create_radial_contours(
            binary_mask, num_bands=C.N_CONTOUR_BANDS, show_plots=True)

    area_nm2 = float(binary_mask.sum()) * C.MASK_BIN_NM ** 2
    area_um2 = area_nm2 / 1e6

    me3_masked.to_parquet(out / "me3_masked.parquet", index=False)
    ac_masked.to_parquet(out / "ac_masked.parquet", index=False)
    np.savez_compressed(out / "mask.npz", binary_mask=binary_mask, contour_bands=contour_bands,
                        x_min=x_min, y_min=y_min, bin_nm=C.MASK_BIN_NM, area_nm2=area_nm2)

    summary = {
        "cell_id": cell_id, "condition": row["condition"], "mask_thresh": thresh,
        "n_me3_raw": len(me3_raw), "n_ac_raw": len(ac_raw),
        "n_me3": len(me3_masked), "n_ac": len(ac_masked),
        "mask_area_um2": area_um2,
        "me3_density_per_um2": len(me3_masked) / area_um2,
        "ac_density_per_um2": len(ac_masked) / area_um2,
        "me3_ac_ratio": len(me3_masked) / max(len(ac_masked), 1),
    }
    save_json(summary, out / "preprocess.json")

    lo, hi = C.MASK_AREA_WARN_UM2
    flag = "" if lo <= area_um2 <= hi else f"  <-- WARNING mask area outside {lo:g}-{hi:g} um^2, check qc_mask"
    log(f"{cell_id} [{row['condition']}]: mask {area_um2:.0f} um^2 | "
        f"me3 {len(me3_masked):,}/{len(me3_raw):,} | ac {len(ac_masked):,}/{len(ac_raw):,} kept{flag}")
    return summary


def main():
    args = parse_args(__doc__.splitlines()[1])
    manifest = discover_cells()
    if args.cells:
        manifest = manifest[manifest["cell_id"].isin(args.cells)]
    log(f"{len(manifest)} cells in {C.DATA_DIR}: "
        f"{manifest['condition'].value_counts().to_dict()}")

    for _, row in manifest.iterrows():
        done = cell_dir(row["cell_id"]) / "preprocess.json"
        if done.exists() and not args.overwrite:
            log(f"{row['cell_id']}: already preprocessed, skipping (--overwrite to redo)")
            continue
        preprocess_cell(row)

    # The manifest always lists every discovered cell, even when --cells ran a subset.
    full = discover_cells()
    full.to_csv(f"{C.OUTPUT_DIR}/manifest.csv", index=False)
    log(f"manifest -> {C.OUTPUT_DIR}/manifest.csv")


if __name__ == "__main__":
    main()
