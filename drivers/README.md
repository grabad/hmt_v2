# sSMLM histone-modification analysis drivers

Scripts that run the manuscript's analysis ("Quantifying spatial interactions of
post-translational histone modifications using sSMLM") end to end on real data,
using the `hmt_v3` library. Set paths and parameters in `config.py`; the step
scripts don't need editing.

## Input
`DATA_DIR` holds RainbowSTORM exports already split by spectral classification:
`<cell_id>_me3.csv` and `<cell_id>_ac.csv` (columns `x [nm]`, `y [nm]`, `z [nm]`).
The condition comes from the cell_id via `CONDITION_PATTERNS`. By default, any
cell_id containing `act` is ActD and everything else is NT.

## Steps
| Step | Script | What it does | Manuscript |
|---|---|---|---|
| 01 | `01_preprocess.py` | axial filter, nucleus mask (auto threshold), 100 radial contours | Methods; Fig. 3b |
| 02 | `02_calibrate.py`  | inverse Monte Carlo: simulate a matching nucleus with known domains, subsample to 10 densities, fit DBSCAN eps(density) and min_samples(density) | Fig. 4 |
| 03 | `03_cluster.py`    | density-corrected DBSCAN per cell and mark, then per-nanodomain size/count | Fig. 2d/e |
| 04 | `04_spatial.py`    | radial density profiles, SNCR curves, SNCR exponential decay rate | Figs. 3c/f |
| 05 | `05_compare.py`    | pool cells, NT vs ActD statistics, condition figures | Fig. 6, Fig. S3 |

```
python drivers/run_all.py --quick     # smoke test (small calibration, separate output dir)
python drivers/run_all.py             # full run
python drivers/03_cluster.py --cells k27_k27_thaw009 --overwrite   # redo one cell
```
Every step skips cells whose outputs already exist unless you pass `--overwrite`.
Step 02 is slow. Its cost is roughly (#density fractions x #min_samples x ~30
DBSCAN runs) per mark per calibration.

## Check before trusting results
- `cells/<id>/qc_mask.png`: the mask should cover the whole nucleus. Fix bad
  masks with `MASK_THRESH_OVERRIDES`.
- `calibration/<key>/<ch>_pair_correlation.png`: the simulated g(r) and L(r)
  should overlay the real ones.
- `calibration/<key>/<ch>_size_vs_density.png`: corrected sizes should stay
  flat across density, while fixed parameters should drift (Fig. 4b vs 4c).
- `cluster.json` `extrapolated: true`: this cell's density was outside the
  calibrated range.
