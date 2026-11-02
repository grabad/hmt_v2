# legacy

Archived code and notebooks, kept so earlier results can be traced and re-run. Nothing here is
maintained. New work should use `hmt_v3` and the two notebooks in the repo root
(`full_pipeline.ipynb`, `single_channel_nanodomains.ipynb`).

## Code

| File | What it was |
|---|---|
| `hmt_functions.py` | The original single-file toolbox, imported as `hmt` by `notebooks/test_script.ipynb` and `notebooks/cluster_and_view.ipynb` |
| `simulate_OLD.py` | Earlier simulation / calibration functions, previously `hmt_v3.simulate_OLD` |
| `viewer_test.py` | napari viewer test for two raw localization files (run from the repo root: `python legacy/viewer_test.py`) |

Import them as `from legacy import hmt_functions as hmt`, with the repo root on `sys.path`. Both
modules import napari at the top level.

## Notebooks

All of these were in the repo root until October 2026. Each one now starts with a setup cell
that puts the repo root on `sys.path`. Their inputs point back to the repo root (`../../test_data/`,
`../../hmt_v1/`, `../../sr-ev_raw/`), and anything they write lands in `notebooks/` next to
them: `simulated_data/`, the `sr_ev_domain_*` tables and figures, `ovito_export/`,
`ssmlm_export/` and `sr_ev_pipeline_demo_output/`.

| Notebook | Uses | Purpose |
|---|---|---|
| `manuscript_plotting.ipynb` | `hmt_v3` | Manuscript figures: preprocessing, simulated nucleus, DBSCAN calibration, size vs. density |
| `sensitivity analysis.ipynb` | `hmt_v3` | Simulation and clustering sensitivity checks; writes `simulated_data/*_real.csv` |
| `nanodomain analysis wip.ipynb` | `hmt_v3` | Toy-model validation of the two-mark spatial metrics; writes `simulated_data/sweeps/` |
| `hmt_v3 test_script.ipynb` | `hmt_v3` | Scratch tests of the simulation building blocks |
| `sr_ev_pipeline_demo.ipynb` | `hmt_v3.srev` | SR-EV config → labelled, rendered sSMLM movie (single and dual colour) |
| `sr_ev_domain_analysis.ipynb` | `hmt_v3.srev` | Packing-domain counts / radii across SR-EV runs; wrote the `sr_ev_domain_*` files |
| `sr_ev_visualization.ipynb` | `hmt_v3.srev` | OVITO / napari / matplotlib views of SR-EV structures |
| `test_script.ipynb`, `cluster_and_view.ipynb` | `legacy.hmt_functions` | Original exploratory notebooks |

`test_script.ipynb` calls `hmt.spawn_nanodomains` and `hmt.extract_empirical_parameters`, which
`hmt_functions.py` doesn't define (only the `_NEW` / `_OLD` variants exist), so those cells were
already broken before the move.

## Replacements in hmt_v3

| Legacy function | Use instead |
|---|---|
| `filter_quantile`, `filter_axial`, `binarize_nucleus`, `create_radial_contours`, `sample_lower_densities` | `hmt_v3.preprocess.*` (same names) |
| `cluster_dbscan` | `hmt_v3.cluster.cluster_dbscan` (adds `use_z`) |
| `calc_nanodomain_size`, `calc_loc_density`, `calculate_nanodomain_characteristics` | `hmt_v3.cluster.calculate_nanodomain_characteristics` |
| `plot_clusters_napari`, `plot_seeds_napari` | `hmt_v3.visualize.*` (same names) |
| `remove_nucleoli` | none (it was an empty stub) |
| `build_probability_map`, `plant_seeds`, `plant_seeds_POISSON` | `hmt_v3.simulate.place_seeds` / `place_seeds_thomas`, or `generate_nucleus` |
| `extract_empirical_parameters_NEW` / `_OLD` | `hmt_v3.measure.extract_empirical_parameters` |
| `spawn_nanodomains_NEW` / `_OLD` | `hmt_v3.simulate.spawn_nanodomains` |
| `epsilon_cost`, `epsilon_cost_OLD` | `hmt_v3.calibrate.domain_size_cost` |
| `optimize_epsilon_hybrid`, `optimize_epsilon_grid_OLD` | `hmt_v3.calibrate.optimize_dbscan_size_matching` (also searches min_samples) |
| `calculate_geometric_gt_target`, `calculate_simulated_gt_radial` | `hmt_v3.measure.calculate_geometric_gt_target` |
