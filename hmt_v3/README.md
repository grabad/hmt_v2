# hmt_v3

Python toolbox for quantifying the spatial distribution of histone post-translational
modifications (H3K27me3 and H3K27ac) from two-colour 3D spectroscopic SMLM (sSMLM)
localizations. It masks the nucleus, calibrates DBSCAN against simulated ground truth so
that nanodomain sizes don't depend on labelling density, clusters each mark into
nanodomains, and measures how the two marks are arranged relative to each other.

To run the analysis, use the notebooks in the repo root: `full_pipeline.ipynb` (every cell,
both marks, via [`drivers/`](../drivers/README.md)) and `single_channel_nanodomains.ipynb`
(one channel of one cell). This page covers the library itself.

## Setup

```bash
conda env create -f environment.yml    # from the repo root
conda activate hmt_smlm
```

Import the package with the repo root as the working directory, or on `sys.path`:

```python
import hmt_v3 as hmt
hmt.preprocess.filter_axial(df)
```

`import hmt_v3` loads every module, but not napari. napari is imported only when you open
a viewer, so the package also works headless.

## Modules

The modules are listed in the order the analysis uses them.

| Module | Purpose |
|---|---|
| [`preprocess`](preprocess.py) | Axial / quantile filters, nucleus mask, radial contour bands, density subsampling |
| [`measure`](measure.py) | Density-independent statistics of one channel: radial profile, RDF/ADF, pair correlation g(r) / Ripley's K, noise fraction, domain diameter |
| [`simulate`](simulate.py) | Synthetic nuclei with known nanodomains, built from `measure` (`generate_nucleus`) |
| [`calibrate`](calibrate.py) | Fit DBSCAN `eps(density)` and `min_samples(density)` on simulated ground truth (Fig. 4) |
| [`cluster`](cluster.py) | Density-corrected DBSCAN, per-nanodomain size, count and density (Figs. 2, 6) |
| [`postprocess`](postprocess.py) | Two-mark spatial metrics: SNCR, cross g(r), Ripley's cross-K, Delaunay / MST / mesh mixing (Figs. 3, 6, S3) |
| [`sweep`](sweep.py) | Toy two-mark centroid model, used to validate the `postprocess` metrics |
| [`export`](export.py) | Write simulated localizations as ThunderSTORM CSVs |
| [`visualize`](visualize.py) | matplotlib figures and napari viewers |
| [`srev`](srev/__init__.py) | SR-EV chromatin polymer → labelled, blinking, rendered sSMLM movie simulator |

### How the pieces fit

The outline below shows how the modules connect. For every individual calculation, with its
inputs and code location, open [`pipeline_flow.html`](../pipeline_flow.html) in a browser.

```
real <cell>_me3.csv / <cell>_ac.csv
   │  preprocess.filter_axial ─► auto_mask_threshold ─► binarize_nucleus ─► create_radial_contours
   │                                                       mask_origin (x_min, y_min)
   ▼
masked me3 / ac ──► simulate.generate_nucleus  (uses measure.*: radial profile, RDF/ADF, g(r))
   │                     │  simulated nucleus with ground-truth cluster_label
   │                     ▼
   │               preprocess.sample_lower_densities ─► calibrate.optimize_dbscan_size_matching
   │                                                     eps(density), min_samples(density)
   ▼                                                              │
cluster.cluster_dbscan_density_corrected ◄────────────────────────┘
   │  cluster.calculate_nanodomain_characteristics   (size, count, centroids)
   ▼
postprocess.nanodomain_sncr ─► postprocess.fit_exponential_decay      (spatial relationship)
measure.extract_radial_density_profile                                (radial distribution)
```

## Quick start (one cell)

This mirrors `drivers/01`–`04`, with the driver defaults.

```python
import numpy as np
import pandas as pd
import hmt_v3 as hmt

me3 = hmt.preprocess.filter_axial(pd.read_csv("cell_me3.csv"))
ac  = hmt.preprocess.filter_axial(pd.read_csv("cell_ac.csv"))

# 1. nucleus mask + radial contours
thresh = hmt.preprocess.auto_mask_threshold(me3, ac)
mask, me3_m, ac_m = hmt.preprocess.binarize_nucleus(me3, ac, thresh)
x_min, y_min = hmt.preprocess.mask_origin(me3, ac)          # PRE-mask inputs!
_, bands = hmt.preprocess.create_radial_contours(mask, num_bands=100)
area_nm2 = mask.sum() * 50.0 ** 2

# 2. calibrate DBSCAN on a simulated copy of this cell (slow; drivers/02 crops first)
np.random.seed(0)                                             # see "Randomness" below
_, sim, info = hmt.simulate.generate_nucleus(me3_m, bands, x_min, y_min, domain_scale=200.0,
                                             rng=np.random.default_rng(0))
fractions = [sim] + [sim.sample(frac=f, random_state=i) for i, f in enumerate(np.linspace(0.9, 0.1, 9))]
cal = hmt.calibrate.optimize_dbscan_size_matching(fractions, area_nm2, min_samples_range=range(5, 75, 5),
                                                  use_z=True)

# 3. density-corrected clustering + nanodomain sizes
clustered, params = hmt.cluster.cluster_dbscan_density_corrected(me3_m, cal, area_nm2, use_z=True)
domains = hmt.cluster.calculate_nanodomain_characteristics(clustered)

# 4. spatial relationship (do the same clustering for ac first)
# sncr = hmt.postprocess.nanodomain_sncr(me3_centroids, me3_locs, ac_locs, np.arange(10, 310, 10))
```

## Conventions

- **Columns.** Localization tables use the RainbowSTORM / ThunderSTORM names `x [nm]`,
  `y [nm]`, `z [nm]`. Simulated and clustered tables add `cluster_label` (-1 = noise /
  background); `cluster.cluster_dbscan` also adds `cluster_color`.
- **Mask grid.** `binarize_nucleus` bins both marks on a `bin_size` (default 50 nm) grid whose
  origin is the minimum x / y of the *pre-mask* inputs. Every function that maps points onto
  the mask takes that origin as `x_min, y_min`, so always get it from
  `preprocess.mask_origin(me3_df, ac_df)` on the same inputs you gave `binarize_nucleus`.
  If you use the masked outputs' minimum instead, points are shifted into the wrong bands.
  Mask arrays are indexed `[x_bin, y_bin]`, so transpose before `imshow`.
- **Contour bands.** `0` = innermost core … `num_bands - 1` = outer periphery, `-1` = outside.
- **Density** is localizations per nm² of *mask* area (not bounding box), the same in
  calibration and in clustering.
- **Domain size.** `geometric_diameter_2D` is 2 × the 75th percentile of the xy distance from
  points to the centroid. The calibration target (`measure.calculate_geometric_gt_target`), the
  calibration cost (`calibrate.domain_size_cost`) and the reported nanodomain size all use this
  same statistic. `bb_size_2D` (mean bounding-box extent) is the size quoted in the manuscript.
- **Randomness.** Functions that take `rng` expect a `np.random.Generator`. Some of the
  simulation steps (`spawn_nanodomains`, `add_noise_locs`, `add_sim_noise`, `place_seeds`,
  `make_grid_seeds`) draw from NumPy's *global* random state instead. As a result,
  `generate_nucleus(rng=...)` alone does not give a reproducible nucleus: call
  `np.random.seed(...)` as well. `drivers/02_calibrate.py` and `calibrate.sweep_domain_scale` do
  this for you.
- **Plots** call `plt.show()`; under a non-interactive backend the figures stay open, which is
  how `drivers/common.capture_figures` saves them.

## Function reference

Every function has a docstring with `Args:` / `Returns:` sections; use `help(...)` or
`hmt.<module>.<function>?` in Jupyter for full details.

### preprocess

| Function | What it does |
|---|---|
| `filter_quantile` | Keep rows inside a column's [p, 1-p] quantile range |
| `filter_axial` | Drop localizations clamped to the edges of the axial fit range |
| `mask_origin` | (x_min, y_min) of the mask grid: the minimum of the pre-mask inputs |
| `auto_mask_threshold` | Density-adaptive mask threshold: 0.5 × Otsu of the smoothed density map |
| `binarize_nucleus` | Threshold the combined, smoothed density into a filled nucleus mask; return the masked points |
| `create_radial_contours` | Split the mask into concentric distance bands |
| `assign_contour_bands` | Add each localization's contour band as a column |
| `sample_lower_densities` | Thin both marks to a ladder of labelling densities |

### measure

| Function | What it does |
|---|---|
| `extract_radial_density_profile` | Localizations per pixel in each contour band (periphery vs. interior enrichment) |
| `extract_empirical_parameters` | Within-domain radial (RDF) and axial (ADF) neighbour distributions |
| `extract_n_locs_from_rdf` | Expected localizations per domain from the RDF integral |
| `estimate_noise_fraction` | Noise fraction from a local-density threshold (heuristic) |
| `measure_pair_correlation` | Edge-corrected Ripley's K, L(r) and g(r) of one channel inside the mask |
| `estimate_n_domains_from_g` | Domain count and scale from g(r), assuming a Poisson-cluster model |
| `estimate_noise_fraction_from_pcf` | Noise fraction from the g(r) excess relative to a noise-free reference |
| `fit_thomas_parameters` | Thomas-process parent intensity and spread from the inter-domain g(r) |
| `calculate_geometric_gt_target` | Mean ground-truth domain diameter (2 × p75 radius) |
| `measure_seed_spacing` | DBSCAN domain centres and nearest-neighbour spacing (density-dependent, kept for reference) |

### simulate

| Function | What it does |
|---|---|
| `generate_nucleus` | **Main entry.** Simulate a nucleus matching a real cell's radial profile, domains, g(r) and noise |
| `simulate_single_domain_scene` | One isolated domain in uniform noise, with an exactly known diameter |
| `build_z_band_pools` | Real z values grouped by contour band, so seed z can follow radial position |
| `place_seeds` | Radially weighted Poisson seed placement |
| `place_seeds_matched` | Radially weighted placement with inhibition (hard-core or empirical spacing) |
| `place_seeds_thomas` | Radially weighted Thomas cluster process (domains of domains) |
| `make_grid_seeds` | Regular grid of seeds, for sensitivity tests |
| `spawn_nanodomains` | Localizations around each seed, sampled from the RDF / ADF |
| `add_noise_locs` | Background localizations inside the mask (uniform or band-weighted) |
| `add_sim_noise` | Uniform noise in a simulation's bounding box (no mask needed) |

### calibrate

| Function | What it does |
|---|---|
| `optimize_dbscan_size_matching` | **Used by the drivers.** Best (eps, min_samples) per density by size matching, then fit eps(density) and min_samples(density) |
| `domain_size_cost` | \|mean DBSCAN cluster diameter − ground-truth diameter\|, with a penalty for merged "continents" |
| `diagnose_domain_packing` | Fraction of ground-truth domains that physically overlap their nearest neighbour |
| `sweep_domain_scale` | Try several `generate_nucleus` domain scales and check separability vs. DBSCAN recoverability |
| `optimize_dbscan_params` | Earlier ARI-driven joint (eps, min_samples) search, with a domain-F1 comparison |
| `optimize_epsilon` | Earlier ARI-driven eps-only search |
| `epsilon_cost_cluster` | 1 − ARI of a DBSCAN run against ground truth |
| `domain_detection_f1` | Per-domain detection precision / recall / F1 at an IoU threshold |

### cluster

| Function | What it does |
|---|---|
| `cluster_dbscan_density_corrected` | DBSCAN with eps / min_samples evaluated from a calibration at this cell's density |
| `cluster_dbscan` | DBSCAN with fixed parameters; adds `cluster_label` and `cluster_color` |
| `calculate_nanodomain_characteristics` | One row per nanodomain: count, centroid, hull / bounding-box / geometric size, density |

### postprocess

| Function | What it does |
|---|---|
| `nanodomain_sncr` | Self:nonself contact ratio around each nanodomain of a real nucleus (Fig. 3d–f) |
| `fit_exponential_decay` | Decay rate k of an SNCR curve, fitted from its peak (Fig. S3) |
| `colocalization_fraction` | Fraction of A points with a B point within r |
| `cross_nn_distances` | Nearest-B distance for each A point |
| `cross_pair_correlation` | Cross g_ab(r) in a square field (border exclusion) |
| `self_nonself_contact_ratio` | Raw and abundance-normalised self/nonself ratio vs. r in a square field |
| `ripleys_k_cross` | Bivariate Ripley's K / L in both directions |
| `delaunay_channel_mixing` | Heterotypic edge fraction and assortativity of the pooled Delaunay graph |
| `friedman_rafsky_test` | MST two-sample test (cross edges vs. a label-permutation null) |
| `mst_merge_curve` | Cumulative heterotypic merge fraction vs. merge radius (H0 persistence) |
| `mesh_overlap` | How often the two marks' own Delaunay meshes cross, with an optional permutation null |

### sweep

| Name | What it does |
|---|---|
| `simulate_toy_centroids` | me3 / ac centroids with a tunable fraction of jittered pairs (symmetric integration) |
| `generate_integration_sweep` | `simulate_toy_centroids` at several integration levels from one RNG |
| `simulate_pair_centroids` | The same model, parametrised by total count + me3 fraction |
| `simulate_one_sided_centroids` | Asymmetric model: one channel clusters onto an independent other |
| `summarize_metrics` | Every `postprocess` metric on one point pair, as a flat dict (`METRIC_COLUMNS`) |
| `oat_parameter_sweep` | One-at-a-time sweep of the five toy parameters (`BASELINE`, `DEFAULT_SWEEP_GRID`) |
| `one_sided_clustering_sweep` | Sweep of the asymmetric model's cluster fraction |
| `density_subsample_experiment` | Metric drift under random thinning of one fixed dataset |
| `perturbation_experiment` | Metric stability under Gaussian position jitter |
| `summarize_stability` | Per-group mean / std / CV of every metric |
| `plot_sncr_sweep`, `plot_ripley_k_sweep`, `plot_delaunay_sweep`, `plot_mst_sweep`, `plot_merge_curve_sweep`, `plot_mesh_overlap_sweep` | One metric across a `generate_integration_sweep` result |

### export

| Function | What it does |
|---|---|
| `export_to_thunderstorm` | Write a simulation as a ThunderSTORM CSV with z-matched uncertainty from a real channel |
| `sample_conditional_on_z` | Bootstrap a column from real data conditioned on z (k nearest in z) |

### visualize

| Function | What it does |
|---|---|
| `plot_clusters_napari` | napari 3D view of both marks' nanodomains and noise |
| `plot_seeds_napari` | napari 3D view of simulated seeds |
| `plot_nanodomain_2d` / `plot_nanodomain_3d` | Scatter of labelled localizations (and seeds, hulls) |
| `plot_rdf_adf` | me3 vs. ac RDF and ADF |
| `plot_pair_correlation` | Real vs. simulated L(r) − r and g(r) |
| `plot_seed_spacing` | Real vs. simulated inter-domain spacing histograms |
| `plot_domain_packing` | Histogram of overlap margins from `calibrate.diagnose_domain_packing` |
| `plot_dbscan_calibration` | me3 and ac eps / min_samples calibration curves |
| `plot_size_vs_density` | Measured single-domain size vs. density: fixed vs. corrected parameters (Fig. 4b/c) |

### srev

SR-EV simulation: chromatin polymer → antibody / dye labelling → blinking → rendered sSMLM
movie (single channel, or a dual-dye spatial + spectral pair for RainbowSTORM). The
orchestrators are `srev.pipeline.run_pipeline` and `srev.pipeline.run_multichannel_pipeline`;
`legacy/notebooks/sr_ev_pipeline_demo.ipynb` walks through both. Each `srev` module opens with
a docstring explaining its model and sources.

| Module | Contents |
|---|---|
| `io` | `parse_config_txt` (SR-EV `config-N.txt` → monomer table), `resolve_srev_run_params`, folder-name / `copy-parameters.txt` parsers |
| `compartments` | `compute_coordination_number`, `classify_compartments` (hetero- / euchromatin by coordination number) |
| `domains` | `identify_packing_domains` (density image → domain centres → radii), plus its steps |
| `labeling` | `LabelingConfig`, `label_monomers` (primary → secondary → fluorophore placement with steric rejection) |
| `photophysics` | `PhotophysicsConfig`, `simulate_blinking` (per-fluorophore Markov chain) |
| `render` | `CameraConfig`, `RenderConfig`, astigmatic PSF rendering, EMCCD noise, `render_stack_to_tiff` |
| `spectral` | `DyeSpectrum`, `SpectralConfig`, dispersed first-order rendering, `render_dual_stack_to_tiff`, calibration-lamp pair, centroid check |
| `ground_truth` | `resolved_localizations_gt`, `save_run` / `load_run`, `check_referential_integrity` |
| `pipeline` | `SrEvSimConfig`, `simulate_channel`, `run_pipeline`, `merge_channels`, `run_multichannel_pipeline` |
| `visualize` | Particle tables for OVITO (`export_for_ovito`), napari (`view_napari`) and matplotlib (`plot_structure_3d`) |

## Where things moved (October 2026 reorganisation)

`simulate.py` was split into `measure`, `simulate`, `calibrate` and `export`, and the toy
model moved from `postprocess` to `sweep`. Function names are unchanged; only the module
prefix differs. The drivers and every notebook in this repo (including `legacy/notebooks/`)
have been updated.

| Old | New |
|---|---|
| `simulate.extract_*`, `estimate_*`, `measure_*`, `fit_thomas_parameters`, `calculate_geometric_gt_target` | `measure.*` |
| `simulate.optimize_*`, `domain_size_cost`, `epsilon_cost_cluster`, `domain_detection_f1`, `diagnose_domain_packing`, `sweep_domain_scale` | `calibrate.*` |
| `simulate.export_to_thunderstorm`, `sample_conditional_on_z` | `export.*` |
| `simulate.assign_contour_bands` | `preprocess.assign_contour_bands` |
| `postprocess.simulate_toy_centroids`, `generate_integration_sweep`, `plot_*_sweep` | `sweep.*` |
| `cluster.calc_nanodomain_size`, `calc_loc_density` | removed (superseded by `calculate_nanodomain_characteristics`) |
| `hmt_v3.simulate_OLD`, root `hmt_functions.py` | [`legacy/`](../legacy/README.md) |
