# hmt_v2

Code for quantifying the spatial distribution of histone post-translational modifications
(H3K27me3 / H3K27ac) from two-colour 3D spectroscopic SMLM, and for simulating ground truth to
validate it.

## Start here

| Notebook | What it runs |
|---|---|
| [`full_pipeline.ipynb`](full_pipeline.ipynb) | The whole two-mark analysis on a folder of cells: preprocess, DBSCAN calibration, density-corrected clustering, radial profiles / SNCR, condition comparison (same steps as `drivers/run_all.py`, with QC shown after each one) |
| [`single_channel_nanodomains.ipynb`](single_channel_nanodomains.ipynb) | One channel of one cell: preprocess, inverse Monte Carlo calibration, size-vs-density check, density-corrected DBSCAN, nanodomain characteristics |

Open either notebook with the repo root as the working directory.

To see what every step calculates before running it, open
[`pipeline_flow.html`](pipeline_flow.html) in a browser: an interactive flowchart of the whole
pipeline (steps 01–05). Click any box to see how it is computed, what it depends on, and the
code that does it.

To see what the pipeline produces, read the [worked example](example/README.md): the full
analysis run on the five `test_data/` cells, following one nucleus through every step with
all of its tables and figures.

## Layout

| Path | What it is |
|---|---|
| [`hmt_v3/`](hmt_v3/README.md) | **The toolbox** (Python package): preprocessing, simulation, DBSCAN calibration, clustering, spatial metrics, SR-EV movie simulator. Its README has a function reference. |
| [`drivers/`](drivers/README.md) | Scripts behind `full_pipeline.ipynb`; also runnable from the command line (`python drivers/run_all.py`). All analysis parameters are in `drivers/config.py`. |
| [`pipeline_flow.html`](pipeline_flow.html) | Interactive flowchart of every calculation in the drivers' steps 01–05 (open in a browser) |
| [`example/`](example/README.md) | Worked example: the full pipeline on `test_data/`, with every output for one nucleus and the pooled results |
| `notebooks/` | Notebooks built on `hmt_v3`: [`manuscript_figures.ipynb`](notebooks/manuscript_figures.ipynb) draws the methods paper's four figures (data and parameters, simulation, calibration, clustering one nucleus) from a driver run; [`edu_s_phase_ripley.ipynb`](notebooks/edu_s_phase_ripley.ipynb) compares Ripley's K / L and g(r) of EdU replication patterns across S-phase stages (`meng_data/`) |
| [`legacy/`](legacy/README.md) | Archived code and notebooks: the pre-`hmt_v3` functions, earlier exploratory / figure / SR-EV notebooks, and the data they wrote |
| [`hmt_v1/`](hmt_v1/DOCUMENTATION.md) | Original MATLAB toolbox and SR-EV movie simulation (also the source of the SR-EV `config-1.txt` and dye blinking statistics). [`hmt_v1/DOCUMENTATION.md`](hmt_v1/DOCUMENTATION.md) has its original instructions, data flow, missing inputs and known issues |
| `reference_manuscripts/` | Papers the methods are based on |
| `environment.yml` | Conda environment (`conda env create -f environment.yml`) |

Data and output folders (`test_data/`, `raw_data/`, `meng_data/`, `sr-ev_raw/`, `analysis_output*/`) are
git-ignored.
