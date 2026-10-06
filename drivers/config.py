"""
Configuration for the sSMLM histone-modification analysis drivers.

Edit THIS file (not the step scripts) to point at data or change a parameter.
Every step imports it, so a value is set once and read everywhere.

Run order (from the repo root):
    python drivers/01_preprocess.py   # mask nucleus, radial contours           (fast)
    python drivers/02_calibrate.py    # inverse-MC DBSCAN eps/min_samples curves  (slow)
    python drivers/03_cluster.py      # density-corrected DBSCAN -> nanodomains   (fast)
    python drivers/04_spatial.py      # radial profiles + SNCR per cell           (fast)
    python drivers/05_compare.py      # pool cells, condition stats, figures      (fast)
or all of them:  python drivers/run_all.py
Every step accepts --quick (smoke-test settings below), --cells ID [ID ...] to
restrict to some cells, and --overwrite to redo cells whose outputs already exist.
"""
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]

# --- Data --------------------------------------------------------------------
# Folder of RainbowSTORM exports, already split by spectral classification into
# one CSV per mark: <cell_id>_me3.csv and <cell_id>_ac.csv. Columns used: x [nm],
# y [nm], z [nm]. Wherever both marks are needed together (nucleus mask, mask
# origin) the pipeline concatenates the me3 and ac tables itself; the
# RainbowSTORM *_comb.csv exports are never read and don't need to be present.
DATA_DIR = REPO_ROOT / "raw_data"
OUTPUT_DIR = REPO_ROOT / "analysis_output"

ME3_SUFFIX = "_me3.csv"
AC_SUFFIX = "_ac.csv"

# Condition is assigned from the cell_id: the FIRST (name, regex) whose regex is
# found in the cell_id wins (re.search, case-insensitive). The catch-all "" must be last.
CONDITION_PATTERNS = [
    ("ActD", r"act"),
    ("NT", r""),
]
CONDITION_ORDER = ["NT", "ActD"]  # plotting order; the first entry is the reference group

EXCLUDE_CELLS = []  # cell_ids to skip entirely (e.g. failed acquisitions)

CHANNELS = ("me3", "ac")
CHANNEL_LABELS = {"me3": "H3K27me3", "ac": "H3K27ac"}
CHANNEL_COLORS = {"me3": "#77DD76", "ac": "#FF6962"}

# --- 01 preprocess: nucleus mask + radial contours -----------------------------
AXIAL_PADDING_NM = 1      # filter_axial: drop locs piled up at the z-range edges
MASK_BIN_NM = 50          # histogram pixel for the mask / contours (also the contour pixel)
MASK_SIGMA = 4.0          # gaussian smoothing of the density map, in pixels
# Smoothed locs-per-pixel threshold for "inside the nucleus". "auto" =
# MASK_AUTO_FACTOR x Otsu of the smoothed map (preprocess.auto_mask_threshold),
# which adapts to each cell's labelling density; a fixed number (the old
# notebooks used 2) only fits cells of similar density -- on sparse ActD cells
# thresh=2 kept just the heterochromatin rim and dropped 45-80% of localizations.
MASK_THRESH = "auto"
MASK_AUTO_FACTOR = 0.5
MASK_THRESH_OVERRIDES = {}  # {cell_id: number} for cells the default masks badly (check the QC pngs)
MASK_AREA_WARN_UM2 = (30.0, 600.0)  # warn if the mask area falls outside this range
# create_radial_contours hardcodes "100 - band", so it only numbers bands
# correctly (0 = core ... 99 = periphery) for exactly 100 bands.
N_CONTOUR_BANDS = 100

# --- 02 calibrate: inverse Monte Carlo DBSCAN calibration (Fig. 4) --------------
# "per_condition": calibrate once per condition on its densest cell (Results text)
# "per_cell":      calibrate every cell on itself (slowest; Results "for each imaging experiment")
# "pooled":        one calibration from the densest cell overall
CALIBRATION_MODE = "per_condition"

USE_Z = True  # cluster in 3D (XYZ); must be the same in calibration and clustering

# generate_nucleus knobs. domain_scale is the within-domain cutoff radius (nm);
# None lets generate_nucleus read it off g(r), which tends to underestimate.
DOMAIN_SCALE_NM = {"me3": 200.0, "ac": 250.0}
NOISE_FRACTION = {"me3": None, "ac": None}  # None = estimate from g(r)
SEED_PLACEMENT = "thomas"  # "thomas" | "inhibition" | "poisson"

# The calibration runs on a square crop of the simulated nucleus (DBSCAN over a
# whole nucleus x every eps/min_samples x every fraction is too slow). The crop
# is placed where its density matches the simulated nucleus's MEAN density, and
# its area is the in-mask area inside it -- so its locs/nm^2 is measured the same
# way as a real cell's (locs / mask area) in step 03.
CROP_SIZE_NM = 5000.0
CROP_MIN_INSIDE = 0.95  # minimum fraction of the crop that must lie inside the mask

DENSITY_FRACTIONS = np.round(np.linspace(1.0, 0.1, 10), 3)  # manuscript: 1.0 .. 0.1 in 0.1 steps
MIN_SAMPLES_RANGE = range(5, 75, 5)
SIZE_COST_THRESHOLD_NM = 50.0  # fractions whose best size error exceeds this are left out of the fit
SIZE_PERCENTILE = 75.0         # diameter statistic used for calibration and reporting

RANDOM_SEED = 0

# --- 03 cluster ---------------------------------------------------------------
MIN_EPS_NM = 10.0
MIN_MIN_SAMPLES = 3
# A cell is flagged "extrapolated" if its density is outside the calibrated
# density range by more than this fraction. The calibration source cell itself
# sits a few % above the top fraction because the simulation reproduces its
# localization count only approximately.
EXTRAPOLATION_TOLERANCE = 0.05

# --- 04 spatial: radial profiles + SNCR (Figs. 3, 6) ----------------------------
SNCR_RADII_NM = np.arange(10, 310, 10)
# What counts as a neighbour when computing SNCR around a nanodomain:
#   "nanodomains":   other nanodomain centroids (Fig. 3d/e illustration)
#   "localizations": localizations, including the reference domain's own (Methods text, Fig. 3f)
SNCR_NEIGHBORS = ("nanodomains", "localizations")  # both are computed
# The one used for the condition figures / decay-rate stats. "localizations"
# reproduces the shape of Fig. 3f (rise to a peak at tens of nm, then decay);
# with the current calibration's large domains, "nanodomains" curves only rise
# over 0-300 nm (neighbouring centroids are hundreds of nm apart), so no decay to fit.
SNCR_PRIMARY = "localizations"
SNCR_SUMMARY = "mean_ratio"    # "mean_ratio" (per-domain average, manuscript) or "pooled_ratio"
SNCR_USE_Z = False             # 2D (xy) disks, as in the manuscript figures
SNCR_FIT_MAX_NM = 200.0        # exponential decay fit range (Fig. S3b/c)

# --- quick mode (--quick): small settings to smoke-test the whole chain -------
QUICK_OVERRIDES = {
    "OUTPUT_DIR": REPO_ROOT / "analysis_output_quick",
    "CROP_SIZE_NM": 3000.0,
    "DENSITY_FRACTIONS": np.round(np.linspace(1.0, 0.2, 5), 3),
    "MIN_SAMPLES_RANGE": range(5, 45, 10),
}
