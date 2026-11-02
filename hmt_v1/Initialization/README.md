# Initialization

Simulation-based DBSCAN calibration (steps 3–9 of the
[original instructions](../DOCUMENTATION.md#original-instructions)). It produces
`methyl_epsilon_surround2.mat` and `acetyl_epsilon_surround2.mat`, which hold the fit
ε = coef · density^exp that `Analysis/analyze_STORM.m` uses. For the overall data flow, missing
inputs and known issues, see [../DOCUMENTATION.md](../DOCUMENTATION.md).

Run the scripts from this folder, in this order:

```
denstest_1 → simulation_optimization_2 → sim_densities_3 → split_script_4 → param_opt2_5
```

Units: positions are in nm. Images and masks use 10 nm pixels (`px = 10`), so a mask area
`sa` is in µm² and densities are in localizations / µm².

## Scripts

### `denstest_1.m`
Makes 10 thinned copies of one real cell by keeping every (11 − k)-th localization: `dens1`
keeps every 10th row (10 %) and `dens10` keeps every row (100 %).

- **Edit:** `fnam1`, `fnam0` (lines 8–9), the me3 and ac file names without `.csv`.
- **Reads:** `<fnam1>.csv`, `<fnam0>.csv` (header row skipped).
- **Writes:** `dens1_me3.csv … dens10_ac.csv` with the 17-column RainbowSTORM header.
- The header comment ("match first order localizations…") is left over from another script.
  No other script reads the `dens*` files. `pxsz` is unused.

### `simulation_optimization_2.m`
Fits the nanodomain simulation to one real cell and derives the ground-truth domain size.

1. **Mask:** histograms both channels on a 10 nm grid, then blurs (σ = 1 px), thresholds at
   `thresh`, blurs again (σ = `dil` px), fills holes and keeps objects > 10⁶ px. `find_area` then
   refines the mask; `sa` is its area (µm²).
2. **Axial profile:** histograms every z within ±500 nm of the most common z in 10 nm bins.
   `N` is the fraction per bin and `frame` the bin centres (nm). These become the z-slices of
   the simulated nucleus.
3. **Real statistics:** `find_within2` → per-ring counts and |Δz| around me3 and ac domain
   centroids.
4. **Fit loop:** simulate (`sim_histones_3d`), split (`split_me_ac`), measure the simulation
   (`find_within`), update the per-ring correction factors (`create_variables`). It repeats up to
   10 times until every factor is within 30 % of 1 (`dif ≤ 0.3`).
5. **Ground truth:** `find_gt_width` gives the me3 and ac domain diameters.

- **Edit:** `simnames` (cell file prefix, line 2), `thresh` (mask threshold, line 3), `dil`
  (blur after thresholding, line 4).
- **Reads:** `<cell>_me3.csv`, `<cell>_ac.csv`; `initial_variables.mat`; everything
  `sim_histones_3d` loads.
- **Writes:** `Simulation.tif`, `Simulation.csv`, `Simulation_me3/_ac.csv` (overwritten every
  iteration), `sim_param2.mat` (via `sim_histones_3d`), `sim_param.mat` (`nmeth, nac, ratiom,
  ratioa, ratiomz, ratioaz, nseedm, nseeda, sa, N, frame, meth_between, ac_between_ac`) and
  `gt_values.mat` (`ac_width, m_width, dens, densa, thresh, dil`).
- **Note:** this copy reads z from **column 5**; `find_area`, the demo copy and the original instructions use
  column 12.

### `sim_densities_3.m`
Loads `sim_param.mat` and simulates 10 density levels with `sim_histones_3d_2`. Level n
multiplies the per-ring correction factors by n, and the simulator divides per-ring counts by
those factors, so `simdens1` is the densest and `simdens10` the sparsest.

- **Writes:** `simdens1surround2.csv … simdens10surround2.csv` (6-column simulation format).
- It builds a fresh nucleus with `sim_nucleus_3d`, but `sim_histones_3d_2`'s
  `load('sim_param2')` replaces it with the saved one. In this copy the `load` also undoes the
  ×n scaling, so the 10 levels may all have the same density. The demo copy first strips those
  variables from `sim_param2.mat`. See Known issue 2 in the main documentation.

### `split_script_4.m`
Splits each `simdens*surround2.csv` by spectral centroid (column 6): below 680 nm → `_me3.csv`,
above 680 nm → `_ac.csv`. Simulated centroids are 670 nm (me3) and 700 nm (ac).

### `param_opt2_5.m`
Finds the density-dependent DBSCAN ε.

1. Builds the nucleus mask from `simdens1` (threshold 0.17, objects > 10⁵ px); `sa` = its area.
2. For the simulated densities 1, 2, 5, 9 and 10, it searches for the ε (`factor`, nm; 2D
   DBSCAN, minPts 8) whose mean domain diameter (`epsilon_cost`) is within 1.3 nm of the
   ground truth. The search starts at 30 nm (later levels start from the previous result −
   2 nm), steps 5 nm, halves the step every time it changes direction, and stops after 20
   iterations.
3. Density = simulated localizations / `sa`. It fits log₁₀ ε against log₁₀ density, separately
   for me3 and ac.

- **Reads:** `simdens*surround2_me3/_ac.csv`, `gt_values.mat` (`m_width`, `ac_width`).
- **Writes:** `methyl_epsilon_surround2.mat` (`coef`, `exp`), `acetyl_epsilon_surround2.mat`
  (`coefa`, `expa`). Copy both to `Analysis/`.
- Also computes SNCR decay rates for all 10 levels with `calc_decay` (`decaym`, `decaya`); these
  stay in the workspace and are not saved. The me3/ac arguments are swapped at line 141.
- The masking block (lines 38–53) has no effect, since the loop reloads unmasked data
  (`exclude_bg` is commented out). The variable `exp` shadows MATLAB's `exp`.

## Functions

### `find_area(px, title, thresh, dil, cell_img2bin) → [cell_img2bin, sa]`
Re-reads `<title{1}>_me3/_ac.csv` (z from column 12) and rebuilds the summed density image. It
thresholds the image at `thresh`, multiplies by the incoming mask, blurs (σ = `dil`), fills holes
and keeps objects > 10⁵ px. Returns the refined mask and its area `sa` (µm²), and displays the
mask.

### `find_within(xm, ym, zm, xa, ya, za, px) → [meth_withinc, mdistz, ac_within_ac_c, adistz]`
Within-domain statistics of a cell (used on the simulation). It runs DBSCAN in 2D on each
channel (ε = 30 nm, minPts = 8) and takes cluster centroids. Then, for every centroid and ring
r = 10, 20, … 1000 nm:

- `meth_withinc` / `ac_within_ac_c`: cumulative count of same-channel localizations within r
  (2D), one row per domain
- `mdistz` / `adistz`: mean |Δz| from the centroid of the localizations in ring r

Plots each channel's clusters. Cross-channel counts are computed but not returned.

### `find_within2(xm, ym, zm, xa, ya, za, px, cell_img2bin, xm2, ym2, xa2, ya2)`
Same outputs as `find_within`, after first dropping localizations outside `cell_img2bin`
(`xm2` … are the pixel coordinates used for the lookup). **Bug:** x and y are dropped but z is
not, so z is misaligned after the first dropped point.

### `create_variables(meth_withinc, mdistz, ac_within_ac_c, adistz, sim_*, ratiom_i, ratioa_i, ratiomz_i, ratioaz_i) → [ratiom_fin, ratioa_fin, ratiomz_fin, ratioaz_fin, meth_between, ac_between_ac, dif]`
Compares the simulation with the real cell over the first 20 rings (200 nm). It turns
cumulative counts into per-ring counts (`*_between`, negatives set to 0) and takes the mean per
ring. Then it forms simulation/real ratios for counts (`ratiom`, `ratioa`) and for |Δz|
(`ratiomz`, `ratioaz`).

- The new factors are the old factors × these ratios. The simulators divide by the factors, so
  a simulation that is too dense or too spread in z gets pulled back.
- `dif` = max |ratio − 1| is the loop's convergence test.
- Also returns the real per-ring count tables (`meth_between`, `ac_between_ac`) that the
  simulators sample from.

### `find_gt_width(meth_between, nmeth, ac_between_ac, nac) → [a_width, m_width, dens, densa]`
Ground-truth domain diameter from the real per-ring counts. For rings to 200 nm, it takes the
histogram-weighted mean count and divides by ring area to get a radial density profile, then
subtracts the profile's minimum. Width = 2 × the first radius at which density falls below 20 %
of its maximum (nm). `dens`, `densa` are the profiles. `nmeth` and `nac` cancel out.

### `sim_nucleus_3d(nNucleoli) → nucleus`
A random synthetic nucleus: a 1400 × 1400 × 100 binary array with 10 nm pixels (a 14 µm field,
1 µm deep).

- **Outline:** a 250-vertex polygon whose radius follows a random walk, stretched by random x/y
  factors (aspect 0.5–2) and scaled to fit. A second random outline is drawn and the outline is
  interpolated linearly between the two through z, then filled and dilated by 5 px.
- **Nucleoli:** `nNucleoli` random 75-vertex polygons at random positions, redrawn until they lie
  inside the nucleus on the first and last slice, then subtracted.

### `sim_histones_3d(name, nmeth, nac, ratcool, ratcoola, ratcoolz, ratcoolaz, nseedm, nseeda, sa, N, frame, meth_between, ac_between_ac)`
Simulates one two-mark nucleus and writes `<name>.csv`.

1. **Nucleus:** `sim_nucleus_3d(2)`, cut to `length(N)` slices and saved as `<name>.tif`.
2. **Weight maps per slice:**
   - 50 nm rings are peeled inward from the edge (`imerode`, 5 px disk). Each ring is weighted
     by `mean_ratio^1.5` (me3:ac ratio, `mean_ratio.mat`) and by density `mdens`
     (`mean_mdens.mat`).
   - 50 nm rings are grown outward from the nucleoli; the inner quarter is weighted by
     `ratio/2.5 + 0.45` (`ratio.mat`) and `mdens` (`mdens2.mat`).
   - Each pixel takes the larger of the two.
   - The whole workspace is then saved to `sim_param2.mat`.
3. **Seeds:** each nucleus pixel draws one random number. Below `meth_thresh` it becomes an me3
   domain seed; above `ac_thresh` an ac seed. The thresholds scale with `nseedm` / `nseeda`, the
   local weights, `N(slice)` and the area ratio `sa / sa_sim`, using hardcoded normalisation
   constants 561 000 000 (me3) and 425 110 000 (ac). Seed z = `frame(slice)`.
4. **Domains:** for each seed and each 10 nm ring to 200 nm, it draws a count from the real
   per-ring distribution (`meth_between` / `ac_between_ac`) and divides it by 6 × `ratcool(ring)`.
   That many points are placed at random angles within the ring, with z offset
   ±2.8 · r / √2 / `ratcoolz(ring)` (r jittered ± 10 nm).
5. **Output:** `<name>.csv` with id, frame (= 1), x, y, z (nm) and centroid (670 me3, 700 ac).

Loads `mean_ratio.mat`, `mean_adens.mat` (unused), `mean_mdens.mat`, `ratio.mat` and
`mdens2.mat`. `nmeth` and `nac` only appear in a normalisation that cancels.

### `sim_histones_3d_2(name, …, ac_between_ac, nucleus)`
Steps 3–5 of `sim_histones_3d` only. It reuses the weight maps from `sim_param2.mat`, so every
density level shares one nucleus geometry. Because `load('sim_param2')` restores every variable
saved there, it also overwrites the function's arguments. This copy re-sets the output name
afterwards; see Known issue 2 for the rest.

### `split_me_ac()`
Function form of `split_script_4.m`, hardcoded to `Simulation.csv` → `Simulation_me3.csv` /
`Simulation_ac.csv`.

### `epsilon_cost(factor, xm, ym, zm, gt) → cost`
DBSCAN on (x, y) with ε = `factor` nm and minPts = 8 (coordinates are divided by px = 10
internally). Cost = mean over clusters of the mean x/y extent, minus `gt`, in nm. Returns −100
(and prints it) if DBSCAN finds no clusters. `zm` is unused.

### `exclude_bg(xm, ym, zm, cell_img2bin, minx, miny, px) → [xm, ym, zm]`
Drops localizations outside the mask, keeping x, y and z aligned. Only referenced in
commented-out lines of `param_opt2_5.m`.

### `calc_decay(xm, ym, zm, xa, ya, za) → [meth_pos, ac_pos]`
3D self:nonself contact ratio (SNCR) curves. For every localization it counts same- and
other-channel neighbours in spheres of r = 10 … 1000 nm. `meth_ac_con` = mean me3 neighbours /
mean ac neighbours around me3 localizations (self excluded), and `ac_meth_con` the reverse.
`compile_function_in` then fits the decay rates.

### `compile_function_in(ac_meth_con, meth_ac_con) → [meth_pos1, ac_pos1]`
Min-max normalises the first 20 radii (10–200 nm) of each SNCR curve. It skips the rising part
up to the peak, keeps the decreasing run after it, and fits `a·exp(b·r)` (`fit … 'exp1'`),
returning b (1/nm). A first value > 10⁶ (division by zero at small r) is dropped first. Prints
`IND` while searching. The ac window ends one point later than the me3 window.
