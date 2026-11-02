# Analysis

Per-cell analysis and pooling (steps 10–16 of the
[original instructions](../DOCUMENTATION.md#original-instructions)). It needs `methyl_epsilon_surround2.mat`
and `acetyl_epsilon_surround2.mat` from [`../Initialization`](../Initialization/README.md), copied
into this folder with the cells' CSVs. For the overall data flow, missing inputs and known
issues, see [../DOCUMENTATION.md](../DOCUMENTATION.md).

```
analyze_STORM  (once per cell, writes <cell>.mat)  →  compile_data  (pools the .mat files)
```

## `analyze_STORM.m` (script)

Runs the full analysis on each cell in `simnames` and saves `<cell>.mat`.

1. **Load** `<cell>_me3.csv` and `<cell>_ac.csv` (x column 3, y column 4, **z column 5 in this
   copy**; the demo copy reads `_C1`/`_C2` files with z in column 12).
2. **Localization-level contacts (3D):** for every localization it counts me3 and ac neighbours
   within spheres of r = 10, 20, … 1000 nm:
   - `meth_within`, `ac_within` around me3 localizations
   - `ac_within_ac`, `meth_within_ac` around ac localizations
   - the mean distance to the 10 nearest neighbours of each channel (`pix_*`)

   From these it computes the self:nonself contact ratios `meth_ac_con` and `ac_meth_con`
   (SNCR), plus "RDF" curves `meth_domain`, `ac_domain` = mean count / (4/3)π(r/32)³. The 32 is
   not explained in the code.
3. **Nucleus mask** on a 10 nm grid: blur σ = 1 px, threshold **0.05** (line 155), blur σ = 2,
   fill holes, keep objects > **100 000 px** (line 158). A figure shows the mask; adjust these
   two values if it is wrong.
4. **Radial profile:** 50 nm rings are peeled inward from the edge of the first mask object.
   For every ring holding > 85 me3 localizations it records the me3:ac ratio and both densities
   (localizations / µm²). These are `ratio`, `mdens`, `adens`, listed periphery first. The
   ring masks are written to `ring_binary1.tif`.
5. **Nanodomains:**
   - The sparsest 5 % of each channel (largest mean 10-NN distance) and all points outside the
     mask are dropped.
   - DBSCAN runs in 2D with minPts 8 and ε = `coef · density^exp` (me3) or
     `coefa · density^expa` (ac), from the two calibration files.
   - Density = all of that channel's localizations / mask area (µm²).
   - Each cluster gives a centroid.
6. **Domain-centred contacts:** neighbour counts and |Δz| per ring around each domain centroid
   (`meth_withinc`, `ac_withinc`, `ac_within_ac_c`, `meth_within_ac_c`, `mdistz`, `adistz`), and
   the centroid SNCR curves `meth_ac_con_c`, `ac_meth_con_c`.
7. **Domain size:** `diam_m`, `diam_a` = mean of each cluster's x and y extent, **in 10 nm pixel
   units** (`compile_data.m` multiplies by 10).
8. **Save** `<cell>.mat` with `cell_img2bin, sdis, meth_within, ac_within, ac_within_ac,
   meth_within_ac, meth_ac_con, ac_meth_con, meth_domain, ac_domain, meth_ac_con_c,
   ac_meth_con_c, meth_withinc, ac_withinc, ac_within_ac_c, meth_within_ac_c, ratio, diam_m,
   diam_a, mdens, adens, adistz, mdistz`.

**Edit:** `simnames` on **both** line 3 and line 6. The loop clears the workspace each
iteration, so line 6 re-creates the list.

**Caveats** (details in the main documentation):
- Points outside the mask are dropped from x and y but not z (lines 221–236), so domain-centroid
  z values are misaligned.
- In the centroid contact counts, same-channel counts use 2D discs and other-channel counts 3D
  spheres, and "−1 for self" is subtracted even though a centroid is not a localization.
- Every localization runs its own range searches, so this is slow for large cells.

## `compile_data.m` (script)

Pools the per-cell `.mat` files. For each cell it loads the file and calls `compile_function`,
keeping that cell's results under a numbered suffix (`…1`, `…2`, …) and appending them to
`sample.mat`. Then it pools:

- **Plots:**
  - mean normalised SNCR curves (10–200 nm) per channel
  - mean me3:ac ratio vs. radius (%)
  - mean me3 and ac density vs. radius (%)
- **Saves to `sample.mat`:**
  - `meth_pos`, `ac_pos`: per-cell SNCR decay rates
  - `meth_n`, `ac_n`: normalised SNCR curves
  - `diam_methylation`, `diam_acetylation`: mean domain diameter, nm
  - `mean_ratio`, `mdens`, `adens`: radial profiles
  - `abundance`: me3:ac localization ratio
  - `no_meth`, `no_ac`: localization counts
  - `nclusm`, `nclusa`: domain counts
- **Computed but not saved:** `mpeak` / `apeak` (peak SNCR ÷ mean SNCR at 800–1000 nm),
  `surroundme` / `surroundac`, `densitym` / `densitya`, `sa`.

**Edit:** the file name in every `load(...)` and `save(...)`, and add or remove blocks to match
the number of cells. As committed, all four blocks load `well_1_1.mat`, so the "pool" is one
cell four times.

## `Copy_of_compile_data.m` (script)

`compile_data.m` set up for two cells, `well_2_1.mat` and `well_2_2.mat`, saving to `well_2.mat`.

## `compile_function(ac_meth_con, ac_within_ac, diam_a, diam_m, meth_ac_con, meth_within, ratio, sdis, mdens, adens)` (function)

Summarises one cell's `.mat` contents. Returns 16 values:

| Output | Meaning |
|---|---|
| `surroundme1`, `surroundac1` | mean cumulative same-channel count curve ÷ number of localizations |
| `meth_pos1`, `ac_pos1` | SNCR decay rate b (1/nm) from an `exp1` fit to the falling part after the peak, within 10–200 nm (same procedure as `Initialization/compile_function_in.m`) |
| `meth_n1`, `ac_n1` | SNCR at 10–200 nm, min-max normalised |
| `meth_ab1` | me3 : ac localization count ratio |
| `ratio1`, `mdens1`, `adens1` | radial ratio and density profiles resampled to 100 points (radius %) |
| `diam_m1`, `diam_a1` | mean domain diameter (10 nm pixel units) |
| `no_meth1`, `no_ac1` | localization counts |
| `nclusm1`, `nclusa1` | domain counts |
