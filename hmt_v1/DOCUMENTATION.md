# hmt_v1 documentation

`hmt_v1` is the original MATLAB implementation of the toolbox from *"A Toolbox to Quantify the
Spatial Distribution of Post-translational Histone Modifications from 3D Spectroscopic SMLM
Images"*. The authors' original step-by-step instructions (formerly `README.md`) are reproduced
verbatim under [Original instructions](#original-instructions). This file also explains what each step does, which files flow between them, and what to watch out for.
Each folder also has its own README with a per-file reference.

These notes were written in October 2026 by reading the code; no `.m` file was modified. The
Python package [`hmt_v3`](../hmt_v3/README.md) re-implements and extends this toolbox (see
[Where each piece lives in hmt_v3](#where-each-piece-lives-in-hmt_v3)).

## Contents

| Folder | What it does | Reference |
|---|---|---|
| `Initialization/` | Simulation-based DBSCAN calibration: simulate nuclei that match a real cell, find the DBSCAN ε that recovers the true domain size at each labelling density, fit ε(density) | [Initialization/README.md](Initialization/README.md) |
| `Analysis/` | Per-cell analysis (nucleus mask, radial me3:ac profile, density-corrected DBSCAN nanodomains, self:nonself contact ratios), then pooling across cells | [Analysis/README.md](Analysis/README.md) |
| `Simulation/` | Separate pipeline: SR-EV chromatin polymer → antibody-labelled, blinking, two-dye sSMLM camera movies (zeroth + first order) with noise, for validating the toolbox on known ground truth | [Simulation/README.md](Simulation/README.md) |
| `demo/` | An alternate copy of `Initialization/` + `Analysis/`, the version that matches the original instructions' input format | [demo/README.md](demo/README.md) |

## Requirements

- **MATLAB R2019a or later**: `dbscan` and `writecell` were introduced in R2019a.
- **Toolboxes:**
  - Statistics and Machine Learning: `dbscan`, `rangesearch`, `knnsearch`, `kmeans`, `random`
  - Image Processing: `imbinarize`, `imgaussfilt`, `imfill`, `bwareaopen`, `bwconncomp`, `strel`, `imerode`, `imdilate`, `imtranslate`
  - Curve Fitting: `fit(..., 'exp1')`
  - Signal Processing: `resample`
- **Third-party functions:** `txt2mat` (MATLAB File Exchange, © Andres) is included in
  `Simulation/`. `saveastiff` (File Exchange "Multipage TIFF stack"), which `master2.m` calls, is
  **not** included.
- Every script reads and writes files in the **current folder**, so run each one from its own
  folder with its inputs copied next to it.

## Original instructions

The authors' original `README.md`, reproduced verbatim (only its title is demoted to a
sub-heading). File names and line numbers in it describe the code as it was then; the notes
after it say where the committed code differs.

### Histone-Modification-Toolbox
Contains all processing code from the paper "A Toolbox to Quantify the Spatial Distribution of Post-translational Histone Modifications from 3D Spectrscopic SMLM Images". This code is intended for use on multicolor images of nuclear proteins such as histone modifications or RNA polymerase
1.	Process all 3D multicolor SMLM datasets and save the localization lists as CSV files with the x-locations in the 3rd column, the y-locations in the 4th column, and the z-location in the 12th column.
2.	Split the results into the two color channels (1 and 2) with the channel 1 filename appended as “_C1” and the channel 2 filename appended as “_C2”.
3.	Find the largest file and place it and its “_C1” and “_C2” files in the “Initialization” folder.
4.	Open the MATLAB file “denstest_1.m” and In lines 8-9, change the “_C1” and “_C2” filenames to match the CSV files you placed in the initialization folder. Run the code.
5.	Open the file “simulation_optimization¬_2.m” and change line one to match the filename of original file you placed in the Initialization folder (i.e. not the “_C1” or “_C2” appended files). Run the code
6.	Open and run the code “sim_densities_3.m”
7.	Open and run “split_script_4.m”
8.	Open and run “param_opt2_5.m”
9.	Transfer the output files “methyl_epsilon_surround2” and “acetyl_epsilon_surround2” to the “Analysis” folder
10.	Place all 2 color images and “_C1” and “_C2” appended files for a single experimental condition in the analysis folder.
11.	Open “analyze_STORM.m”
12.	 In line 3 and line 6 change the list of names (saved to the variable “simnames”) to the names of the csv files.
13.	Run the code. View the resulting figures and ensure that the nuclear binarization looks correct.
14.	If the binarized nucleus includes too much signal, lower the threshold in line 155 (set to 0.17). If it includes too little, raise it.
15.	If no objects appear, lower the threshhold set in line 158 (set to 500000). If too many objects appear, raise it.
16.	Open the file “compile_data.m” and change each “load” and “save” filename to match the file names of the localization csv files. Adjust the number of variables in lines 128-170 accordingly.

**Where these steps no longer match the code**

- **Step 1:** RainbowSTORM exports do have z in column 12, but the main copies of
  `analyze_STORM.m` and `simulation_optimization_2.m` read z from column 5. The `demo/` copies
  read column 12 (see [Input file formats](#input-file-formats)).
- **Steps 2–4, 10:** the code expects `_me3` / `_ac` file suffixes, not `_C1` / `_C2`. The
  exception is `demo/Analysis/analyze_STORM.m`, which reads `_C1` / `_C2`. In step 4 the lines
  to edit are still 8–9 of `denstest_1.m` (`fnam1`, `fnam0`).
- **Step 5:** the script is `simulation_optimization_2.m` (the `¬` in the name above is a stray
  character). The cell name is on line 2, not line 1. The step also needs `.mat` inputs that
  must first be restored from git (see [Missing inputs](#missing-inputs)).
- **Steps 14–15:** in the main copy the mask threshold is 0.05 (line 155) and the minimum object
  size 100000 (line 158). The values quoted, 0.17 and 500000, are the `demo/` copy's, at its
  lines 181 and 184.
- **Step 16:** as committed, `compile_data.m` loads `well_1_1.mat` in all four blocks.

## Workflow 1: calibrate, then analyse cells (original steps 1–16)

### What the calibration does

DBSCAN with a fixed ε merges or splits nanodomains differently depending on how densely a cell
is labelled. The toolbox therefore calibrates ε as a function of localization density on
simulated cells whose true domain size is known:

1. Measure, around DBSCAN domain centroids of one real cell (ε = 30 nm, minPts = 8), how many
   localizations fall in each 10 nm ring and how far apart they are in z.
2. Simulate a nucleus with randomly seeded domains that reproduce those per-ring statistics,
   iterating until the simulation matches (up to 10 rounds, each per-ring factor within 30 %).
3. Take the ground-truth domain diameter from the real ring-density profile: 2 × the radius
   where density falls below 20 % of its peak.
4. Simulate 10 density levels. For five of them, find the ε at which DBSCAN's mean domain
   diameter equals the ground truth (within 1.3 nm).
5. Fit ε = coef · density^exp, with density in localizations / µm² of nucleus.

`analyze_STORM.m` then clusters each real cell with ε evaluated at that cell's own density.

### Steps and files

Run in this order. Inputs are copied into the folder by hand. Files marked *(in git history)*
were removed from the working tree but can be restored (see [Missing inputs](#missing-inputs)).

| # | Run | Reads | Writes |
|---|---|---|---|
| 1 | `Initialization/denstest_1.m` | `<cell>_me3.csv`, `<cell>_ac.csv` | `dens1…dens10_me3/_ac.csv` (10 %…100 % subsamples; not used later) |
| 2 | `Initialization/simulation_optimization_2.m` | `<cell>_me3/_ac.csv`, `initial_variables.mat`, and via `sim_histones_3d`: `mean_ratio.mat`, `mean_adens.mat`, `mean_mdens.mat`, `ratio.mat`, `mdens2.mat` *(all in git history)* | `Simulation.tif/.csv`, `Simulation_me3/_ac.csv`, `sim_param2.mat`, **`sim_param.mat`**, **`gt_values.mat`** |
| 3 | `Initialization/sim_densities_3.m` | `sim_param.mat`, `sim_param2.mat` | `simdens1surround2.csv … simdens10surround2.csv` |
| 4 | `Initialization/split_script_4.m` | `simdens*surround2.csv` | `simdens*surround2_me3.csv`, `_ac.csv` |
| 5 | `Initialization/param_opt2_5.m` | `simdens*surround2_me3/_ac.csv`, `gt_values.mat` | **`methyl_epsilon_surround2.mat`** (`coef`, `exp`), **`acetyl_epsilon_surround2.mat`** (`coefa`, `expa`) |
| 6 | copy the two `*_epsilon_surround2.mat` files to `Analysis/` | | |
| 7 | `Analysis/analyze_STORM.m` | `<cell>_me3/_ac.csv` for each cell, the two ε files | `<cell>.mat` per cell, `ring_binary1.tif` |
| 8 | `Analysis/compile_data.m` | the `<cell>.mat` files | `sample.mat` (pooled results) + figures |

An example ε fit ships in `demo/Analysis/acetyl_epsilon_surround2.mat`: `coefa = 51.1`,
`expa = -0.077`, i.e. ε ≈ 51 nm · density^-0.077.

## Workflow 2: SR-EV movie simulation (`Simulation/`)

```
SR-EV generator (external) ─► config-N.txt  (monomer x, y, z; only config-1.txt is included)
   ─► vis_SREV_int.m or vis_SREV_seg.m       label monomers euchromatin / heterochromatin
        ─► config_N_<int|seg>_alph*.mat  (+ .xyz for OVITO via gen_ovito_fin)
   ─► master2.m                              per configuration:
        sim_anti                  antibody labelling → fluorophore positions
        blinking_sim_unrealistic  zeroth-order movie + ground-truth table
        make_first_order          dispersed dye spectra
        make_real_fist            first-order (spectral) movie
        add_noise                 shot + background + read noise
        ─► <name>_z.tif, _f.tif (noise-free), <name>_zerothh.tif, _firsth.tif (noisy)
   ─► fit the movies externally (ThunderSTORM), then QD_analyze_A (missing) ─► <name>_fin.csv
   ─► splitmeac.m ─► <name>_fin_me3.csv / _fin_ac.csv  → input for Workflow 1
```

Heterochromatin is labelled with Alexa Fluor 647 (the "me3" channel) and euchromatin with CF680
(the "ac" channel).

## Input file formats

| Format | Columns | Used by |
|---|---|---|
| RainbowSTORM export (17 columns) | id, frame, x [nm], y [nm], sigmax, sigmay, intensity, offset, uncertainty, **centroid [nm]** (10), bgstd, **z [nm] (12)**, z uncertainty, ydif, centroid2, fwhm0, fwhm1 | real cells; `denstest_1.m` and `splitmeac.m` write this header |
| Simulation output (6 columns) | id, frame, x [nm], y [nm], **z [nm] (5)**, **centroid [nm] (6)** | everything written by `sim_histones_3d*.m` and the split scripts |

Every script reads x from column 3 and y from column 4. **Which column holds z differs between
copies of the same script.** The original instructions say column 12 (RainbowSTORM), but:

| Script | z column |
|---|---|
| `Analysis/analyze_STORM.m`, `Initialization/simulation_optimization_2.m` | **5** (the simulation format; on a RainbowSTORM file this reads `sigmax`) |
| `demo/Analysis/analyze_STORM.m`, `demo/Initialization/simulation_optimization_2.m`, `find_area.m` (both copies) | **12** |
| `param_opt2_5.m` (reads only simulated files) | 5 (correct) |

Check this before running either copy on your data. Channel file suffixes also differ: the original
instructions use `_C1` / `_C2`; the code uses `_me3` / `_ac`, except `demo/Analysis/analyze_STORM.m`, which
uses `_C1` / `_C2`.

## Missing inputs

### Removed, but recoverable from git

Commit `ce89e0c` ("added only code from hmt_v1", 2026-04-06) deleted the data files of a
complete Workflow 1 run. Its parent commit still has them:

- `Initialization/`: `initial_variables.mat`, `mean_ratio.mat`, `mean_mdens.mat`,
  `mean_adens.mat`, `ratio.mat`, `mdens2.mat`, `gt_values.mat`, both `*_epsilon_surround2.mat`,
  the `dens*`, `simdens*surround2*`, `Simulation*` and real-cell `well_*` CSVs
- `Analysis/`: per-cell results (`well_1_1.mat` … `well_2_2.mat`, `1141_com057.mat`), both ε
  files, `sample.mat`, `ring_binary1.tif`, and the matching CSVs
- `demo/Initialization/`: the same six input `.mat` files plus `gt_values.mat`

To restore one file without touching anything else:

```bash
git show ce89e0c~1:hmt_v1/Initialization/initial_variables.mat > initial_variables.mat
git ls-tree -r --name-only ce89e0c~1 hmt_v1    # list everything the parent commit had
```

| File | Needed by | Contents |
|---|---|---|
| `initial_variables.mat` | `simulation_optimization_2.m` | `ratiom_i`, `ratioa_i` (1 × 20) and `ratiomz_i`, `ratioaz_i` (1 × 100): starting per-ring correction factors. `meth_betweeni` (459 × 100), `ac_between_aci` (642 × 100): starting per-ring count tables |
| `mean_ratio.mat`, `mean_mdens.mat`, `mean_adens.mat` | `sim_histones_3d.m` | `mean_ratio`, `mdens`, `adens` (1 × 100): me3:ac ratio and densities from the nuclear edge inward, averaged over cells, which is the shape `compile_data.m` produces |
| `ratio.mat`, `mdens2.mat` | `sim_histones_3d.m` | `ratio`, `mdens`, `adens` (1 × 146): one cell's 50 nm-ring profiles, as written by `analyze_STORM.m`; used for the rings around nucleoli |

### Never in the repository

| File | Needed by | Holds (from how it is used) |
|---|---|---|
| `count_dist.mat` | `blinking_sim_unrealistic.m` | blink-count (`y647`, `y680`), blink-frame (`Nf647`, `Nf680`) and photon (`cumint647`, `cumint680`) distributions per dye. `Simulation/blink_stats/count_blinks.mat` is an empty file |
| `fwhm.mat` | `blinking_sim_unrealistic.m` | `zeroth_simulated`, `first_simulated`: PSF width vs. z (indexed at z + 4100 nm) |
| `new_cal.mat` | `make_first_order.m`, `make_real_fist.m` | `p1`, `p2`, `p3`: camera column → wavelength calibration (λ = p1·d² + p2·d + p3) |
| `Alexa_Fluor_647.csv`, `CF 680.csv` | `make_first_order.m` | dye emission spectra (wavelength in column 1, intensity in column 3) |
| `config_N_<int\|seg>_alph*.mat` | `master2.m` | `xpos`, `ypos`, `zpos`, `mask1`, written by `vis_SREV_int.m` / `vis_SREV_seg.m` (rename the output to the suffix `master2` expects) |
| `config-2.txt … config-20.txt` | `vis_SREV_*.m` | further SR-EV configurations; the SR-EV generator itself is not included either |
| `QD_analyze_A.m` | `master2.m` | matches zeroth- and first-order fits to each other and to ground truth, writing `<name>_fin.csv` |
| `saveastiff.m` | `master2.m` | File Exchange TIFF writer |

## Known issues

Found while documenting; nothing was changed. Line numbers refer to the files as committed.

**Affect results**

1. **z misaligned after masking.** `Initialization/find_within2.m` (lines 3–18) and
   `Analysis/analyze_STORM.m` (lines 221–236; same block in `demo/`) delete out-of-mask points
   from x and y but not from z. Every z after the first deleted point belongs to a different
   localization. This affects the axial statistics the simulation is fitted to (`mdistz`,
   `adistz`) and the z of every nanodomain centroid in `analyze_STORM.m`, which feeds its 3D
   centroid contact counts. `exclude_bg.m` shows the correct pattern.
2. **Density levels may all be the same** (`Initialization/sim_densities_3.m`).
   `sim_histones_3d.m` saves its whole workspace to `sim_param2.mat` (line 136), and
   `sim_histones_3d_2.m` loads it back (line 5). That overwrites the function's own arguments,
   including the ×n-scaled `ratcool*` factors that are supposed to make the 10 densities
   differ. `demo/Initialization/sim_densities_3.m` removes those variables from
   `sim_param2.mat` first (its lines 1–3); the main copy doesn't. Conversely, the main
   `sim_histones_3d_2.m` protects the output file name from the same overwrite and the demo copy
   doesn't, so in `demo/` all 10 levels are written to `Simulation.csv`, one over the next. Each
   copy therefore needs the other's fix. Check that the `simdens*` files exist and really differ
   in localization count.
3. **Dye identities shift** (`Simulation/sim_anti.m`). Sterically rejected fluorophores are
   removed from `f1x/f1y/f1z` (lines 68–70) but not from `identity2`, so after the first
   rejection every later fluorophore gets the AF647/CF680 identity of a different antibody.
4. **Antibody offset always points into +x +y +z** (`sim_anti.m` lines 35–41: all three components
   are `rand` ∈ [0, 1]).
5. **SR-EV coordinates scaled by 4.9 instead of 10** (`Simulation/master2.m` line 38 and its
   copies). SR-EV coordinates are in units of the minimum bond length, 10 nm; 4.9 nm is the
   bead's excluded-volume radius. Every simulated distance is about 2× too small (worked out
   in `hmt_v3/srev/io.py`, `parse_config_txt`).
6. **CF680 spectral window too red for real data** (`Simulation/splitmeac.m` line 7: 704.7–737.4 nm).
   The single-dye CF680 controls in `Simulation/blink_stats/` center at 699–701 nm, and only
   ~20 % of their localizations fall in this window. It suits the simulated movies, which put
   CF680 at 705 nm, but not real CF680 data.
7. **z column differs between copies** (see [Input file formats](#input-file-formats)).
8. **`Analysis/compile_data.m` loads `well_1_1.mat` four times**, so as committed it averages
   one cell with itself. Its filenames must be edited for each dataset (original step 16).
9. **Mixed 2D / 3D neighbour counts** in `analyze_STORM.m`'s centroid contact ratios:
   same-channel counts use 2D discs (lines 304, 339) while other-channel counts use 3D spheres
   (lines 318, 330). The same-channel counts also subtract 1 for "self" (lines 310, 346), although
   a centroid is not a localization.

**Make scripts fail or do nothing**

10. `demo/Analysis/analyze_STORM.m` saves `C1_within`, `C2_within`, … (line 509), which are
    never created, so the final `save` errors. `demo/Analysis/compile_data.m` expects the main
    copy's `meth_*` / `ac_*` names.
11. `demo/Initialization/` has no `split_me_ac.m`, which its `simulation_optimization_2.m` calls.
12. `master2.m` parameters that do nothing: `anti_length` (`sim_anti.m` hardcodes 9 nm and a
    25 % labelling efficiency), `SH` (passed to `make_real_fist` and ignored), `lambda` (the dye
    centroid is used instead), `SD`, `dens`, `numChains`, `tub_diam`, `numVerts`, `cutoff`. These
    are left over from a microtubule simulator.
13. `Simulation/add_noise.m` hardcodes `EM = 1`, so the EMCCD excess-noise factor is 1 and read
    noise is not divided by gain.
14. `master2.m` renders separate AF647-only and CF680-only movies (lines 82–83) with fresh random
    blinking, then never uses them. Its last cell uses `gt` after the loop has `clear`ed it.
15. `blinking_sim_unrealistic.m`: every blink lasts exactly one frame, and blinks drawn after
    `nframes` stay in the ground-truth table but are never rendered.
16. `param_opt2_5.m`: the masking block (lines 38–53) has no effect because the loop reloads
    unmasked data. Line 141 passes ac before me3 to `calc_decay`, so `decaym` / `decaya` are
    swapped; they are not saved. The variable `exp` shadows MATLAB's `exp`.
17. The original instructions' line references for the mask threshold (line 155, "0.17") and object size
    (line 158, "500000") match the main copy's line numbers but the demo copy's values. The main
    copy uses 0.05 and 100000; the demo copy has 0.17 and 500000 at lines 181 and 184.

## Duplicate and variant files

| Files | Relationship |
|---|---|
| `demo/Initialization/*` vs `Initialization/*` | 12 identical files; `denstest_1`, `simulation_optimization_2`, `sim_densities_3`, `sim_histones_3d_2` differ (see [demo/README.md](demo/README.md)); `split_me_ac`, `split_script_4` only in main |
| `demo/Analysis/*` vs `Analysis/*` | `compile_function.m` identical; `analyze_STORM.m` and `compile_data.m` differ |
| `Analysis/Copy_of_compile_data.m` | `compile_data.m` for two cells (`well_2_1`, `well_2_2` → `well_2.mat`) |
| `Initialization/split_me_ac.m` vs `split_script_4.m` | same split; function for `Simulation.csv` vs script for the 10 `simdens` files |
| `Initialization/find_within.m` vs `find_within2.m` | `find_within2` first restricts to the nucleus mask |
| `Initialization/sim_histones_3d.m` vs `sim_histones_3d_2.m` | `_2` skips building the nucleus and reuses the weighting maps saved by the first |
| `Simulation/master2.m`, `Copy_of_master2.m`, `Copy_2_of_master2.m`, `Copy_3_of_master2.m`, `master_sim_tub_1frame.m` | the same script pointed at different segmentations (`int_alph1p0`, `seg_alph1p1`, `int_alph1p2`, `seg_alph1p2`, `seg_alph1p1`) with different batch start indices; see [Simulation/README.md](Simulation/README.md) |
| `Simulation/gen_ovito_fin.m` vs `gen_ovito_fin2.m` | `_2` also writes fluorophore positions; not called anywhere |

## Where each piece lives in hmt_v3

| hmt_v1 | hmt_v3 |
|---|---|
| nucleus mask in `analyze_STORM.m` / `simulation_optimization_2.m` | `preprocess.binarize_nucleus`, `preprocess.auto_mask_threshold` |
| 50 nm edge rings, me3:ac ratio vs. radius | `preprocess.create_radial_contours`, `measure.extract_radial_density_profile` |
| `find_within2`, `create_variables`, `sim_histones_3d`, `sim_nucleus_3d` | `simulate.generate_nucleus` (fitted to the real cell's mask, radial profile, RDF/ADF and pair correlation; no synthetic nucleus outline) |
| `find_gt_width` | `measure.calculate_geometric_gt_target` (2 × 75th-percentile radius) |
| `epsilon_cost`, `param_opt2_5` | `calibrate.domain_size_cost`, `calibrate.optimize_dbscan_size_matching` (also fits min_samples) |
| density-corrected DBSCAN in `analyze_STORM` | `cluster.cluster_dbscan_density_corrected` |
| contact ratios, `calc_decay`, `compile_function` exp fit | `postprocess.nanodomain_sncr`, `postprocess.fit_exponential_decay` |
| `vis_SREV_int`, `vis_SREV_seg` | `srev.compartments.classify_compartments` (coordination number) |
| `sim_anti` | `srev.labeling` |
| `blinking_sim_unrealistic`, `add_noise` | `srev.photophysics`, `srev.render` |
| `make_first_order`, `make_real_fist` | `srev.spectral` |
| `gen_ovito_fin` | `srev.visualize.export_for_ovito` |
| `Simulation/blink_stats/` | not read at run time; the dye centroids, PSF widths and zeroth:first photon split measured from these files are `srev.spectral`'s defaults |
