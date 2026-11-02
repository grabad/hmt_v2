# Simulation

Builds synthetic two-colour spectroscopic SMLM (sSMLM) camera movies from an SR-EV chromatin
polymer configuration, so the analysis can be tested against known ground truth.
Heterochromatin is labelled with Alexa Fluor 647 (the "me3" channel) and euchromatin with CF680
(the "ac" channel). Each movie has a zeroth-order (spatial) stack and a first-order (spectrally
dispersed) stack, with noise.

Several inputs are not in the repository (`count_dist.mat`, `fwhm.mat`, `new_cal.mat`, the dye
spectra, `QD_analyze_A.m`, `saveastiff.m`), so this folder does not run end to end as committed.
See [../DOCUMENTATION.md](../DOCUMENTATION.md) for those and for the known issues.

## Pipeline

```
config-N.txt ─► vis_SREV_int.m | vis_SREV_seg.m ─► config_N_<int|seg>_alph*.mat (+ .xyz for OVITO)
             ─► master2.m ─► <name>_z.tif, _f.tif (noise-free), <name>_zerothh.tif, _firsth.tif (noisy)
             ─► fit the movies externally (ThunderSTORM) ─► QD_analyze_A (missing) ─► <name>_fin.csv
             ─► splitmeac.m ─► <name>_fin_me3.csv, <name>_fin_ac.csv
```

## Data in this folder

| File | Contents |
|---|---|
| `config-1.txt` | One SR-EV configuration: an atom count (186 742), the comment `Atoms. Timestep: 1300`, then `type x y z` per monomer. Coordinates are in SR-EV units of 10 nm (see `hmt_v3/srev/io.py`). Only config-1 is included; the scripts list config-1…20. |
| `blink_stats/` | Single-dye control acquisitions, three runs each of AF647 and CF680. `<dye>_<n>.csv`: RainbowSTORM matched localizations (15 columns, including `centroid [nm]` and `z [nm]`). `<dye>_<n>_zeroth.csv`, `_first.csv`: ThunderSTORM fits of the zeroth- and first-order images. Run 2 has only the zeroth/first files. Median centroids are 670 nm (AF647) and 699–701 nm (CF680). `count_blinks.mat` is empty. No MATLAB script here reads these files; `hmt_v3.srev.spectral` takes its defaults from them. |
| `txt2mat.m` | Third-party ASCII table reader (MATLAB File Exchange, © 2005–2021 Andres), used to read `config-N.txt` |

## Compartment labelling

### `vis_SREV_int.m` (script)
For each configuration (read with `txt2mat`; columns type, x, y, z), it splits the monomers into
compartments by k-means:

- k-means with k = 10 on monomer positions; each cluster's density = count / bounding-box volume.
- The 5 sparsest clusters → euchromatin (`mask1 = 1`); the rest → heterochromatin.
- `integration` (0–1, line 6) then flips the label of a random 0.5 × `integration` fraction of
  monomers, to mix the compartments.

Writes `config_N_seg_alph1p2.xyz` (via `gen_ovito_fin`) and `config_N_seg_alph1p2.mat` (`eux`,
`euy`, `euz`, `hetx`, `hety`, `hetz`, `xpos`, `ypos`, `zpos`, `mask1`). Despite the `seg` in the
output names, this is the k-means ("int") variant; `master2.m` expects `config_N_int_alph1p0.mat`,
so the suffix was evidently edited per run. `config-13` is listed twice; `ndomain` is unused.

### `vis_SREV_seg.m` (script)
Local-packing variant. For each monomer it takes the mean distance to its 10 nearest neighbours,
searching a box of 1/70 of the configuration's extent and doubling the box until it holds
≥ 10 monomers. Monomers with mean distance > 1.165 (SR-EV units) → euchromatin. Same
`integration` option and outputs (`config_N_seg_alph1p1.*`, plus `mdistn`). As committed the
loop runs only the first file (`for file=1:1`).

Both scripts call the sparser chromatin euchromatin.

### `gen_ovito_fin(eux, euy, euz, hetx, hety, hetz, outputFileName)`
Writes an XYZ file for OVITO: the particle count, a comment line, then `<type> x y z 0.8` per
monomer. Heterochromatin is type `He` and euchromatin `Li`; element names serve as type labels,
and the 5th column is transparency.

### `gen_ovito_fin2(…, f1x, f1y, f1z, outputFileName)`
Same, plus fluorophore positions as type `Be` (transparency 0). Not called anywhere.

## Movie generation

### `master2.m` (script)
For each configuration in `simnames` / `simnames2` (the loop starts at entry 17 in this copy,
evidently resuming a batch):

1. Loads `config_N_int_alph1p0.mat`, multiplies the coordinates by 4.9 to get nm (see Known
   issue 5), writes `coordinates.txt` and centres z.
2. `identity = 2 − mask1`: heterochromatin → 2 (AF647), euchromatin → 1 (CF680).
3. `sim_anti` → fluorophore positions.
4. `blinking_sim_unrealistic` → zeroth-order movie `data` and ground-truth table `gt`. It also
   renders AF647-only and CF680-only movies, which are independent random draws and never used.
5. `make_first_order` → dispersed dye spectra. `make_real_fist` once per dye → first-order
   movie, summed and cropped to the zeroth-order size.
6. Saves `<name>_z.tif` and `<name>_f.tif`, runs `add_noise`, then saves `<name>_zerothh.tif`
   and `<name>_firsth.tif` (with `saveastiff`, not included).

Parameters that are used: `px` = 100 nm, `NA` = 1.4, `FOV` = 5 µm, `nframes` = 25 000, `Ibgp` = 10
and `Ibg` = 15 (background photons/px, zeroth and first order), `RN` = 1 (read noise). Unused,
left over from a microtubule simulator: `lambda`, `SH`, `SD`, `dens`, `numChains`, `anti_length`,
`tub_diam`, `numVerts`, `cutoff`.

The cells after the loop are run separately, once the TIFFs have been fitted outside MATLAB:
`QD_analyze_A` (missing) for each movie → `<name>_fin.csv`, then `splitmeac`. The last cell is
scratch analysis of frames with exactly two emitters; it needs `gt`, which the loop has cleared.

| Variant | Input segmentation | Differences |
|---|---|---|
| `master2.m` | `config_N_int_alph1p0` | loop from 17; `QD_analyze_A` for movies 1–4 |
| `Copy_of_master2.m` | `config_N_seg_alph1p1` | loop over 2–7; noise-free TIFFs not saved |
| `Copy_2_of_master2.m` | `config_N_int_alph1p2` | loop from 4 |
| `Copy_3_of_master2.m` | `config_N_seg_alph1p2` | loop from 16 |
| `master_sim_tub_1frame.m` | `config_N_seg_alph1p1` | loop from 1; `QD_analyze_A` from movie 4. Despite the name, this is not a tubulin or single-frame simulation |

### `sim_anti(tempx, tempy, zpos, identity) → [gtx, gty, gtz, f1x, f1y, f1z, identity2]`
Antibody labelling.

- Each monomer is labelled with probability 0.25.
- The fluorophore position is built from:
  - the monomer plus a 5.83 nm anchor offset (√(3² + 5²))
  - a 9 nm primary antibody (random direction, elevation within ±135°)
  - a 9 nm secondary antibody (random direction)
- A fluorophore is discarded if its antibody path (10 sample points per segment) passes within
  4.9 nm of any other monomer.
- **Returns:** every monomer (`gt*`), the kept fluorophores (`f1*`), and the dye identity of each
  labelled monomer (`identity2`).
- **Issues:**
  - Rejected fluorophores are not removed from `identity2`, so later identities are shifted.
  - The anchor offset always points into +x +y +z.
  - Labelling efficiency and antibody length are hardcoded; `master2.m`'s `anti_length` does
    nothing.

### `blinking_sim_unrealistic(xpos, ypos, zpos, px, NA, nframes, size, identity2, fnameout) → [blinking_data, gt]`
Zeroth-order movie.

- **Per fluorophore:**
  - The number of blinks, the frame of each blink (in 100-frame bins) and photons per blink are
    drawn from distributions in `count_dist.mat`.
  - Spectral centroid: 671 nm (AF647, `identity2 = 2`) or 705 nm (CF680).
  - PSF widths `sig1` / `sig2` (zeroth / first order) come from `fwhm.mat`, at index z + 4100 nm,
    divided by 2.75.
- **Rendering:**
  - Each blink lasts one frame.
  - Each frame sums diffraction-limited PSFs: the FFT of a circular pupil with cutoff
    NA / (λ · σ(z) / σ_min), where λ is the centroid. Each PSF is shifted to the emitter with
    `imtranslate` and scaled by its photon count.
  - The image is `size` · 1000 / `px` pixels square (50 × 50 for 5 µm at 100 nm), centred on the
    mean emitter position.
  - Blinks after `nframes` are kept in `gt` but not rendered.
- **`gt` columns:** id, frame, x, y, z, centroid, sig1, sig2, photons.
- **Notes:** the argument `size` shadows MATLAB's `size()`. `fnameout` is unused because the
  `saveastiff` call is commented out.

### `make_first_order` (script, runs inside `master2.m`'s workspace)
Reads the dye spectra `Alexa_Fluor_647.csv` and `CF 680.csv` (wavelength in column 1, intensity
in column 3) and the calibration `new_cal.mat` (`p1`, `p2`, `p3`). It maps camera columns
900–1100 to wavelength with λ = p1·d² + p2·d + p3, then resamples each spectrum onto those
columns and normalises it to sum 1, giving `sptimg647` and `sptimg680`.

### `make_real_fist(first_ord, spec_het, gt, FOV, px, NA, nframes) → blinking_first`
First-order (spectral) movie. For each blink in `gt`:

- builds a PSF with the first-order width (`sig2`)
- shifts it like the zeroth order in y, and in x by (d − 950) px, where d is the camera column
  of the blink's centroid wavelength (the inverse of the calibration)
- convolves it with the dye's dispersed spectrum, scales it by photons and adds it to the
  blink's frame

`spec_het` is unused. The documented second output `gtfist` is never assigned. Lines 14–16 are
leftover debugging.

### `add_noise` (script, runs inside `master2.m`'s workspace)
Camera noise on both movies: Poisson shot noise on the signal and on a flat background (`Ibgp`
photons/px zeroth, `Ibg` first), plus Gaussian read noise σ = `RN`. These are combined as
√(2 − 1/EM) · shot noise + read noise / EM + signal. `EM = 1` is hardcoded, so there is no
EM excess-noise factor. Produces `data_noise` and `first_noise`; `split` and `eta` are unused.

### `splitmeac(fnam)`
Reads `<fnam>_fin.csv` (17-column format) and splits it by spectral centroid (column 10):

- 655.62–688.35 nm → `<fnam>_fin_me3.csv` (AF647)
- 704.715–737.445 nm → `<fnam>_fin_ac.csv` (CF680)
- everything else is dropped

The CF680 window suits the simulation (705 nm) but holds only ~20 % of the real CF680 control
localizations in `blink_stats/`, which center at ~700 nm.
