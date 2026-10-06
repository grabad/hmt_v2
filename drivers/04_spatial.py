"""
Step 04 -- per-cell spatial analysis (Methods: "Histone spatial association analysis").

For each cell:
  1. radial profile (Fig. 3c, 6c/f): localization density of each mark in the
     100 contour bands following the nuclear outline, 0 = core ... 99 = periphery
  2. SNCR (Fig. 3d-f, 6b/e): around every nanodomain of one mark, the ratio of
     same-mark to other-mark neighbours within a growing radius r, for both
     neighbour definitions in SNCR_NEIGHBORS (nanodomain centroids, or
     localizations), each summarized as the per-domain mean and the pooled ratio
  3. exponential decay rate of the SNCR curve (Fig. S3b/c), fit to the
     SNCR_SUMMARY curve for every (reference mark, neighbour definition)

Writes cells/<cell_id>/radial_profile.csv, sncr.csv, sncr_fit.csv, spatial.png.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C
from common import cell_dir, hmt, load_cell, load_manifest, log, parse_args, xyz_cols

CENTROID_COLS = {2: ["centroid_x", "centroid_y"], 3: ["centroid_x", "centroid_y", "centroid_z"]}


def radial_profile(cell, cell_id, condition):
    rows = []
    for channel in C.CHANNELS:
        profile = hmt.simulate.extract_radial_density_profile(
            cell[channel], cell["contour_bands"], cell["x_min"], cell["y_min"], bin_size=cell["bin_nm"])
        peak = profile.max() if profile.max() > 0 else 1.0
        for band, value in enumerate(profile):
            rows.append({"cell_id": cell_id, "condition": condition, "channel": channel, "band": band,
                         "radius_pct": band + 0.5,  # band b spans b..b+1 % of the core->rim distance
                         "density_per_um2": value / cell["bin_nm"] ** 2 * 1e6,
                         "density_norm": value / peak})
    return pd.DataFrame(rows)


def sncr(cell_id, condition, out):
    dims = 3 if C.SNCR_USE_Z else 2
    nanodomains = {ch: pd.read_csv(out / f"{ch}_nanodomains.csv") for ch in C.CHANNELS}
    locs = {ch: pd.read_parquet(out / f"{ch}_clustered.parquet") for ch in C.CHANNELS}
    me3, ac = C.CHANNELS

    curves, fits = [], []
    for ref_ch, other_ch in ((me3, ac), (ac, me3)):
        ref = nanodomains[ref_ch][CENTROID_COLS[dims]].to_numpy()
        for mode in C.SNCR_NEIGHBORS:
            if mode == "nanodomains":
                self_pts = ref
                nonself_pts = nanodomains[other_ch][CENTROID_COLS[dims]].to_numpy()
                exclude_ref = True
            elif mode == "localizations":
                self_pts = locs[ref_ch][xyz_cols(C.SNCR_USE_Z)].to_numpy()
                nonself_pts = locs[other_ch][xyz_cols(C.SNCR_USE_Z)].to_numpy()
                exclude_ref = False
            else:
                raise ValueError(f"unknown SNCR_NEIGHBORS entry {mode!r}")

            res = hmt.postprocess.nanodomain_sncr(ref, self_pts, nonself_pts, C.SNCR_RADII_NM,
                                                  exclude_ref=exclude_ref)
            curves.append(pd.DataFrame({
                "cell_id": cell_id, "condition": condition, "ref_channel": ref_ch, "neighbors": mode,
                "r_nm": res["r"], "mean_ratio": res["mean_ratio"], "pooled_ratio": res["pooled_ratio"],
                "chance_ratio": res["chance_ratio"], "mean_self": res["mean_self"],
                "mean_nonself": res["mean_nonself"], "frac_no_nonself": res["frac_no_nonself"],
                "n_ref": res["n_ref"]}))
            fit = hmt.postprocess.fit_exponential_decay(res["r"], res[C.SNCR_SUMMARY], r_max=C.SNCR_FIT_MAX_NM)
            fits.append({"cell_id": cell_id, "condition": condition, "ref_channel": ref_ch, "neighbors": mode,
                         "summary": C.SNCR_SUMMARY, "fit_r_max_nm": C.SNCR_FIT_MAX_NM, **fit})
    return pd.concat(curves, ignore_index=True), pd.DataFrame(fits)


def plot_cell(radial, curves, cell_id, condition, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    prim = curves[curves["neighbors"] == C.SNCR_PRIMARY]
    for ch in C.CHANNELS:
        g = prim[prim["ref_channel"] == ch]
        y = g[C.SNCR_SUMMARY].to_numpy()
        if np.any(np.isfinite(y)) and np.nanmax(y) > 0:
            axes[0].plot(g["r_nm"], y / np.nanmax(y), color=C.CHANNEL_COLORS[ch], label=C.CHANNEL_LABELS[ch])
        r = radial[radial["channel"] == ch]
        axes[1].plot(r["radius_pct"], r["density_norm"], color=C.CHANNEL_COLORS[ch], label=C.CHANNEL_LABELS[ch])
    axes[0].set(xlabel="Distance (nm)", ylabel="SNCR (normalized to max)", ylim=(0, 1.05),
                title=f"SNCR ({C.SNCR_PRIMARY}, {C.SNCR_SUMMARY})")
    axes[1].set(xlabel="Radius (%)  core → periphery", ylabel="Normalized localization density",
                title="Radial contour density")
    for ax in axes:
        ax.legend()
    fig.suptitle(f"{cell_id} [{condition}]")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    args = parse_args(__doc__.splitlines()[1])
    manifest = load_manifest(args)
    for _, row in manifest.iterrows():
        cell_id, condition = row["cell_id"], row["condition"]
        out = cell_dir(cell_id)
        if not (out / "cluster.json").exists():
            raise FileNotFoundError(f"{out / 'cluster.json'} missing -- run 03_cluster.py")
        if (out / "sncr_fit.csv").exists() and not args.overwrite:
            log(f"{cell_id}: spatial analysis exists, skipping (--overwrite to redo)")
            continue

        radial = radial_profile(load_cell(cell_id), cell_id, condition)
        curves, fits = sncr(cell_id, condition, out)
        radial.to_csv(out / "radial_profile.csv", index=False)
        curves.to_csv(out / "sncr.csv", index=False)
        fits.to_csv(out / "sncr_fit.csv", index=False)
        plot_cell(radial, curves, cell_id, condition, out / "spatial.png")

        prim = fits[fits["neighbors"] == C.SNCR_PRIMARY].set_index("ref_channel")["k"]
        log(f"{cell_id} [{condition}]: SNCR decay rate ({C.SNCR_PRIMARY}) "
            + ", ".join(f"{C.CHANNEL_LABELS[ch]} {prim.get(ch, np.nan):.4f}/nm" for ch in C.CHANNELS))


if __name__ == "__main__":
    main()
