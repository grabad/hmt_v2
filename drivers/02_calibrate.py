"""
Step 02 -- inverse Monte Carlo DBSCAN calibration (manuscript Fig. 4).

For each calibration source (CALIBRATION_MODE: the densest cell of each
condition, every cell, or the densest cell overall) and each mark:
  1. generate_nucleus: rebuild a simulated nucleus that matches the real one's
     radial distribution, within-domain point distribution (RDF/ADF) and
     inter-domain arrangement (g(r)/Ripley's L) -- with known per-localization
     ground-truth domain labels
  2. crop a CROP_SIZE_NM square inside the nucleus whose density matches the
     simulated nucleus's mean density (so crop density and a real cell's
     locs/mask-area are measured the same way)
  3. subsample the crop to DENSITY_FRACTIONS (fewer localizations, same domains)
  4. optimize_dbscan_size_matching: per fraction, the (eps, min_samples) whose
     DBSCAN clusters best match the ground-truth domain size; fit
     eps = a * density^b and min_samples = m * density + c
  5. Fig. 4b/c check: cluster every fraction with FIXED params (the corrected
     params at full density) vs. DENSITY-CORRECTED params, and report how
     strongly measured size tracks density (R^2) for each

Writes calibration/<key>/<ch>_calibration.json (what step 03 reads), plus the
simulated nucleus, the size-vs-density table and diagnostic figures.

This is the slow step: DBSCAN runs ~(#fractions x #min_samples x ~30) times per
mark. Use --quick first to check the chain end to end.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

import config as C
from common import (calibration_dir, calibration_key, capture_figures, hmt, load_cell, load_json,
                    load_manifest, log, parse_args, save_json)


def calibration_sources(manifest):
    """{calibration key: cell_id to calibrate on}, picking the densest cell where needed."""
    dens = []
    for _, row in manifest.iterrows():
        pre = load_json(f"{C.OUTPUT_DIR}/cells/{row['cell_id']}/preprocess.json")
        dens.append(pre["me3_density_per_um2"] + pre["ac_density_per_um2"])
    manifest = manifest.assign(density=dens)

    if C.CALIBRATION_MODE == "per_cell":
        return dict(zip(manifest["cell_id"], manifest["cell_id"]))
    if C.CALIBRATION_MODE == "pooled":
        return {"pooled": manifest.loc[manifest["density"].idxmax(), "cell_id"]}
    return {cond: grp.loc[grp["density"].idxmax(), "cell_id"]
            for cond, grp in manifest.groupby("condition")}


def choose_crop(sim_df, binary_mask, x_min, y_min, bin_nm, size_nm, min_inside):
    """
    Square window, at least `min_inside` inside the mask, whose localization
    density (per in-mask area) is closest to the whole simulated nucleus's mean.
    Shrinks the window if the nucleus is too small to hold it.

    Returns:
        (bounds, area_nm2): bounds = (xmin, xmax, ymin, ymax) in nm; area_nm2 = the
        in-mask area inside the window (the density denominator)
    """
    mask = binary_mask.astype(float)
    nx, ny = mask.shape
    x_edges = x_min + np.arange(nx + 1) * bin_nm
    y_edges = y_min + np.arange(ny + 1) * bin_nm
    counts, _, _ = np.histogram2d(sim_df["x [nm]"], sim_df["y [nm]"], bins=(x_edges, y_edges))
    target = len(sim_df) / (mask.sum() * bin_nm ** 2)

    def box_sums(img, w):
        s = np.pad(img, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
        return s[w:, w:] - s[:-w, w:] - s[w:, :-w] + s[:-w, :-w]  # [i, j] = window starting at pixel (i, j)

    size = size_nm
    while size >= 1500:
        w = int(round(size / bin_nm))
        if w < nx and w < ny:
            inside = box_sums(mask, w)
            n = box_sums(counts, w)
            ok = inside >= min_inside * w * w
            if ok.any():
                density = np.where(ok, n / np.maximum(inside * bin_nm ** 2, 1), np.nan)
                i, j = np.unravel_index(np.nanargmin(np.abs(density - target)), density.shape)
                bounds = (x_edges[i], x_edges[i + w], y_edges[j], y_edges[j + w])
                return bounds, float(inside[i, j] * bin_nm ** 2)
        size *= 0.8
    raise ValueError("nucleus too small for any calibration crop >= 1.5 um")


def crop(df, bounds):
    xmin, xmax, ymin, ymax = bounds
    return df[df["x [nm]"].between(xmin, xmax, inclusive="left") &
              df["y [nm]"].between(ymin, ymax, inclusive="left")].copy()


def mean_sizes(clustered_df):
    """Mean nanodomain size (manuscript bounding-box size and calibration diameter)."""
    nd = hmt.cluster.calculate_nanodomain_characteristics(clustered_df, percentile=C.SIZE_PERCENTILE)
    nd = nd[nd["n_locs"] >= 10]
    if nd.empty:
        return np.nan, np.nan, 0
    return nd["bb_size_2D"].mean(), nd["geometric_diameter_2D"].mean(), len(nd)


def size_vs_density_check(fraction_dfs, area_nm2, calibration):
    """Fig. 4b/c: fixed vs density-corrected DBSCAN across the density fractions."""
    _, fixed = hmt.cluster.cluster_dbscan_density_corrected(
        fraction_dfs[0], calibration, area_nm2, use_z=C.USE_Z,
        min_eps=C.MIN_EPS_NM, min_min_samples=C.MIN_MIN_SAMPLES)
    rows = []
    for frac, df in zip(C.DENSITY_FRACTIONS, fraction_dfs):
        fixed_df = hmt.cluster.cluster_dbscan(df, eps=fixed["eps"], min_samples=fixed["min_samples"],
                                              use_z=C.USE_Z)
        corr_df, corr = hmt.cluster.cluster_dbscan_density_corrected(
            df, calibration, area_nm2, use_z=C.USE_Z, min_eps=C.MIN_EPS_NM, min_min_samples=C.MIN_MIN_SAMPLES)
        for method, clustered, params in (("fixed", fixed_df, fixed), ("corrected", corr_df, corr)):
            bb, diam, n = mean_sizes(clustered)
            rows.append({"fraction": frac, "density_per_um2": len(df) / area_nm2 * 1e6, "method": method,
                         "eps": params["eps"], "min_samples": params["min_samples"],
                         "mean_bb_size_nm": bb, "mean_diameter_nm": diam, "n_domains": n})
    table = pd.DataFrame(rows)

    r2 = {}
    for method, grp in table.groupby("method"):
        g = grp.dropna(subset=["mean_bb_size_nm"])
        r2[method] = linregress(g["density_per_um2"], g["mean_bb_size_nm"]).rvalue ** 2 if len(g) > 2 else np.nan
    return table, r2


def plot_size_vs_density(table, r2, gt_diameter, channel, path):
    color = C.CHANNEL_COLORS[channel]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, method, title in ((axes[0], "fixed", "Fixed eps/min_samples (Fig. 4b)"),
                              (axes[1], "corrected", "Density-corrected (Fig. 4c)")):
        g = table[table["method"] == method]
        ax.plot(g["density_per_um2"], g["mean_bb_size_nm"], "o", color=color, label="bounding box")
        ax.plot(g["density_per_um2"], g["mean_diameter_nm"], "s", mfc="none", color=color, label="p75 diameter")
        ax.axhline(gt_diameter, color="0.4", ls=":", label=f"GT p75 diameter ({gt_diameter:.0f} nm)")
        ax.set_title(f"{title}\nR$^2$ (bb size vs density) = {r2.get(method, np.nan):.3f}")
        ax.set_xlabel("Localization density (locs/um^2)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(f"{C.CHANNEL_LABELS[channel]} nanodomain size (nm)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def calibrate_channel(cell, channel, out, rng_seed):
    log(f"  [{channel}] generate_nucleus (domain_scale={C.DOMAIN_SCALE_NM[channel]})")
    _, sim_df, info = hmt.simulate.generate_nucleus(
        cell[channel], cell["contour_bands"], cell["x_min"], cell["y_min"], px_size=cell["bin_nm"],
        domain_scale=C.DOMAIN_SCALE_NM[channel], noise_fraction=C.NOISE_FRACTION[channel],
        seed_placement=C.SEED_PLACEMENT, rng=np.random.default_rng(rng_seed), verbose=True)
    sim_df.to_parquet(out / f"{channel}_simulated.parquet", index=False)
    with capture_figures(out / f"{channel}_pair_correlation"):
        hmt.visualize.plot_pair_correlation(info["pcf_real"], info["pcf_noisy"], r_domain=info["r_domain_nm"],
                                            color="green" if channel == "me3" else "red")

    bounds, area_nm2 = choose_crop(sim_df, cell["binary_mask"], cell["x_min"], cell["y_min"], cell["bin_nm"],
                                   C.CROP_SIZE_NM, C.CROP_MIN_INSIDE)
    crop_df = crop(sim_df, bounds)
    log(f"  [{channel}] crop {bounds[1] - bounds[0]:.0f} nm square: {len(crop_df):,} locs, "
        f"{len(crop_df) / area_nm2 * 1e6:.0f} locs/um^2 (nucleus mean {len(sim_df) / cell['area_nm2'] * 1e6:.0f})")

    packing = hmt.simulate.diagnose_domain_packing(crop_df, domain_scale=C.DOMAIN_SCALE_NM[channel], use_z=C.USE_Z)
    log(f"  [{channel}] {packing['n_domains']} GT domains in crop, "
        f"{packing['frac_overlapping']:.0%} overlap their nearest neighbour")

    rng = np.random.default_rng(rng_seed)
    fraction_dfs = [crop_df if f >= 1.0 else crop_df.sample(frac=f, random_state=int(rng.integers(2**31)))
                    for f in C.DENSITY_FRACTIONS]

    log(f"  [{channel}] optimize_dbscan_size_matching: {len(fraction_dfs)} fractions x "
        f"{len(C.MIN_SAMPLES_RANGE)} min_samples")
    with capture_figures(out / f"{channel}_calibration_fit"):
        calibration = hmt.simulate.optimize_dbscan_size_matching(
            fraction_dfs, area_nm2, min_samples_range=C.MIN_SAMPLES_RANGE, percentile=C.SIZE_PERCENTILE,
            cost_threshold=C.SIZE_COST_THRESHOLD_NM, use_z=C.USE_Z, show_plot=True, verbose=False)
    n_fit = len(calibration["densities"])
    log(f"  [{channel}] eps = {calibration['eps_coef']:.3g} * density^{calibration['eps_exp']:.3f} "
        f"(R^2 {calibration['eps_r_squared']:.3f}); min_samples = {calibration['min_samples_slope']:.4g} * density + "
        f"{calibration['min_samples_intercept']:.2f}; {n_fit}/{len(fraction_dfs)} fractions fit, "
        f"GT diameter {calibration['gt_target']:.0f} nm")

    calibration.update({
        "calibration_area_nm2": area_nm2,
        "crop_bounds_nm": bounds,
        "use_z": C.USE_Z,
        "density_fractions": C.DENSITY_FRACTIONS,
        "min_samples_range": list(C.MIN_SAMPLES_RANGE),
        "generate_nucleus": {k: v for k, v in info.items() if np.isscalar(v) or v is None},
        "domain_scale_nm": C.DOMAIN_SCALE_NM[channel],
        "gt_frac_domains_overlapping": packing["frac_overlapping"],
        "converged": n_fit >= 2,
    })

    if n_fit >= 2:
        table, r2 = size_vs_density_check(fraction_dfs, area_nm2, calibration)
        table.to_csv(out / f"{channel}_size_vs_density.csv", index=False)
        plot_size_vs_density(table, r2, calibration["gt_target"], channel, out / f"{channel}_size_vs_density.png")
        calibration["size_vs_density_r2"] = r2
        failed = table[table["mean_bb_size_nm"].isna()].groupby("method").size().to_dict()
        note = (f" (no nanodomain found at {failed} fraction(s) -- NaN rows, left out of R^2)"
                if failed else "")
        log(f"  [{channel}] size-vs-density R^2: fixed {r2.get('fixed', np.nan):.3f} -> "
            f"corrected {r2.get('corrected', np.nan):.3f}{note}")
    else:
        log(f"  [{channel}] WARNING: calibration did not converge (<2 fractions under "
            f"{C.SIZE_COST_THRESHOLD_NM:.0f} nm size error) -- see {channel}_calibration_fit.png; "
            f"try a different DOMAIN_SCALE_NM or a larger SIZE_COST_THRESHOLD_NM")
    return calibration


def main():
    args = parse_args(__doc__.splitlines()[1])
    manifest = load_manifest()  # always pick sources from ALL cells
    sources = calibration_sources(manifest)
    if args.cells:  # only the calibrations those cells need
        wanted = {calibration_key(r["cell_id"], r["condition"])
                  for _, r in manifest[manifest["cell_id"].isin(args.cells)].iterrows()}
        sources = {k: v for k, v in sources.items() if k in wanted}
    log(f"CALIBRATION_MODE={C.CALIBRATION_MODE}: {sources}")

    for i, (key, cell_id) in enumerate(sources.items()):
        out = calibration_dir(key)
        if all((out / f"{ch}_calibration.json").exists() for ch in C.CHANNELS) and not args.overwrite:
            log(f"calibration '{key}' exists, skipping (--overwrite to redo)")
            continue
        log(f"calibration '{key}' from {cell_id}")
        cell = load_cell(cell_id)
        for j, channel in enumerate(C.CHANNELS):
            calibration = calibrate_channel(cell, channel, out, rng_seed=C.RANDOM_SEED + 100 * i + j)
            calibration["source_cell"] = cell_id
            save_json(calibration, out / f"{channel}_calibration.json")


if __name__ == "__main__":
    main()
