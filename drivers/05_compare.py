"""
Step 05 -- pool cells and compare conditions (manuscript Fig. 6 and Fig. S3).

Collects every cell's outputs from steps 01-04 into results/:
    per_cell_summary.csv    one row per cell: localization counts, ME:AC ratio, per-mark
                            nanodomain count / size / DBSCAN params, SNCR decay rates
    nanodomains_all.csv     every nanodomain of every cell
    sncr_all.csv, radial_profiles_all.csv
    stats.csv               each condition vs CONDITION_ORDER[0]: n, mean, sd, median,
                            Mann-Whitney U and Welch t p-values, per metric
and figures:
    fig6g_nanodomain_size.png, fig6h_nanodomain_count.png   per-cell means by condition
    fig6_sncr_by_condition.png                              mean +/- SEM SNCR curves (6b/6e)
    fig6_radial_by_condition.png                            mean +/- SEM radial profiles (6c/6f)
    figS3_ratio_and_decay.png                               ME:AC ratio, SNCR decay rates (S3a-c)
Each dot in the per-cell plots is one nucleus.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import config as C
from common import cells_root, load_json, load_manifest, log, parse_args, results_dir


def collect(manifest):
    per_cell, nanodomains, sncr, fits, radial = [], [], [], [], []
    for _, row in manifest.iterrows():
        cell_id, d = row["cell_id"], cells_root() / row["cell_id"]
        needed = [d / "preprocess.json", d / "cluster.json", d / "sncr_fit.csv"]
        if not all(p.exists() for p in needed):
            log(f"WARNING: {cell_id} is missing step 01-04 outputs; left out")
            continue

        pre, clu = load_json(d / "preprocess.json"), load_json(d / "cluster.json")
        cell_fits = pd.read_csv(d / "sncr_fit.csv")
        rec = {k: pre[k] for k in ("cell_id", "condition", "mask_area_um2", "n_me3", "n_ac", "me3_ac_ratio")}
        rec["calibration"] = clu["calibration"]
        for ch in C.CHANNELS:
            for k in ("n_nanodomains", "mean_bb_size_nm", "median_bb_size_nm", "mean_diameter_nm",
                      "frac_locs_clustered", "density_per_um2", "eps", "min_samples", "extrapolated"):
                rec[f"{ch}_{k}"] = clu[ch][k]
            k = cell_fits[(cell_fits["ref_channel"] == ch) & (cell_fits["neighbors"] == C.SNCR_PRIMARY)]["k"]
            rec[f"{ch}_sncr_decay_rate"] = float(k.iloc[0]) if len(k) else np.nan
        per_cell.append(rec)

        for ch in C.CHANNELS:
            nd = pd.read_csv(d / f"{ch}_nanodomains.csv")
            nanodomains.append(nd.assign(cell_id=cell_id, condition=row["condition"], channel=ch))
        sncr.append(pd.read_csv(d / "sncr.csv"))
        fits.append(cell_fits)
        radial.append(pd.read_csv(d / "radial_profile.csv"))

    if not per_cell:
        raise RuntimeError("no cell has complete outputs -- run steps 01-04 first")
    cat = lambda frames: pd.concat(frames, ignore_index=True)
    return pd.DataFrame(per_cell), cat(nanodomains), cat(sncr), cat(fits), cat(radial)


def compare_conditions(per_cell, metrics):
    ref = C.CONDITION_ORDER[0]
    rows = []
    for metric in metrics:
        a = per_cell.loc[per_cell["condition"] == ref, metric].dropna().astype(float)
        for cond in [c for c in C.CONDITION_ORDER[1:] if c in set(per_cell["condition"])]:
            b = per_cell.loc[per_cell["condition"] == cond, metric].dropna().astype(float)
            enough = len(a) >= 2 and len(b) >= 2
            rows.append({
                "metric": metric, "reference": ref, "condition": cond,
                "n_ref": len(a), "n_cond": len(b),
                "mean_ref": a.mean(), "mean_cond": b.mean(), "sd_ref": a.std(), "sd_cond": b.std(),
                "median_ref": a.median(), "median_cond": b.median(),
                "mannwhitney_p": stats.mannwhitneyu(a, b, alternative="two-sided").pvalue if enough else np.nan,
                "welch_t_p": stats.ttest_ind(a, b, equal_var=False).pvalue if enough else np.nan,
            })
    return pd.DataFrame(rows)


def p_label(p):
    if not np.isfinite(p):
        return "n/a"
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def strip_by_condition(ax, per_cell, metric_by_channel, ylabel, stats_df):
    """Per-cell dots, grouped by condition, one colour per mark (Fig. 6g/h layout)."""
    conds = [c for c in C.CONDITION_ORDER if c in set(per_cell["condition"])]
    rng = np.random.default_rng(0)
    width = 0.35
    for i, cond in enumerate(conds):
        for j, (ch, metric) in enumerate(metric_by_channel.items()):
            vals = per_cell.loc[per_cell["condition"] == cond, metric].dropna().astype(float)
            x = i + (j - (len(metric_by_channel) - 1) / 2) * width
            ax.scatter(x + rng.uniform(-0.06, 0.06, len(vals)), vals, s=18, color=C.CHANNEL_COLORS[ch],
                       label=C.CHANNEL_LABELS.get(ch, ch) if i == 0 else None, zorder=3)
            ax.hlines(vals.mean(), x - 0.1, x + 0.1, color="k", lw=1.5, zorder=4)
            if i > 0:
                s = stats_df[(stats_df["metric"] == metric) & (stats_df["condition"] == cond)]
                if len(s):
                    ax.annotate(p_label(s["mannwhitney_p"].iloc[0]), (x, 1.0), xycoords=("data", "axes fraction"),
                                ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(conds)), conds)
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8)


def mean_sem_curves(ax, df, x, y, group_col, normalize):
    for key, g in df.groupby(group_col):
        piv = g.pivot_table(index=x, columns="cell_id", values=y)
        if normalize:
            piv = piv / piv.max()
        m, sem = piv.mean(axis=1), piv.sem(axis=1)
        ax.plot(m.index, m, color=C.CHANNEL_COLORS[key], label=f"{C.CHANNEL_LABELS[key]} (n={piv.shape[1]})")
        ax.fill_between(m.index, m - sem, m + sem, color=C.CHANNEL_COLORS[key], alpha=0.25, lw=0)


def main():
    parse_args(__doc__.splitlines()[1])
    manifest = load_manifest()
    per_cell, nanodomains, sncr, fits, radial = collect(manifest)
    out = results_dir()
    conds = [c for c in C.CONDITION_ORDER if c in set(per_cell["condition"])]

    per_cell.to_csv(out / "per_cell_summary.csv", index=False)
    nanodomains.to_csv(out / "nanodomains_all.csv", index=False)
    sncr.to_csv(out / "sncr_all.csv", index=False)
    fits.to_csv(out / "sncr_fits_all.csv", index=False)
    radial.to_csv(out / "radial_profiles_all.csv", index=False)

    metrics = ["me3_ac_ratio"] + [f"{ch}_{m}" for ch in C.CHANNELS for m in
                                  ("mean_bb_size_nm", "n_nanodomains", "sncr_decay_rate", "frac_locs_clustered")]
    stats_df = compare_conditions(per_cell, metrics)
    stats_df.to_csv(out / "stats.csv", index=False)

    n_extrap = sum(per_cell[f"{ch}_extrapolated"].sum() for ch in C.CHANNELS)
    if n_extrap:
        log(f"NOTE: {n_extrap} (cell, mark) pairs were clustered outside their calibrated density range")

    # Fig. 6g / 6h
    for fname, metric, ylabel in (("fig6g_nanodomain_size", "mean_bb_size_nm", "Nanodomain size (nm)"),
                                  ("fig6h_nanodomain_count", "n_nanodomains", "Nanodomain count")):
        fig, ax = plt.subplots(figsize=(5, 4.5))
        strip_by_condition(ax, per_cell, {ch: f"{ch}_{metric}" for ch in C.CHANNELS}, ylabel, stats_df)
        fig.tight_layout()
        fig.savefig(out / f"{fname}.png", dpi=200)
        plt.close(fig)

    # Fig. S3a-c
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    strip_by_condition(axes[0], per_cell, {"me3": "me3_ac_ratio"}, "ME:AC localization ratio", stats_df)
    axes[0].get_legend().remove()
    for ax, ch in zip(axes[1:], C.CHANNELS):
        strip_by_condition(ax, per_cell, {ch: f"{ch}_sncr_decay_rate"}, f"{C.CHANNEL_LABELS[ch]} SNCR decay rate (1/nm)",
                           stats_df)
        ax.axhline(0, color="0.5", ls=":", lw=1)
    fig.tight_layout()
    fig.savefig(out / "figS3_ratio_and_decay.png", dpi=200)
    plt.close(fig)

    # Fig. 6b/e and 6c/f
    prim = sncr[sncr["neighbors"] == C.SNCR_PRIMARY]
    for fname, df, x, y, xlabel, ylabel in (
            ("fig6_sncr_by_condition", prim.rename(columns={"ref_channel": "channel"}), "r_nm", C.SNCR_SUMMARY,
             "Distance (nm)", f"Mean SNCR, normalized ({C.SNCR_PRIMARY})"),
            ("fig6_radial_by_condition", radial, "radius_pct", "density_norm",
             "Radius (%)  core → periphery", "Normalized localization density")):
        fig, axes = plt.subplots(1, len(conds), figsize=(5.5 * len(conds), 4.2), sharey=True, squeeze=False)
        for ax, cond in zip(axes[0], conds):
            mean_sem_curves(ax, df[df["condition"] == cond], x, y, "channel", normalize=(x == "r_nm"))
            ax.set(title=cond, xlabel=xlabel)
            ax.legend()
        axes[0][0].set_ylabel(ylabel)
        fig.tight_layout()
        fig.savefig(out / f"{fname}.png", dpi=200)
        plt.close(fig)

    log(f"{len(per_cell)} cells pooled ({per_cell['condition'].value_counts().to_dict()}) -> {out}")
    show = stats_df[["metric", "condition", "n_ref", "n_cond", "mean_ref", "mean_cond", "mannwhitney_p"]]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
