#!/usr/bin/env python
"""
Is up-regulation of the N-acquisition gene families SPECIFIC to nitrogen stress?
POOLED across all genomes that have DE data, for three families: urt*, amt*, cynS.

Self-contained: downloads the data from the KG via the multiomics_explorer
Python API (differential_expression_by_ortholog -> to_dataframe), caches it to
data/pooled_de_downloaded.csv, then computes the test and the figure. Rerun with:

    uv run python analyses/.../scripts/pooled_n_specificity.py

Test (per family): Fisher's exact on a 2x2 of
    (experiment is N-stress vs not) x (family significantly UP vs not-up)
- Observation unit = experiment (each experiment is one genome; timepoints +
  family subunits collapsed: "up" if ANY member gene is significantly up in ANY
  timepoint; "tested" if any member is detected in that experiment).
- Uses the DE / not-DE call only -- NO fold-change magnitude -- so it is robust to
  platform dynamic-range differences (RNA-seq vs microarray vs proteomics).
- One-sided (greater) = up enriched under N; BH-FDR across the three families.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

from multiomics_explorer import differential_expression_by_ortholog as deo, to_dataframe

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "figures")

# Three families -> their curated cyanorak ortholog groups (cross-genus Pro+Syn).
FAMILIES = {
    "urt*": ["cyanorak:CK_00000076", "cyanorak:CK_00001365", "cyanorak:CK_00001366",
             "cyanorak:CK_00001367", "cyanorak:CK_00008074"],   # urtA-E (urea ABC transporter)
    "amt*": ["cyanorak:CK_00000244", "cyanorak:CK_00008701"],   # amt1, amt2 (ammonium transporter)
    "cynS": ["cyanorak:CK_00001552"],                           # cyanate hydratase
}
FAMILY_DESC = {"urt*": "urea transporter (urtA-E)",
               "amt*": "ammonium transporter (amt1/2)",
               "cynS": "cyanate hydratase"}
ORDER = ["urt*", "amt*", "cynS"]


def download():
    """Pull all DE-by-ortholog rows (all genomes, incl. not-significant) for every group."""
    all_groups = [g for gs in FAMILIES.values() for g in gs]
    res = deo(group_ids=all_groups, significant_only=False)  # package default limit=None -> all rows
    df = to_dataframe(res)
    df.to_csv(os.path.join(DATA, "pooled_de_downloaded.csv"), index=False)
    print(f"downloaded {len(df)} rows across {df['organism_name'].nunique()} genomes "
          f"and {df['experiment_id'].nunique()} experiments")
    return df


def is_nitrogen(tt) -> bool:
    return "nitrogen" in str(tt).split("|")


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def fisher_family(df_fam):
    """Collapse to one row per experiment, then 2x2 Fisher (up enriched under N)."""
    g = (df_fam.assign(is_N=df_fam["treatment_type"].map(is_nitrogen),
                       up=df_fam["significant_up"] > 0)
               .groupby("experiment_id")
               .agg(up=("up", "max"), is_N=("is_N", "max")).reset_index())
    N, O = g[g.is_N], g[~g.is_N]
    a, c = int(N.up.sum()), len(N) - int(N.up.sum())   # N: up, not-up
    b, d = int(O.up.sum()), len(O) - int(O.up.sum())   # non-N: up, not-up
    _, p_greater = fisher_exact([[a, c], [b, d]], alternative="greater")
    or_raw, _ = fisher_exact([[a, c], [b, d]], alternative="two-sided")
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_h = (a_ * d_) / (b_ * c_)
    se = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    return dict(up_N=a, tested_N=len(N), up_nonN=b, tested_nonN=len(O),
                odds_ratio=or_raw, or_haldane=or_h,
                ci_low=np.exp(np.log(or_h) - 1.96 * se),
                ci_high=np.exp(np.log(or_h) + 1.96 * se),
                p_greater=p_greater)


def main():
    df = download()
    g2fam = {g: fam for fam, gs in FAMILIES.items() for g in gs}
    df["family"] = df["group_id"].map(g2fam)

    rows = []
    for fam in ORDER:
        sub = df[df.family == fam]
        rows.append({"family": fam, "description": FAMILY_DESC[fam], **fisher_family(sub)})
    res = pd.DataFrame(rows)
    res["p_bh"] = bh_fdr(res["p_greater"].values)
    res.to_csv(os.path.join(DATA, "n_specificity_pooled_results.csv"), index=False)
    print(res[["family", "up_N", "tested_N", "up_nonN", "tested_nonN",
               "odds_ratio", "p_greater", "p_bh"]].to_string(index=False))

    # ---- forest plot (3 rows) ----
    RED, GREY = "#d73027", "#9aa3a8"
    res_o = res.set_index("family").loc[ORDER].reset_index()
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.4),
                                   gridspec_kw={"width_ratios": [7, 5]})
    S = 1.6  # row spacing
    ys = (np.arange(len(res_o)) * S)[::-1]
    for y, r in zip(ys, res_o.itertuples()):
        sig = r.p_bh < 0.05
        col = RED if sig else GREY
        axL.plot([max(r.ci_low, 0.05), r.ci_high], [y, y], color=col, lw=2.5, zorder=2)
        axL.scatter([r.or_haldane], [y], s=110, color=col, zorder=3,
                    edgecolor="white", linewidth=1.3)
        axL.text(0.052, y + 0.30, f"{r.family}", fontsize=12.5, weight="bold",
                 va="bottom", ha="left", transform=axL.get_yaxis_transform())
        axL.text(0.052, y - 0.30, FAMILY_DESC[r.family], fontsize=8.5,
                 color="#666", va="top", ha="left", transform=axL.get_yaxis_transform())
        star = "  *" if sig else ""
        axR.text(0.00, y, f"{r.up_N}/{r.tested_N}", fontsize=10.5, va="center", ha="center")
        axR.text(0.22, y, f"{r.up_nonN}/{r.tested_nonN}", fontsize=10.5, va="center", ha="center")
        axR.text(0.46, y, f"{r.odds_ratio:.1f}" if np.isfinite(r.odds_ratio) else "∞",
                 fontsize=10.5, va="center", ha="center")
        axR.text(0.74, y, f"{r.p_bh:.3g}{star}", fontsize=10.5, va="center", ha="center",
                 color=RED if sig else "#333", weight="bold" if sig else "normal")

    axL.axvline(1, ls="--", color="#444", lw=1, zorder=1)
    axL.set_xscale("log"); axL.set_xlim(0.05, 3000)
    axL.xaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1000]))
    axL.xaxis.set_major_formatter(FixedFormatter(["0.1", "1", "10", "100", "1000"]))
    axL.set_xlabel("odds ratio   (up-regulation enriched under N  →)", fontsize=10.5)
    ytop = (len(res_o) - 1) * S
    axL.set_yticks([]); axL.set_ylim(-0.95, ytop + 1.0)
    for s in ("top", "right", "left"):
        axL.spines[s].set_visible(False)
    axR.set_xlim(-0.1, 0.95); axR.set_ylim(-0.95, ytop + 1.0); axR.axis("off")
    hy = ytop + 0.85
    for x, lab in [(0.0, "up/tested\n(N)"), (0.22, "up/tested\n(other)"),
                   (0.46, "odds\nratio"), (0.74, "q (BH)")]:
        axR.text(x, hy, lab, fontsize=9, ha="center", va="center", weight="bold", color="#555")

    fig.suptitle("Is up-regulation specific to nitrogen stress?  (pooled across genomes)\n"
                 "Fisher exact on DE / not-DE calls — fold-change magnitude not used",
                 fontsize=12.5, weight="bold")
    fig.text(0.5, 0.01,
             "OR > 1 (right of dashed line) = up-regulation more frequent under N than other conditions.   "
             "Red = BH q < 0.05.   Whiskers = 95% CI.", fontsize=8.5, color="#555", ha="center")
    fig.tight_layout(rect=[0, 0.04, 1, 0.88])
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIGS, f"fig8_n_specificity_forest.{ext}"),
                    dpi=200, bbox_inches="tight")
    print("\nwrote figures/fig8_n_specificity_forest.{png,svg}, "
          "data/n_specificity_pooled_results.csv, data/pooled_de_downloaded.csv")


if __name__ == "__main__":
    main()
