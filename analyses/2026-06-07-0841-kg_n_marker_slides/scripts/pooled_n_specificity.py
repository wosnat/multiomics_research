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

    res_o = res.set_index("family").loc[ORDER].reset_index()

    # ---- (1) the table -> CSV (slide-ready, friendly columns) ----
    tbl = pd.DataFrame({
        "family": res_o["family"],
        "description": res_o["description"],
        "up_per_tested_N": res_o["up_N"].astype(str) + "/" + res_o["tested_N"].astype(str),
        "up_per_tested_other": res_o["up_nonN"].astype(str) + "/" + res_o["tested_nonN"].astype(str),
        "odds_ratio": res_o["odds_ratio"].round(1),
        "q_BH": res_o["p_bh"].round(4),
        "significant_q<0.05": res_o["p_bh"] < 0.05,
    })
    tbl.to_csv(os.path.join(DATA, "fig8_table.csv"), index=False)

    # ---- (2) forest plot, large fonts for PPT (no table; OR labelled on points) ----
    plt.rcParams.update({"font.size": 18, "svg.fonttype": "none"})
    RED, GREY = "#d73027", "#9aa3a8"
    S = 1.6  # row spacing
    ys = (np.arange(len(res_o)) * S)[::-1]
    fig, ax = plt.subplots(figsize=(13, 6.5))
    for y, r in zip(ys, res_o.itertuples()):
        sig = r.p_bh < 0.05
        col = RED if sig else GREY
        ax.plot([max(r.ci_low, 0.05), r.ci_high], [y, y], color=col, lw=5, zorder=2,
                solid_capstyle="round")
        ax.scatter([r.or_haldane], [y], s=420, color=col, zorder=3,
                   edgecolor="white", linewidth=2)
        ax.text(0.05, y + 0.34, r.family, fontsize=30, weight="bold",
                va="bottom", ha="left", transform=ax.get_yaxis_transform())
        ax.text(0.05, y - 0.30, FAMILY_DESC[r.family], fontsize=17, color="#666",
                va="top", ha="left", transform=ax.get_yaxis_transform())
        star = " *" if sig else ""
        or_txt = f"OR {r.odds_ratio:.1f}{star}" if np.isfinite(r.odds_ratio) else f"OR ∞{star}"
        ax.text(r.ci_high * 1.4, y, or_txt, fontsize=21, va="center", ha="left",
                color=col, weight="bold")

    ax.axvline(1, ls="--", color="#444", lw=2, zorder=1)
    ax.set_xscale("log"); ax.set_xlim(0.05, 9000)
    ax.xaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1000]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.1", "1", "10", "100", "1000"]))
    ax.tick_params(axis="x", labelsize=20, length=7, width=1.5)
    ax.set_xlabel("odds ratio   (up-regulation enriched under N  →)", fontsize=22)
    ytop = (len(res_o) - 1) * S
    ax.set_yticks([]); ax.set_ylim(-0.95, ytop + 1.05)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)

    ax.set_title("Is up-regulation specific to nitrogen stress?",
                 fontsize=27, weight="bold", pad=26)
    fig.text(0.5, 0.015,
             "Pooled across genomes · Fisher exact on DE/not-DE calls (no fold-change) · "
             "red = significant (q<0.05) · whiskers = 95% CI",
             fontsize=15, color="#555", ha="center")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIGS, f"fig8_n_specificity_forest.{ext}"),
                    dpi=200, bbox_inches="tight")
    print("\nwrote figures/fig8_n_specificity_forest.{png,svg}, "
          "data/fig8_table.csv, data/n_specificity_pooled_results.csv, data/pooled_de_downloaded.csv")


if __name__ == "__main__":
    main()
