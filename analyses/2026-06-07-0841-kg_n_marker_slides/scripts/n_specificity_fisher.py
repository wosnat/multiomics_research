#!/usr/bin/env python
"""
Is up-regulation of the N-acquisition genes SPECIFIC to nitrogen stress?

Test (per gene): Fisher's exact on a 2x2 contingency of
    (experiment is N-stress vs not)  x  (gene significantly UP vs not-up)
Unit of observation = experiment (timepoints collapsed: "up" if the gene is
significantly up in >=1 timepoint of that experiment). "Tested" = the gene is
DETECTED (has a DE row) in that experiment. Fold-change MAGNITUDE is never used
-- only the DE/not-DE call -- so the test is robust to platform dynamic-range
differences (RNA-seq vs microarray vs proteomics).

Scopes:
  - MED4 : per-gene panel (cynS, cynD, urtA-E, amt1) within Prochlorococcus MED4
  - POOLED: gene-family level (urtA, cynS ortholog groups) pooled across all
            genomes that have DE data; observation = (genome x experiment).

Inputs (staged from the KG, frozen):
  data/med4_obs.csv     gene_name, locus_tag, experiment_id, treatment_types, expression_status
  data/pooled_obs.csv   group_id, gene, organism, experiment_id, treatment_types, significant_up/down, not_significant
Outputs:
  data/n_specificity_results.csv
  figures/fig8_n_specificity_forest.{png,svg}
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "figures")

PANEL_ORDER = ["cynS", "cynD", "amt1", "urtA", "urtB", "urtC", "urtD", "urtE"]


def is_nitrogen(tt: str) -> bool:
    return "nitrogen" in str(tt).split("|")


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downward
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def fisher_row(label, up_N, tested_N, up_O, tested_O):
    """2x2: rows = [N, non-N], cols = [up, not-up]. One-sided 'greater' = up enriched under N."""
    a, c = up_N, tested_N - up_N           # N: up, not-up
    b, d = up_O, tested_O - up_O           # non-N: up, not-up
    table = [[a, c], [b, d]]
    # one-sided greater test for up-enrichment under N
    _, p_greater = fisher_exact(table, alternative="greater")
    or_raw, _ = fisher_exact(table, alternative="two-sided")  # sample OR
    # Haldane-Anscombe corrected OR + 95% CI (handles zero cells, used for plotting)
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_h = (a_ * d_) / (b_ * c_)
    se = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    lo, hi = np.exp(np.log(or_h) - 1.96 * se), np.exp(np.log(or_h) + 1.96 * se)
    return dict(
        label=label, up_N=a, tested_N=tested_N, up_nonN=b, tested_nonN=tested_O,
        rate_N=a / tested_N if tested_N else np.nan,
        rate_nonN=b / tested_O if tested_O else np.nan,
        odds_ratio=or_raw, or_haldane=or_h, ci_low=lo, ci_high=hi,
        p_greater=p_greater,
    )


# ---------- MED4 per-gene ----------
med4 = pd.read_csv(os.path.join(DATA, "med4_obs.csv"))
med4["is_N"] = med4["treatment_types"].map(is_nitrogen)
med4["up"] = med4["expression_status"].eq("significant_up")
# collapse timepoints -> per (gene, experiment): up if up in any timepoint; experiment is N if any row N
g = (med4.groupby(["gene_name", "experiment_id"])
          .agg(up=("up", "max"), is_N=("is_N", "max")).reset_index())

med4_rows = []
for gene in PANEL_ORDER:
    sub = g[g["gene_name"] == gene]
    if sub.empty:
        continue
    N, O = sub[sub.is_N], sub[~sub.is_N]
    med4_rows.append({"scope": "MED4", **fisher_row(
        gene, int(N.up.sum()), len(N), int(O.up.sum()), len(O))})

# ---------- POOLED family level ----------
pool = pd.read_csv(os.path.join(DATA, "pooled_obs.csv"))
pool["is_N"] = pool["treatment_types"].map(is_nitrogen)
pool["up"] = pool["significant_up"] > 0
# collapse timepoints -> per (group, organism, experiment): up if any member up in any timepoint
gp = (pool.groupby(["group_id", "gene", "organism", "experiment_id"])
          .agg(up=("up", "max"), is_N=("is_N", "max")).reset_index())

pool_rows = []
for gid, name in [("cyanorak:CK_00000076", "urtA family"),
                  ("cyanorak:CK_00001552", "cynS family")]:
    sub = gp[gp["group_id"] == gid]
    N, O = sub[sub.is_N], sub[~sub.is_N]
    pool_rows.append({"scope": "POOLED", **fisher_row(
        name, int(N.up.sum()), len(N), int(O.up.sum()), len(O))})

# ---------- combine + BH within each scope ----------
res = pd.DataFrame(med4_rows + pool_rows)
for scope in res["scope"].unique():
    m = res["scope"] == scope
    res.loc[m, "p_bh"] = bh_fdr(res.loc[m, "p_greater"].values)
res.to_csv(os.path.join(DATA, "n_specificity_results.csv"), index=False)
print(res[["scope", "label", "up_N", "tested_N", "up_nonN", "tested_nonN",
           "odds_ratio", "p_greater", "p_bh"]].to_string(index=False))

# ---------- forest plot (single stacked axis + text column) ----------
RED, GREY = "#d73027", "#9aa3a8"
med4_df = res[res.scope == "MED4"].set_index("label").loc[PANEL_ORDER].reset_index()
pool_df = res[res.scope == "POOLED"].reset_index(drop=True)

# build display rows top->bottom: header, MED4 genes, header, pooled families
disp = [("header", "MED4  —  per gene")]
disp += [("row", r) for r in med4_df.to_dict("records")]
disp += [("gap", None), ("header", "Pooled  —  gene family (all genomes)")]
disp += [("row", r) for r in pool_df.to_dict("records")]

n = len(disp)
ys = {i: n - i for i in range(n)}  # top row highest y

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 6),
                               gridspec_kw={"width_ratios": [7, 5]})
for i, (kind, r) in enumerate(disp):
    y = ys[i]
    if kind == "header":
        axL.text(0.035, y, r, fontsize=11.5, weight="bold", va="center", ha="left")
        continue
    if kind == "gap":
        continue
    sig = r["p_bh"] < 0.05
    col = RED if sig else GREY
    lo = max(r["ci_low"], 0.045)
    axL.plot([lo, r["ci_high"]], [y, y], color=col, lw=2.2, zorder=2)
    axL.scatter([r["or_haldane"]], [y], s=85, color=col, zorder=3,
                edgecolor="white", linewidth=1.2)
    axL.text(0.06, y, "   " + r["label"], fontsize=11, va="center", ha="left")
    # text column: rates + OR + q
    star = "  *" if sig else ""
    axR.text(0.00, y, f"{r['up_N']}/{r['tested_N']}", fontsize=10, va="center", ha="center")
    axR.text(0.18, y, f"{r['up_nonN']}/{r['tested_nonN']}", fontsize=10, va="center", ha="center")
    axR.text(0.40, y, f"{r['odds_ratio']:.1f}" if np.isfinite(r['odds_ratio']) else "∞",
             fontsize=10, va="center", ha="center")
    axR.text(0.70, y, f"{r['p_bh']:.3g}{star}", fontsize=10, va="center", ha="center",
             color=RED if sig else "#333", weight="bold" if sig else "normal")

# left axis: log-OR forest
axL.axvline(1, ls="--", color="#444", lw=1, zorder=1)
axL.set_xscale("log")
axL.set_xlim(0.045, 3000)
axL.xaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1000]))
axL.xaxis.set_major_formatter(FixedFormatter(["0.1", "1", "10", "100", "1000"]))
axL.set_xlabel("odds ratio   (up-regulation enriched under N  →)", fontsize=10.5)
axL.set_ylim(0.3, n + 0.7)
axL.set_yticks([])
for s in ("top", "right", "left"):
    axL.spines[s].set_visible(False)

# right axis: text table header
axR.set_xlim(-0.08, 0.92)
axR.set_ylim(0.3, n + 0.7)
axR.axis("off")
hdr_y = n + 0.6
axR.text(0.00, hdr_y, "up/​tested\n(N)", fontsize=9, ha="center", va="center", weight="bold", color="#555")
axR.text(0.18, hdr_y, "up/​tested\n(other)", fontsize=9, ha="center", va="center", weight="bold", color="#555")
axR.text(0.40, hdr_y, "odds\nratio", fontsize=9, ha="center", va="center", weight="bold", color="#555")
axR.text(0.70, hdr_y, "q (BH)", fontsize=9, ha="center", va="center", weight="bold", color="#555")

fig.suptitle("Is up-regulation specific to nitrogen stress?\n"
             "Fisher exact on DE / not-DE calls (fold-change magnitude not used)",
             fontsize=13.5, weight="bold")
fig.text(0.5, 0.015,
         "OR > 1 (right of dashed line) = significant up-regulation is more frequent under N than under other conditions.   "
         "Red = BH q < 0.05.   Whiskers = 95% CI (Haldane-corrected).",
         fontsize=8.5, color="#555", ha="center")
fig.tight_layout(rect=[0, 0.03, 1, 0.93])
for ext in ("png", "svg"):
    fig.savefig(os.path.join(FIGS, f"fig8_n_specificity_forest.{ext}"),
                dpi=200, bbox_inches="tight")
print("\nwrote figures/fig8_n_specificity_forest.{png,svg} and data/n_specificity_results.csv")
