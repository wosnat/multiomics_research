"""Paired RNA/protein FC overview: a global, descriptive look at the data.

For each of the 11 paired observations (organism × condition × timepoint), pull
per-gene log2FC + significance for RNAseq and Proteomics, inner-join on locus_tag
(the paired-detection pool), and report:

  - Pearson r and Spearman rho of (log2FC_RNA vs log2FC_protein) across the
    paired pool.
  - Per-gene category by significance × direction:
      concordant_up      RNA sig up   ∧ prot sig up
      concordant_down    RNA sig down ∧ prot sig down
      rna_only_up        RNA sig up   ∧ prot not sig
      rna_only_down      RNA sig down ∧ prot not sig
      prot_only_up       RNA not sig  ∧ prot sig up
      prot_only_down     RNA not sig  ∧ prot sig down
      opposite_rna_up    RNA sig up   ∧ prot sig down
      opposite_rna_down  RNA sig down ∧ prot sig up
      both_ns            neither significant

  - Compact summary collapsing categories to: concordant_sig, asymmetric_sig,
    opposite_sig, both_ns.

Inputs:
  - data/paired_observations.csv  (from 02_build_paired_observations.py)
Outputs:
  - data/paired_fc_long.csv         long: one row per (paired observation × gene)
  - data/paired_fc_summary.csv      one row per paired observation, with counts and correlations
  - figures/paired_fc_scatter_grid.png   11-panel scatter
  - figures/paired_fc_category_bars.png  stacked bar per observation
  - data/03_paired_fc_overview.log

F1 reminder: MED4 axenic RNAseq has no timepoint metadata in the KG; per
researcher confirmation it maps to 336h (day 14). Handled here by overriding
timepoint_hours on rows from that experiment.

Run from the multiomics_research repo root:
  .venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/03_paired_fc_overview.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from multiomics_explorer import differential_expression_by_gene

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
PAIRED_CSV = DATA_DIR / "paired_observations.csv"
LONG_CSV = DATA_DIR / "paired_fc_long.csv"
SUMMARY_CSV = DATA_DIR / "paired_fc_summary.csv"
LOG_PATH = DATA_DIR / "03_paired_fc_overview.log"
SCATTER_FIG = FIG_DIR / "paired_fc_scatter_grid.png"
BAR_FIG = FIG_DIR / "paired_fc_category_bars.png"

# F1 fix: MED4 axenic RNAseq → 336h (day 14) — see gaps_and_friction.md F1.
F1_EXP_ID = (
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic"
)
F1_TIMEPOINT_HOURS = 336.0


def log(msg: str, *, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def fetch_de(experiment_id: str) -> pd.DataFrame:
    """Pull all DE rows for an experiment (one per gene × timepoint)."""
    result = differential_expression_by_gene(experiment_ids=[experiment_id], limit=None)
    if result.get("truncated"):
        sys.exit(f"ERROR: truncated DE results for {experiment_id}")
    df = pd.DataFrame(result["results"])
    # F1 fix
    if experiment_id == F1_EXP_ID:
        df = df.copy()
        df["timepoint_hours"] = F1_TIMEPOINT_HOURS
    return df


def categorize(rna_status: str, prot_status: str) -> str:
    """Categorize a paired gene by (RNA status, protein status)."""
    rna_sig = rna_status in ("significant_up", "significant_down")
    prot_sig = prot_status in ("significant_up", "significant_down")
    if not rna_sig and not prot_sig:
        return "both_ns"
    if rna_sig and prot_sig:
        if rna_status == prot_status:
            return "concordant_up" if rna_status == "significant_up" else "concordant_down"
        else:
            return (
                "opposite_rna_up"
                if rna_status == "significant_up"
                else "opposite_rna_down"
            )
    if rna_sig:
        return "rna_only_up" if rna_status == "significant_up" else "rna_only_down"
    return "prot_only_up" if prot_status == "significant_up" else "prot_only_down"


CATEGORY_ORDER = [
    "concordant_up", "concordant_down",
    "rna_only_up", "rna_only_down",
    "prot_only_up", "prot_only_down",
    "opposite_rna_up", "opposite_rna_down",
    "both_ns",
]
COMPACT_MAP = {
    "concordant_up": "concordant_sig",
    "concordant_down": "concordant_sig",
    "rna_only_up": "asymmetric_sig",
    "rna_only_down": "asymmetric_sig",
    "prot_only_up": "asymmetric_sig",
    "prot_only_down": "asymmetric_sig",
    "opposite_rna_up": "opposite_sig",
    "opposite_rna_down": "opposite_sig",
    "both_ns": "both_ns",
}
COMPACT_ORDER = ["concordant_sig", "asymmetric_sig", "opposite_sig", "both_ns"]
COMPACT_COLORS = {
    "concordant_sig": "#4C9BCF",
    "asymmetric_sig": "#E89A4F",
    "opposite_sig": "#C0504D",
    "both_ns": "#BDBDBD",
}


def build_long_paired(paired: pd.DataFrame, fh) -> pd.DataFrame:
    """For each paired observation, fetch DE for both omics, inner-join, categorize."""
    rna_cache: dict[str, pd.DataFrame] = {}
    prot_cache: dict[str, pd.DataFrame] = {}

    rows = []
    for _, p in paired.iterrows():
        rna_id = p["rnaseq_experiment_id"]
        prot_id = p["proteomics_experiment_id"]
        tp_h = p["timepoint_hours"]

        if rna_id not in rna_cache:
            rna_cache[rna_id] = fetch_de(rna_id)
        if prot_id not in prot_cache:
            prot_cache[prot_id] = fetch_de(prot_id)

        rna_tp = rna_cache[rna_id]
        prot_tp = prot_cache[prot_id]

        rna_tp = rna_tp[rna_tp["timepoint_hours"] == tp_h].copy()
        prot_tp = prot_tp[prot_tp["timepoint_hours"] == tp_h].copy()

        # Inner-join on locus_tag — paired-detected pool at this TP
        merged = rna_tp.merge(
            prot_tp,
            on="locus_tag",
            suffixes=("_rna", "_prot"),
            how="inner",
        )

        merged["organism_name"] = p["organism_name"]
        merged["condition"] = p["condition"]
        merged["timepoint_hours"] = tp_h
        merged["timepoint_label"] = p["timepoint_label"]
        merged["category"] = [
            categorize(r, q)
            for r, q in zip(merged["expression_status_rna"], merged["expression_status_prot"])
        ]
        merged["category_compact"] = merged["category"].map(COMPACT_MAP)

        # Coalesce gene_name from both omics — RNA-side may have null where
        # protein-side has the annotation, and vice versa. Locus_tag remains the
        # primary identifier throughout the analysis (Rule 2).
        merged["gene_name"] = merged["gene_name_rna"].fillna(merged["gene_name_prot"])
        keep = [
            "organism_name", "condition", "timepoint_label", "timepoint_hours",
            "locus_tag", "gene_name",
            "log2fc_rna", "padj_rna", "expression_status_rna",
            "log2fc_prot", "padj_prot", "expression_status_prot",
            "category", "category_compact",
        ]
        keep = [c for c in keep if c in merged.columns]
        merged = merged[keep]
        rows.append(merged)

        log(
            f"  {p['organism_name']:<32} {p['condition']:<10} "
            f"{p['timepoint_label']:<10} ({int(tp_h):>4}h)  "
            f"paired_pool={len(merged):>4}",
            fh=fh,
        )

    return pd.concat(rows, ignore_index=True)


def per_observation_summary(long: pd.DataFrame) -> pd.DataFrame:
    """Per-observation category counts + correlations."""
    out = []
    for (org, cond, tpl, tph), grp in long.groupby(
        ["organism_name", "condition", "timepoint_label", "timepoint_hours"], sort=False
    ):
        # Drop genes with NaN log2fc in either omic for correlation
        sub = grp.dropna(subset=["log2fc_rna", "log2fc_prot"])
        if len(sub) >= 2:
            pearson_r, pearson_p = stats.pearsonr(sub["log2fc_rna"], sub["log2fc_prot"])
            spearman_r, spearman_p = stats.spearmanr(sub["log2fc_rna"], sub["log2fc_prot"])
        else:
            pearson_r = pearson_p = spearman_r = spearman_p = float("nan")

        cat_counts = grp["category"].value_counts().to_dict()
        compact_counts = grp["category_compact"].value_counts().to_dict()

        row = {
            "organism_name": org,
            "condition": cond,
            "timepoint_label": tpl,
            "timepoint_hours": tph,
            "paired_pool_size": len(grp),
            "n_with_fc": len(sub),
            "pearson_r": pearson_r,
            "pearson_p": pearson_p,
            "spearman_rho": spearman_r,
            "spearman_p": spearman_p,
        }
        for c in CATEGORY_ORDER:
            row[f"n_{c}"] = cat_counts.get(c, 0)
        for c in COMPACT_ORDER:
            row[f"n_{c}"] = compact_counts.get(c, 0)
            row[f"pct_{c}"] = (
                100.0 * compact_counts.get(c, 0) / len(grp) if len(grp) else 0.0
            )
        out.append(row)
    return pd.DataFrame(out)


def make_scatter_grid(long: pd.DataFrame) -> None:
    """11-panel scatter of (log2FC_RNA vs log2FC_protein) per observation."""
    obs_keys = (
        long[["organism_name", "condition", "timepoint_label", "timepoint_hours"]]
        .drop_duplicates()
        .sort_values(["organism_name", "condition", "timepoint_hours"])
        .reset_index(drop=True)
    )
    n = len(obs_keys)
    ncols = 4
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.6, nrows * 3.6),
                             sharex=True, sharey=True)
    axes = axes.flatten()

    for i, (_, k) in enumerate(obs_keys.iterrows()):
        ax = axes[i]
        sub = long[
            (long["organism_name"] == k["organism_name"])
            & (long["condition"] == k["condition"])
            & (long["timepoint_hours"] == k["timepoint_hours"])
        ].dropna(subset=["log2fc_rna", "log2fc_prot"])

        for cat in COMPACT_ORDER:
            ss = sub[sub["category_compact"] == cat]
            ax.scatter(
                ss["log2fc_rna"], ss["log2fc_prot"],
                s=4, alpha=0.5, color=COMPACT_COLORS[cat],
                label=cat if i == 0 else None, rasterized=True,
            )
        ax.axhline(0, color="k", lw=0.5, alpha=0.3)
        ax.axvline(0, color="k", lw=0.5, alpha=0.3)
        # y=x reference
        lim = (-8, 8)
        ax.plot(lim, lim, color="k", lw=0.5, alpha=0.3, linestyle="--")
        ax.set_xlim(lim); ax.set_ylim(lim)

        # Pearson r
        if len(sub) >= 2:
            r, _ = stats.pearsonr(sub["log2fc_rna"], sub["log2fc_prot"])
        else:
            r = float("nan")
        org_short = "MED4" if "MED4" in k["organism_name"] else "HOT1A3"
        ax.set_title(
            f"{org_short} {k['condition']} {k['timepoint_label']}\n"
            f"n={len(sub)}, Pearson r={r:.2f}",
            fontsize=9,
        )
        if i % ncols == 0:
            ax.set_ylabel("log2FC protein")
        if i // ncols == nrows - 1 or i + ncols >= n:
            ax.set_xlabel("log2FC mRNA")

    # Hide unused axes
    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", fontsize=9, frameon=False, ncol=4)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(SCATTER_FIG, dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_category_bars(summary: pd.DataFrame) -> None:
    """Stacked bar of compact categories per observation."""
    summary = summary.sort_values(["organism_name", "condition", "timepoint_hours"]).reset_index(drop=True)
    labels = [
        f"{('MED4' if 'MED4' in r['organism_name'] else 'HOT1A3')} "
        f"{r['condition']} {r['timepoint_label']}"
        for _, r in summary.iterrows()
    ]

    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = np.zeros(len(summary))
    for cat in COMPACT_ORDER:
        vals = summary[f"pct_{cat}"].to_numpy()
        ax.bar(labels, vals, bottom=bottom, color=COMPACT_COLORS[cat], label=cat)
        bottom += vals
    ax.set_ylabel("% of paired-pool genes")
    ax.set_ylim(0, 100)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    ax.set_title("Paired RNA/protein agreement by category, per observation")
    fig.tight_layout()
    fig.savefig(BAR_FIG, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not PAIRED_CSV.exists():
        sys.exit(f"ERROR: {PAIRED_CSV} not found.")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    paired = pd.read_csv(PAIRED_CSV)

    with LOG_PATH.open("w") as fh:
        log(f"Building long paired-FC table for {len(paired)} observations:", fh=fh)

        long = build_long_paired(paired, fh)
        long.to_csv(LONG_CSV, index=False)
        log("", fh=fh)
        log(f"Wrote {LONG_CSV}  ({len(long)} rows)", fh=fh)

        summary = per_observation_summary(long)
        summary.to_csv(SUMMARY_CSV, index=False)
        log(f"Wrote {SUMMARY_CSV}  ({len(summary)} rows)", fh=fh)

        # Pretty-print per-observation summary
        log("", fh=fh)
        log("Per-observation summary:", fh=fh)
        for _, r in summary.iterrows():
            org_short = "MED4" if "MED4" in r["organism_name"] else "HOT1A3"
            log(
                f"  {org_short:<7} {r['condition']:<10} {r['timepoint_label']:<8} "
                f"({int(r['timepoint_hours']):>4}h)  "
                f"n={int(r['paired_pool_size']):>4}  "
                f"Pearson r={r['pearson_r']:>5.2f}  Spearman ρ={r['spearman_rho']:>5.2f}  "
                f"||  conc={int(r['n_concordant_sig']):>3} "
                f"asym={int(r['n_asymmetric_sig']):>4} "
                f"opp={int(r['n_opposite_sig']):>3} "
                f"ns={int(r['n_both_ns']):>4}",
                fh=fh,
            )

        # Compact category percentages
        log("", fh=fh)
        log("Compact category breakdown (% of paired pool):", fh=fh)
        log(
            f"  {'observation':<35} {'concordant_sig':>15} {'asymmetric_sig':>15} "
            f"{'opposite_sig':>15} {'both_ns':>10}",
            fh=fh,
        )
        for _, r in summary.iterrows():
            org_short = "MED4" if "MED4" in r["organism_name"] else "HOT1A3"
            label = f"{org_short} {r['condition']} {r['timepoint_label']}"
            log(
                f"  {label:<35} {r['pct_concordant_sig']:>14.1f}% "
                f"{r['pct_asymmetric_sig']:>14.1f}% "
                f"{r['pct_opposite_sig']:>14.1f}% {r['pct_both_ns']:>9.1f}%",
                fh=fh,
            )

        # Cross-observation overall
        log("", fh=fh)
        log("Cross-observation overall (pooled across all 11):", fh=fh)
        total_n = len(long)
        for cat in COMPACT_ORDER:
            n_cat = (long["category_compact"] == cat).sum()
            log(f"  {cat:<18} = {n_cat:>5} ({100*n_cat/total_n:.1f}% of {total_n} paired-gene-observations)", fh=fh)

        # Figures
        make_scatter_grid(long)
        make_category_bars(summary)
        log("", fh=fh)
        log(f"Wrote {SCATTER_FIG}", fh=fh)
        log(f"Wrote {BAR_FIG}", fh=fh)
        log(f"Wrote {LOG_PATH}", fh=fh)


if __name__ == "__main__":
    main()
