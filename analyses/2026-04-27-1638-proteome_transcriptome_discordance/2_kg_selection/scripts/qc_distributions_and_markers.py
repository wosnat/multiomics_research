"""QC: log2FC distributions, NaN handling, significance rates, marker sanity.

Four QC checks at step 2:

1. **FC distribution + dynamic range** per (experiment, timepoint).
   median, IQR (q1, q3), |max|, std. Confirms RNAseq has wider dynamic range than
   proteomics (compression effect) and rules out scaling bugs / outlier extremes.

2. **NaN / null log2FC** counts in the paired pool. For each paired observation,
   how many genes have NaN log2fc in either omic? These are dropped from the
   correlation but may carry signal in the categorical view.

3. **Significance rate per TP per omic**, cross-checked against the
   experiment-node `tp_significant_up` / `tp_significant_down` counts already
   present in the timepoints CSV. Catches API-vs-node-metadata drift.

4. **Sense-check on canonical N-stress markers** in MED4. Resolves
   ntcA, glnA, glnB, urtA, amt1, plus rpoB / atpA as housekeeping references,
   via `resolve_gene` to get locus_tags (not by gene_name string-match — many
   genes have null gene_name; locus_tag is the only safe identifier). Pulls
   their FC values across all 5 MED4 paired observations and checks expected
   biological behavior (N markers up in both omics where both detected;
   housekeeping near zero).

Inputs:
  - data/paired_fc_long.csv    (from 03)
  - data/timepoints_weissberg2025.csv  (from 01)
Outputs:
  - data/qc_fc_distributions.csv      one row per (experiment, timepoint)
  - data/qc_paired_nan_summary.csv    one row per paired observation
  - data/qc_significance_crosscheck.csv  per-TP sig counts vs experiment-node
  - data/qc_canonical_markers.csv     one row per (marker × paired observation)
  - figures/qc_fc_distribution_box.png  per-experiment per-TP boxplot
  - figures/qc_canonical_markers.png    heatmap of marker FC across observations
  - data/qc_distributions_and_markers.log

Run from the multiomics_research repo root:
  .venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/qc_distributions_and_markers.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from multiomics_explorer import (
    differential_expression_by_gene,
    resolve_gene,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
LONG_CSV = DATA_DIR / "paired_fc_long.csv"
TP_CSV = DATA_DIR / "timepoints_weissberg2025.csv"
OUT_DIST_CSV = DATA_DIR / "qc_fc_distributions.csv"
OUT_NAN_CSV = DATA_DIR / "qc_paired_nan_summary.csv"
OUT_SIG_CSV = DATA_DIR / "qc_significance_crosscheck.csv"
OUT_MARK_CSV = DATA_DIR / "qc_canonical_markers.csv"
OUT_DIST_FIG = FIG_DIR / "qc_fc_distribution_box.png"
OUT_MARK_FIG = FIG_DIR / "qc_canonical_markers.png"
LOG_PATH = DATA_DIR / "qc_distributions_and_markers.log"

# F1 fix: MED4 axenic RNAseq → 336h (day 14)
F1_EXP_ID = (
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic"
)
F1_TIMEPOINT_HOURS = 336.0

# Marker panel for MED4. Resolved via locus_tag, not gene_name (many genes have
# null gene_name in the KG; locus_tag is the only safe identifier per Rule 2).
MED4_MARKERS = {
    "ntcA":  ("concordance_pos_ctrl", "N regulator"),
    "glnA":  ("concordance_pos_ctrl", "glutamine synthetase"),
    "glnB":  ("concordance_pos_ctrl", "PII protein"),
    "urtA":  ("concordance_pos_ctrl", "urea ABC transporter"),
    "amt1":  ("concordance_pos_ctrl", "ammonium transporter"),
    "rpoB":  ("housekeeping",          "RNA polymerase β"),
    "atpA":  ("housekeeping",          "ATP synthase α (imperfect — flagged in prior analysis)"),
}


def log(msg: str, *, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def fetch_de(experiment_id: str) -> pd.DataFrame:
    result = differential_expression_by_gene(experiment_ids=[experiment_id], limit=None)
    if result.get("truncated"):
        sys.exit(f"ERROR: truncated DE results for {experiment_id}")
    df = pd.DataFrame(result["results"])
    if experiment_id == F1_EXP_ID:
        df = df.copy()
        df["timepoint_hours"] = F1_TIMEPOINT_HOURS
    return df


def qc_fc_distributions(experiment_ids: list[str], fh) -> pd.DataFrame:
    """QC #1: log2FC distribution per (experiment, timepoint)."""
    rows = []
    for eid in experiment_ids:
        df = fetch_de(eid)
        # Group by timepoint_hours; keep NaN-hours rows in their own group
        for tp_h, grp in df.groupby("timepoint_hours", dropna=False):
            n_total = len(grp)
            valid = grp.dropna(subset=["log2fc"])
            n_with = len(valid)
            n_null = n_total - n_with
            if n_with > 0:
                lfc = valid["log2fc"].to_numpy()
                stats_row = {
                    "median": float(np.median(lfc)),
                    "q1": float(np.percentile(lfc, 25)),
                    "q3": float(np.percentile(lfc, 75)),
                    "iqr": float(np.percentile(lfc, 75) - np.percentile(lfc, 25)),
                    "abs_max": float(np.max(np.abs(lfc))),
                    "std": float(np.std(lfc, ddof=1)) if n_with > 1 else float("nan"),
                }
            else:
                stats_row = {k: float("nan") for k in ["median", "q1", "q3", "iqr", "abs_max", "std"]}
            rows.append({
                "experiment_id": eid,
                "timepoint_hours": tp_h,
                "n_total_rows": n_total,
                "n_with_log2fc": n_with,
                "n_null_log2fc": n_null,
                **stats_row,
            })
    return pd.DataFrame(rows)


def qc_paired_nan(long: pd.DataFrame, fh) -> pd.DataFrame:
    """QC #2: how many genes have NaN log2fc in either omic per paired observation?"""
    rows = []
    for (org, cond, tpl, tph), grp in long.groupby(
        ["organism_name", "condition", "timepoint_label", "timepoint_hours"], sort=False
    ):
        n_total = len(grp)
        n_rna_null = grp["log2fc_rna"].isna().sum()
        n_prot_null = grp["log2fc_prot"].isna().sum()
        n_either_null = (grp["log2fc_rna"].isna() | grp["log2fc_prot"].isna()).sum()
        n_both_null = (grp["log2fc_rna"].isna() & grp["log2fc_prot"].isna()).sum()
        rows.append({
            "organism_name": org,
            "condition": cond,
            "timepoint_label": tpl,
            "timepoint_hours": tph,
            "n_paired_pool": n_total,
            "n_rna_null_log2fc": n_rna_null,
            "n_prot_null_log2fc": n_prot_null,
            "n_either_null": n_either_null,
            "n_both_null": n_both_null,
            "pct_either_null": 100.0 * n_either_null / n_total if n_total else float("nan"),
        })
    return pd.DataFrame(rows)


def qc_significance_crosscheck(
    experiment_ids: list[str], tp_df: pd.DataFrame, exp_full: pd.DataFrame, fh
) -> pd.DataFrame:
    """QC #3: per-TP sig counts from row-level data vs experiment-node metadata.

    Falls back to experiment-level totals when per-TP fields are NaN (which is
    expected for single-contrast experiments like MED4 axenic RNA-seq, where
    `is_time_course=False` and per-TP fields are not populated).
    """
    rows = []
    for eid in experiment_ids:
        df = fetch_de(eid)
        for tp_h, grp in df.groupby("timepoint_hours", dropna=False):
            row_sig_up = (grp["expression_status"] == "significant_up").sum()
            row_sig_down = (grp["expression_status"] == "significant_down").sum()

            node = tp_df[
                (tp_df["experiment_id"] == eid)
                & (tp_df["timepoint_hours"].fillna(-1).astype(float) == (tp_h if not pd.isna(tp_h) else -1))
            ]
            node_sig_up: int | None = None
            node_sig_down: int | None = None
            crosscheck_source = "per-TP node"
            if len(node) == 1 and pd.notna(node["tp_significant_up"].iloc[0]):
                node_sig_up = int(node["tp_significant_up"].iloc[0])
                node_sig_down = int(node["tp_significant_down"].iloc[0])
            else:
                # Fallback: experiment-level totals (single-contrast experiments)
                exp = exp_full[exp_full["experiment_id"] == eid]
                if len(exp) == 1 and "genes_by_status_significant_up" in exp.columns:
                    val_up = exp["genes_by_status_significant_up"].iloc[0]
                    val_down = exp["genes_by_status_significant_down"].iloc[0]
                    if pd.notna(val_up):
                        node_sig_up = int(val_up)
                        node_sig_down = int(val_down)
                        crosscheck_source = "experiment-level (single-contrast)"

            row_n = len(grp)
            sig_rate = (row_sig_up + row_sig_down) / row_n if row_n else float("nan")

            match_up = (node_sig_up is not None and node_sig_up == row_sig_up)
            match_down = (node_sig_down is not None and node_sig_down == row_sig_down)

            rows.append({
                "experiment_id": eid,
                "timepoint_hours": tp_h,
                "crosscheck_source": crosscheck_source,
                "row_n": row_n,
                "row_sig_up": int(row_sig_up),
                "row_sig_down": int(row_sig_down),
                "node_sig_up": node_sig_up,
                "node_sig_down": node_sig_down,
                "match_up": match_up,
                "match_down": match_down,
                "sig_rate": sig_rate,
            })
    return pd.DataFrame(rows)


def qc_canonical_markers(long: pd.DataFrame, fh) -> pd.DataFrame:
    """QC #4: resolve N-stress markers to locus_tags and pull their FC across MED4 obs."""
    locus_lookup: dict[str, str] = {}
    for name, (role, product) in MED4_MARKERS.items():
        result = resolve_gene(identifier=name, organism="MED4", limit=10)
        results = result.get("results", [])
        # Pick the MED4 hit
        med4_hits = [r for r in results if "MED4" in r.get("organism_name", "")]
        if not med4_hits:
            log(f"  WARNING: no MED4 hit for marker '{name}' — skipping", fh=fh)
            continue
        if len(med4_hits) > 1:
            log(f"  NOTE: '{name}' resolved to {len(med4_hits)} MED4 hits; using first ({med4_hits[0]['locus_tag']})", fh=fh)
        lt = med4_hits[0]["locus_tag"]
        locus_lookup[name] = lt
        log(f"  resolved {name} → {lt}  ({med4_hits[0].get('product', '?')})", fh=fh)

    # Filter long to MED4 + marker locus_tags
    med4 = long[long["organism_name"].str.contains("MED4")].copy()
    rows = []
    for name, lt in locus_lookup.items():
        sub = med4[med4["locus_tag"] == lt]
        for _, r in sub.iterrows():
            rows.append({
                "marker_name": name,
                "marker_role": MED4_MARKERS[name][0],
                "marker_product": MED4_MARKERS[name][1],
                "locus_tag": lt,
                "condition": r["condition"],
                "timepoint_label": r["timepoint_label"],
                "timepoint_hours": r["timepoint_hours"],
                "log2fc_rna": r["log2fc_rna"],
                "expression_status_rna": r["expression_status_rna"],
                "log2fc_prot": r["log2fc_prot"],
                "expression_status_prot": r["expression_status_prot"],
                "category": r["category"],
            })
    return pd.DataFrame(rows)


def make_distribution_box(dist: pd.DataFrame) -> None:
    """Per-experiment per-TP boxplot of log2FC distribution."""
    dist = dist.copy()
    # Build per-experiment label + omics short
    dist["short"] = dist["experiment_id"].apply(
        lambda s: ("MED4" if "med4" in s else "HOT1A3") + "/"
        + ("RNA" if "rnaseq" in s else "Prot") + "/"
        + ("ax" if "axenic" in s else "cc")
    )
    dist["label"] = dist.apply(
        lambda r: f"{r['short']} {int(r['timepoint_hours']) if pd.notna(r['timepoint_hours']) else 'na'}h",
        axis=1,
    )
    dist = dist.sort_values(["short", "timepoint_hours"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(13, 6))
    # boxplot-style from precomputed q1, median, q3; whiskers at ±1.5*IQR
    for i, r in dist.iterrows():
        if pd.isna(r["median"]):
            continue
        # IQR-based whiskers (capped at observed abs_max)
        whisker_lo = max(r["q1"] - 1.5 * r["iqr"], -r["abs_max"])
        whisker_hi = min(r["q3"] + 1.5 * r["iqr"], r["abs_max"])
        # Box
        ax.add_patch(plt.Rectangle((i - 0.35, r["q1"]), 0.7, r["q3"] - r["q1"],
                                   fill=True, facecolor="#cce5ff",
                                   edgecolor="#1f4e79", linewidth=1))
        # Median
        ax.plot([i - 0.35, i + 0.35], [r["median"], r["median"]],
                color="#1f4e79", linewidth=2)
        # Whiskers
        ax.plot([i, i], [whisker_lo, r["q1"]], color="#1f4e79", linewidth=1)
        ax.plot([i, i], [r["q3"], whisker_hi], color="#1f4e79", linewidth=1)
        ax.plot([i - 0.18, i + 0.18], [whisker_lo, whisker_lo], color="#1f4e79", linewidth=1)
        ax.plot([i - 0.18, i + 0.18], [whisker_hi, whisker_hi], color="#1f4e79", linewidth=1)
        # |max| markers
        ax.scatter([i, i], [-r["abs_max"], r["abs_max"]], color="#c0504d", s=10, zorder=3, marker="x")

    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(len(dist)))
    ax.set_xticklabels(dist["label"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("log2FC")
    ax.set_title("Per-experiment per-TP log2FC distribution\n(box = IQR, line = median, whiskers ±1.5×IQR, × = ±|max|)")
    fig.tight_layout()
    fig.savefig(OUT_DIST_FIG, dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_marker_heatmap(markers: pd.DataFrame) -> None:
    """Marker × paired-observation heatmap of log2FC for both omics."""
    if markers.empty:
        return
    # Build a wide table: rows = marker, cols = (cond, tp), values = log2fc
    markers["obs_label"] = markers.apply(
        lambda r: f"{r['condition']} {r['timepoint_label']}",
        axis=1,
    )
    rna_pivot = markers.pivot(index="marker_name", columns="obs_label", values="log2fc_rna")
    prot_pivot = markers.pivot(index="marker_name", columns="obs_label", values="log2fc_prot")
    # Reorder rows by marker order
    marker_order = [m for m in MED4_MARKERS if m in rna_pivot.index]
    rna_pivot = rna_pivot.reindex(marker_order)
    prot_pivot = prot_pivot.reindex(marker_order)
    # Sort columns by timepoint_hours via lookup
    obs_order = (
        markers[["obs_label", "timepoint_hours"]]
        .drop_duplicates().sort_values("timepoint_hours")["obs_label"].tolist()
    )
    rna_pivot = rna_pivot[obs_order]
    prot_pivot = prot_pivot[obs_order]

    fig, axes = plt.subplots(1, 2, figsize=(10, max(3, 0.5 * len(marker_order))), sharey=True)
    vmax = max(np.nanmax(np.abs(rna_pivot.to_numpy())),
               np.nanmax(np.abs(prot_pivot.to_numpy())))
    for ax, df, title in zip(axes, [rna_pivot, prot_pivot], ["RNA log2FC", "Protein log2FC"]):
        im = ax.imshow(df.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(df.columns)))
        ax.set_xticklabels(df.columns, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(df.index)))
        ax.set_yticklabels(df.index, fontsize=9)
        ax.set_title(title)
        # Annotate cells
        for i in range(len(df.index)):
            for j in range(len(df.columns)):
                val = df.iloc[i, j]
                if pd.notna(val):
                    ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                            fontsize=7, color="k" if abs(val) < vmax * 0.6 else "w")
        fig.colorbar(im, ax=ax, fraction=0.04)
    fig.suptitle("Canonical markers: log2FC across MED4 paired observations")
    fig.tight_layout()
    fig.savefig(OUT_MARK_FIG, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not LONG_CSV.exists():
        sys.exit(f"ERROR: {LONG_CSV} not found. Run 03_paired_fc_overview.py first.")
    if not TP_CSV.exists():
        sys.exit(f"ERROR: {TP_CSV} not found. Run 01_select_experiments.py first.")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    long = pd.read_csv(LONG_CSV)
    tp_df = pd.read_csv(TP_CSV)
    exp_full = pd.read_csv(DATA_DIR / "experiments_weissberg2025_full.csv")

    # Distinct experiment_ids in scope (from the per-TP CSV; 8 total)
    experiment_ids = sorted(tp_df["experiment_id"].unique().tolist())

    with LOG_PATH.open("w") as fh:
        log(f"Step-2 QC for {len(experiment_ids)} experiments, {len(long)} paired-gene rows", fh=fh)

        # === QC #1: FC distributions ===
        log("", fh=fh)
        log("--- QC #1: log2FC distributions per (experiment, timepoint) ---", fh=fh)
        dist = qc_fc_distributions(experiment_ids, fh)
        dist.to_csv(OUT_DIST_CSV, index=False)
        log(f"Wrote {OUT_DIST_CSV}  ({len(dist)} rows)", fh=fh)
        for _, r in dist.iterrows():
            short = ("MED4" if "med4" in r["experiment_id"] else "HOT1A3")
            omic = "RNA " if "rnaseq" in r["experiment_id"] else "Prot"
            cond = "ax" if "axenic" in r["experiment_id"] else "cc"
            tp = int(r["timepoint_hours"]) if pd.notna(r["timepoint_hours"]) else "—"
            log(
                f"  {short:<6} {omic} {cond:<2} tp={tp!s:>5}h  "
                f"n={int(r['n_total_rows']):>5}  null={int(r['n_null_log2fc']):>4}  "
                f"med={r['median']:>+6.2f}  IQR={r['iqr']:>5.2f}  "
                f"|max|={r['abs_max']:>6.2f}  σ={r['std']:>5.2f}",
                fh=fh,
            )

        make_distribution_box(dist)
        log(f"Wrote {OUT_DIST_FIG}", fh=fh)

        # === QC #2: NaN counts in paired pool ===
        log("", fh=fh)
        log("--- QC #2: NaN log2FC in paired pool ---", fh=fh)
        nan_summary = qc_paired_nan(long, fh)
        nan_summary.to_csv(OUT_NAN_CSV, index=False)
        log(f"Wrote {OUT_NAN_CSV}  ({len(nan_summary)} rows)", fh=fh)
        for _, r in nan_summary.iterrows():
            org_short = "MED4" if "MED4" in r["organism_name"] else "HOT1A3"
            log(
                f"  {org_short:<6} {r['condition']:<10} {r['timepoint_label']:<10}  "
                f"n={int(r['n_paired_pool']):>4}  rna_null={int(r['n_rna_null_log2fc']):>3}  "
                f"prot_null={int(r['n_prot_null_log2fc']):>3}  "
                f"either={int(r['n_either_null']):>3} ({r['pct_either_null']:.1f}%)  "
                f"both={int(r['n_both_null']):>3}",
                fh=fh,
            )

        # === QC #3: significance crosscheck ===
        log("", fh=fh)
        log("--- QC #3: per-TP significance crosscheck (rows vs experiment-node) ---", fh=fh)
        sig = qc_significance_crosscheck(experiment_ids, tp_df, exp_full, fh)
        sig.to_csv(OUT_SIG_CSV, index=False)
        log(f"Wrote {OUT_SIG_CSV}  ({len(sig)} rows)", fh=fh)
        any_mismatch = False
        for _, r in sig.iterrows():
            tp = int(r["timepoint_hours"]) if pd.notna(r["timepoint_hours"]) else "—"
            ok = "✓" if (r["match_up"] and r["match_down"]) else "✗"
            if not (r["match_up"] and r["match_down"]):
                any_mismatch = True
            log(
                f"  {ok}  {r['experiment_id'][-50:]:<50}  tp={tp!s:>5}h  "
                f"row_up={int(r['row_sig_up']):>4} (node {r['node_sig_up']!s:>4})  "
                f"row_dn={int(r['row_sig_down']):>4} (node {r['node_sig_down']!s:>4})  "
                f"sig_rate={r['sig_rate']:.1%}",
                fh=fh,
            )
        log(f"  {'PASS' if not any_mismatch else 'WARN: mismatches detected'} (sig counts agree row-level vs experiment-node)", fh=fh)

        # === QC #4: canonical markers ===
        log("", fh=fh)
        log("--- QC #4: canonical N-stress markers (MED4) ---", fh=fh)
        markers = qc_canonical_markers(long, fh)
        markers.to_csv(OUT_MARK_CSV, index=False)
        log(f"Wrote {OUT_MARK_CSV}  ({len(markers)} rows)", fh=fh)

        if not markers.empty:
            log("", fh=fh)
            log(f"  {'marker':<6} {'role':<22} {'condition':<10} {'TP':<8}  "
                f"{'rna lfc':>8} {'rna sig':<18}  {'prot lfc':>8} {'prot sig':<18}  cat", fh=fh)
            for _, r in markers.iterrows():
                lfc_rna = f"{r['log2fc_rna']:+.2f}" if pd.notna(r['log2fc_rna']) else "  NaN"
                lfc_prot = f"{r['log2fc_prot']:+.2f}" if pd.notna(r['log2fc_prot']) else "  NaN"
                log(
                    f"  {r['marker_name']:<6} {r['marker_role']:<22} {r['condition']:<10} "
                    f"{r['timepoint_label']:<8}  {lfc_rna:>8} {r['expression_status_rna']:<18}  "
                    f"{lfc_prot:>8} {r['expression_status_prot']:<18}  {r['category']}",
                    fh=fh,
                )

            make_marker_heatmap(markers)
            log(f"Wrote {OUT_MARK_FIG}", fh=fh)

        log("", fh=fh)
        log(f"Wrote {LOG_PATH}", fh=fh)


if __name__ == "__main__":
    main()
