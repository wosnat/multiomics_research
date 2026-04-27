"""Pull DE values for the step-3 control panel from the 5 Pro MED4 experiments.

Inputs:
  - 3_analysis_framing/data/control_panel.csv  (26 genes: 23 positives + 3 negatives)
  - 2_kg_selection/data/experiments_pro_med4.csv (5 experiment IDs)

Outputs:
  - data/control_de_long.csv         one row per (control gene, experiment, timepoint), with log2FC, padj, expression_status, omics, condition, growth_phase
  - data/control_de_summary.csv      one row per control gene (max |log2FC|, n_significant_TPs, etc.)
  - data/01_pull_control_de.log      filter funnel + diagnostics

Run from the multiomics_research repo root:
  .venv/bin/python analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/3_analysis_framing/scripts/01_pull_control_de.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe

ANALYSIS_ROOT = Path(__file__).resolve().parents[2]
CONTROL_CSV = ANALYSIS_ROOT / "3_analysis_framing" / "data" / "control_panel.csv"
EXP_CSV = ANALYSIS_ROOT / "2_kg_selection" / "data" / "experiments_pro_med4.csv"
TP_CSV = ANALYSIS_ROOT / "2_kg_selection" / "data" / "timepoints_pro_med4.csv"

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = OUT_DIR / "01_pull_control_de.log"


def log(msg: str, *, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def main() -> None:
    with LOG_PATH.open("w") as fh:
        controls = pd.read_csv(CONTROL_CSV)
        log(f"Control panel: {len(controls)} genes across "
            f"{controls['axis'].nunique()} axes "
            f"({(controls['role']=='positive').sum()} positives, "
            f"{(controls['role']=='negative').sum()} negatives)", fh=fh)

        experiments = pd.read_csv(EXP_CSV)
        # Drop the cross-condition contrast experiment from validation (uninformative for stress)
        trajectory_exp_ids = experiments[
            ~experiments["experiment_id"].str.contains("coculture_alteromonas_hot1a3_med4")
        ]["experiment_id"].tolist()
        log(f"Trajectory + summary experiments: {len(trajectory_exp_ids)}", fh=fh)
        for eid in trajectory_exp_ids:
            log(f"  - {eid}", fh=fh)

        locus_tags = controls["locus_tag"].tolist()

        # Single API call per experiment — restrict to control genes
        rows: list[dict] = []
        for eid in trajectory_exp_ids:
            res = differential_expression_by_gene(
                experiment_ids=[eid],
                locus_tags=locus_tags,
                limit=None,
                verbose=True,
            )
            df = to_dataframe(res)
            df["experiment_id"] = eid
            log(f"\n[{eid}]", fh=fh)
            log(f"  total_matching={res['total_matching']} returned={len(df)}", fh=fh)
            if len(df):
                log(f"  columns={list(df.columns)[:10]}...", fh=fh)
            rows.append(df)

        all_de = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        log(f"\nCombined DE rows: {len(all_de)}", fh=fh)

        # Join control panel metadata (axis, role) onto DE rows
        de = all_de.merge(
            controls[["locus_tag", "axis", "role", "kg_evidence"]],
            on="locus_tag",
            how="left",
        )

        # NOTE: differential_expression_by_gene already returns a per-row `growth_phase`
        # column, so we don't need the F3 workaround here. Verified that DE-row growth_phase
        # agrees with the manually-zipped tp_growth_phase from 2_kg_selection (when both
        # are present); the only NAs are on the axenic-RNA single-point experiment, where
        # DE-row growth_phase is "nutrient_limited" matching the experiment metadata.

        # Drop the pooled "days 60+89" rows per step 2 decision
        if "timepoint" in de.columns:
            before = len(de)
            de = de[de["timepoint"].astype(str) != "days 60+89"].copy()
            log(f"\nDropped pooled 'days 60+89' rows: {before} -> {len(de)}", fh=fh)

        # Write long table
        de_path = OUT_DIR / "control_de_long.csv"
        de.to_csv(de_path, index=False)
        log(f"Wrote {de_path}", fh=fh)
        log(f"  shape={de.shape}", fh=fh)
        log(f"  columns={list(de.columns)}", fh=fh)

        # Per-gene summary: max |log2FC|, number of significant TPs, axes
        summary_rows = []
        for locus, group in de.groupby("locus_tag"):
            gene_meta = controls[controls["locus_tag"] == locus].iloc[0]
            log2fc = pd.to_numeric(group["log2fc"], errors="coerce") if "log2fc" in group.columns else pd.Series(dtype=float)
            sig = group["expression_status"] if "expression_status" in group.columns else pd.Series([], dtype=object)
            n_total = len(group)
            sig_str = sig.astype(str)
            n_sig = (sig_str.str.contains("significant", case=False, na=False) & ~sig_str.str.contains("not_significant", case=False, na=False)).sum()
            has_any = log2fc.notna().any()
            summary_rows.append({
                "locus_tag": locus,
                "gene_name": gene_meta["gene_name"],
                "axis": gene_meta["axis"],
                "role": gene_meta["role"],
                "n_de_rows": n_total,
                "n_significant": int(n_sig),
                "max_abs_log2fc": float(log2fc.abs().max()) if has_any else float("nan"),
                "max_log2fc": float(log2fc.max()) if has_any else float("nan"),
                "min_log2fc": float(log2fc.min()) if has_any else float("nan"),
                "n_rows_with_log2fc": int(log2fc.notna().sum()),
            })
        summary = pd.DataFrame(summary_rows).sort_values(["axis", "role", "locus_tag"])
        summary_path = OUT_DIR / "control_de_summary.csv"
        summary.to_csv(summary_path, index=False)
        log(f"\nWrote {summary_path}", fh=fh)

        # Diagnostic: per-axis significance rate
        log("\nPer-axis significance rate (positives only):", fh=fh)
        pos = summary[summary["role"] == "positive"]
        for axis, group in pos.groupby("axis"):
            n_genes = len(group)
            n_with_sig = (group["n_significant"] > 0).sum()
            mean_max_abs = group["max_abs_log2fc"].mean()
            log(
                f"  {axis:<14} n={n_genes:>2}  with_any_sig={n_with_sig:>2}/{n_genes}  "
                f"mean_max|log2FC|={mean_max_abs:.2f}",
                fh=fh,
            )

        log("\nNegative controls (should show low DE):", fh=fh)
        neg = summary[summary["role"] == "negative"]
        for _, row in neg.iterrows():
            log(
                f"  {row['gene_name']:<6} ({row['locus_tag']}): "
                f"n_sig={row['n_significant']}/{row['n_de_rows']}  "
                f"max|log2FC|={row['max_abs_log2fc']:.2f}",
                fh=fh,
            )


if __name__ == "__main__":
    main()
