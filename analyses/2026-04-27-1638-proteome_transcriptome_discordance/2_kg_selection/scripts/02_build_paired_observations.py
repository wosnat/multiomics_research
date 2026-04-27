"""Build the paired-observation table: 11 (organism, condition, timepoint) rows
where both RNA-seq and proteomics fold-change data exist at a matched timepoint.

Inputs:
  - data/timepoints_weissberg2025.csv  (from 01_select_experiments.py)

Outputs:
  - data/paired_observations.csv      11 rows: organism × condition × timepoint
  - data/02_build_paired_observations.log

Pairing rules:
  1. Match on (organism, condition, timepoint_hours). condition derived from
     `background_factors`: "axenic" or "coculture".
  2. Drop the pooled "days 60+89" rows (timepoint_hours = NaN). They are a
     statistical-power pooling of d60 and d89; the individual TPs carry the
     same information at finer resolution.
  3. F1 fix — MED4 axenic RNA-seq has `is_time_course=False` and no per-TP row.
     Researcher confirmed (step 1 dialogue) that this single contrast
     corresponds to day 14 (the only timepoint where RNA was extractable from
     axenic MED4; later proteome timepoints had no extractable RNA — see
     gaps_and_friction.md F1, F2). We assert this mapping here.

Run from the multiomics_research repo root:
  .venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/02_build_paired_observations.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
TP_CSV = DATA_DIR / "timepoints_weissberg2025.csv"
OUT_CSV = DATA_DIR / "paired_observations.csv"
LOG_PATH = DATA_DIR / "02_build_paired_observations.log"

# F1 fix: MED4 axenic RNAseq is recorded as a single contrast without a
# timepoint label. Researcher-confirmed pairing (step 1) is to day 14.
MED4_AXENIC_RNASEQ_ID = (
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic"
)
F1_TIMEPOINT_HOURS = 336.0  # day 14 = 336h
F1_TIMEPOINT_LABEL = "day 14"
F1_GROWTH_PHASE = "nutrient_limited"
F1_TP_GENE_COUNT = 1849  # from the experiment's distinct_gene_count


def log(msg: str, *, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def derive_condition(bg: str) -> str | None:
    """Derive 'axenic' / 'coculture' from background_factors string."""
    if not isinstance(bg, str):
        return None
    bg_lower = bg.lower()
    if "axenic" in bg_lower:
        return "axenic"
    if "coculture" in bg_lower:
        return "coculture"
    return None


def main() -> None:
    if not TP_CSV.exists():
        sys.exit(f"ERROR: {TP_CSV} not found. Run 01_select_experiments.py first.")

    tp_df = pd.read_csv(TP_CSV)

    with LOG_PATH.open("w") as fh:
        log(f"Loaded {len(tp_df)} TP rows from {TP_CSV.name}", fh=fh)

        # Derive condition
        tp_df["condition"] = tp_df["background_factors"].map(derive_condition)
        n_no_cond = tp_df["condition"].isna().sum()
        if n_no_cond:
            log(f"WARNING: {n_no_cond} rows have no derivable condition (axenic/coculture)", fh=fh)

        # F1 fix — patch MED4 axenic RNAseq to d14
        f1_mask = tp_df["experiment_id"] == MED4_AXENIC_RNASEQ_ID
        if f1_mask.sum() != 1:
            log(f"WARNING: expected exactly 1 row for MED4 axenic RNAseq, found {f1_mask.sum()}", fh=fh)
        before = tp_df.loc[f1_mask].iloc[0].to_dict() if f1_mask.any() else None
        tp_df.loc[f1_mask, "timepoint"] = F1_TIMEPOINT_LABEL
        tp_df.loc[f1_mask, "timepoint_hours"] = F1_TIMEPOINT_HOURS
        tp_df.loc[f1_mask, "tp_growth_phase"] = F1_GROWTH_PHASE
        tp_df.loc[f1_mask, "tp_gene_count"] = F1_TP_GENE_COUNT
        log("", fh=fh)
        log("F1 fix applied (MED4 axenic RNAseq → day 14):", fh=fh)
        log(f"  Before: {before}", fh=fh)
        if f1_mask.any():
            log(f"  After:  {tp_df.loc[f1_mask].iloc[0].to_dict()}", fh=fh)

        # Drop pooled "days 60+89" rows (timepoint_hours is NaN AFTER the F1 fix)
        before_n = len(tp_df)
        tp_df = tp_df.dropna(subset=["timepoint_hours"])
        log("", fh=fh)
        log(f"Dropped pooled / null-hours rows: {before_n} → {len(tp_df)}", fh=fh)

        # Split into RNAseq and Proteomics; pivot to a paired table
        rna = tp_df[tp_df["omics_type"] == "RNASEQ"].copy()
        prot = tp_df[tp_df["omics_type"] == "PROTEOMICS"].copy()

        rna_keep = [
            "organism_name", "condition", "timepoint_hours", "timepoint",
            "tp_growth_phase", "experiment_id", "tp_gene_count",
        ]
        rna = rna[rna_keep].rename(columns={
            "experiment_id": "rnaseq_experiment_id",
            "tp_gene_count": "rnaseq_tp_gene_count",
        })
        prot_keep = [
            "organism_name", "condition", "timepoint_hours", "timepoint",
            "tp_growth_phase", "experiment_id", "tp_gene_count",
        ]
        prot = prot[prot_keep].rename(columns={
            "timepoint": "timepoint_prot",
            "tp_growth_phase": "tp_growth_phase_prot",
            "experiment_id": "proteomics_experiment_id",
            "tp_gene_count": "proteomics_tp_gene_count",
        })

        # Inner-merge: only timepoints that exist in BOTH omics
        paired = rna.merge(
            prot,
            on=["organism_name", "condition", "timepoint_hours"],
            how="inner",
            validate="one_to_one",
        )

        # Sanity check: per-row, the timepoint label and growth phase agree across omics
        timepoint_mismatches = paired[paired["timepoint"] != paired["timepoint_prot"]]
        if not timepoint_mismatches.empty:
            log("WARNING: timepoint label mismatches between omics:", fh=fh)
            log(timepoint_mismatches.to_string(), fh=fh)
        phase_mismatches = paired[paired["tp_growth_phase"] != paired["tp_growth_phase_prot"]]
        if not phase_mismatches.empty:
            log("WARNING: growth phase mismatches between omics (per matched TP):", fh=fh)
            log(
                phase_mismatches[
                    ["organism_name", "condition", "timepoint_hours", "timepoint",
                     "tp_growth_phase", "tp_growth_phase_prot"]
                ].to_string(),
                fh=fh,
            )

        # Final cleanup: prefer the RNAseq label/phase, drop the duplicate prot cols
        paired = paired.drop(columns=["timepoint_prot", "tp_growth_phase_prot"])
        paired = paired.rename(columns={"timepoint": "timepoint_label"})
        paired = paired.sort_values(
            ["organism_name", "condition", "timepoint_hours"]
        ).reset_index(drop=True)

        ordered = [
            "organism_name", "condition", "timepoint_label", "timepoint_hours",
            "tp_growth_phase",
            "rnaseq_experiment_id", "rnaseq_tp_gene_count",
            "proteomics_experiment_id", "proteomics_tp_gene_count",
        ]
        paired = paired[ordered]
        paired.to_csv(OUT_CSV, index=False)

        log("", fh=fh)
        log(f"Paired observations: {len(paired)} rows", fh=fh)
        for _, r in paired.iterrows():
            log(
                f"  {r['organism_name']:<32} {r['condition']:<10} "
                f"{r['timepoint_label']:<10} ({int(r['timepoint_hours'])}h, "
                f"{r['tp_growth_phase']})  "
                f"rna_genes={int(r['rnaseq_tp_gene_count'])}  "
                f"prot_genes={int(r['proteomics_tp_gene_count'])}",
                fh=fh,
            )

        # Per (organism, condition) summary
        log("", fh=fh)
        log("Pairing coverage by (organism, condition):", fh=fh)
        summary = (
            paired.groupby(["organism_name", "condition"])
            .agg(n_paired_tps=("timepoint_label", "count"))
            .reset_index()
        )
        for _, r in summary.iterrows():
            log(
                f"  {r['organism_name']:<32} {r['condition']:<10} {r['n_paired_tps']} paired TP(s)",
                fh=fh,
            )

        log("", fh=fh)
        log(f"Wrote {OUT_CSV}", fh=fh)
        log(f"Wrote {LOG_PATH}", fh=fh)


if __name__ == "__main__":
    main()
