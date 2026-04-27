"""Select Weissberg 2025 paired RNA-seq + proteomics experiments.

Inputs:  none (queries the multiomics KG via the Python API)
Outputs:
  - data/experiments_weissberg2025.csv       compact view of the 8 selected experiments
  - data/experiments_weissberg2025_full.csv  full metadata (debug / reference)
  - data/timepoints_weissberg2025.csv        per-(experiment, timepoint), enriched with per-TP growth phase
  - data/01_select_experiments.log           filter funnel + per-experiment TP listing

Selection criteria for this analysis (from step 1):
  - publication_doi == "10.1101/2025.11.24.690089"
  - omics_type in {"RNASEQ", "PROTEOMICS"}
  - treatment_type contains "nitrogen"  (PRO99-lowN N-starvation vs exponential)

Excluded:
  - The two single-timepoint coculture-vs-axenic RNA-seq contrasts (treatment_type
    == ["coculture"], not "nitrogen") — no matching proteomics in the KG, out of
    scope per step 1.

Run from the multiomics_research repo root:
  .venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/01_select_experiments.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

from multiomics_explorer import list_experiments
from multiomics_explorer.analysis import experiments_to_dataframe, to_dataframe

PUBLICATION_DOI = "10.1101/2025.11.24.690089"
OMICS_TYPES = {"RNASEQ", "PROTEOMICS"}
TREATMENT_TYPE = "nitrogen"  # excludes the coculture-vs-axenic single-TP RNA-seq

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = OUT_DIR / "01_select_experiments.log"


def log(msg: str, *, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def attach_per_tp_phase(tp_df: pd.DataFrame, raw_results: list[dict]) -> pd.DataFrame:
    """Attach per-TP growth phase by indexing into time_point_growth_phases.

    The KG envelope per experiment carries `time_point_growth_phases` as a list
    parallel to `timepoints`, ordered by `timepoint_order`. experiments_to_dataframe
    drops this nesting; we re-attach by (experiment_id, timepoint_order) lookup.
    """
    phase_lookup: dict[tuple[str, int], str] = {}
    for exp in raw_results:
        eid = exp["experiment_id"]
        tps = exp.get("timepoints") or []
        phases = exp.get("time_point_growth_phases") or []
        if not tps or not phases:
            continue
        for tp, ph in zip(tps, phases):
            phase_lookup[(eid, tp["timepoint_order"])] = ph

    if "timepoint_order" in tp_df.columns:
        tp_df = tp_df.copy()
        tp_df["tp_growth_phase"] = [
            phase_lookup.get((eid, order))
            for eid, order in zip(tp_df["experiment_id"], tp_df["timepoint_order"])
        ]
    else:
        tp_df["tp_growth_phase"] = None
    return tp_df


def main() -> None:
    with LOG_PATH.open("w") as fh:
        log(f"Filter funnel for {PUBLICATION_DOI}:", fh=fh)

        result = list_experiments(
            publication_doi=[PUBLICATION_DOI], verbose=True, limit=None
        )
        n_pub = result["total_matching"]
        log(f"  publication_doi filter ........................ {n_pub} experiments", fh=fh)
        if result["truncated"]:
            sys.exit("ERROR: result truncated; rerun with limit=None")

        raw = result["results"]

        # Filter to RNA-seq + proteomics with N-stress treatment (this excludes
        # the two coculture-vs-axenic single-timepoint RNA-seq experiments).
        raw_paired = [
            r for r in raw
            if r["omics_type"] in OMICS_TYPES
            and TREATMENT_TYPE in [t.lower() for t in r.get("treatment_type", [])]
        ]
        n_paired = len(raw_paired)
        log(f"  omics_type in {sorted(OMICS_TYPES)} & treatment_type contains '{TREATMENT_TYPE}'", fh=fh)
        log(f"  ............................................... {n_paired} experiments", fh=fh)

        # Diagnostic: which experiments were dropped?
        dropped = [r for r in raw if r not in raw_paired]
        log("", fh=fh)
        log(f"Dropped from selection ({len(dropped)} experiments):", fh=fh)
        for r in dropped:
            log(f"  {r['experiment_id']}  omics={r['omics_type']}  tt={r['treatment_type']}", fh=fh)

        # Rebuild result dict for converters
        sel = {**result, "results": raw_paired, "total_matching": n_paired, "returned": n_paired}

        # Flat experiment table (debug)
        exp_df_full = to_dataframe(sel)
        exp_df_full.to_csv(OUT_DIR / "experiments_weissberg2025_full.csv", index=False)

        # Compact experiment table
        compact_cols = [
            "experiment_id",
            "experiment_name",
            "organism_name",
            "omics_type",
            "treatment_type",
            "background_factors",
            "is_time_course",
            "growth_phases",
            "time_point_growth_phases",
            "treatment",
            "control",
            "distinct_gene_count",
            "table_scope",
            "experimental_context",
        ]
        exp_df = exp_df_full[compact_cols].copy()
        exp_df.to_csv(OUT_DIR / "experiments_weissberg2025.csv", index=False)

        # Per-TP table, enriched with per-TP growth phase
        tp_df = experiments_to_dataframe(sel)
        tp_df = attach_per_tp_phase(tp_df, raw_paired)
        ordered = [
            "experiment_id",
            "organism_name",
            "omics_type",
            "background_factors",
            "is_time_course",
            "timepoint",
            "timepoint_order",
            "timepoint_hours",
            "tp_growth_phase",
            "tp_gene_count",
            "tp_significant_up",
            "tp_significant_down",
            "tp_not_significant",
        ]
        ordered = [c for c in ordered if c in tp_df.columns] + [
            c for c in tp_df.columns if c not in ordered
        ]
        tp_df = tp_df[ordered]
        tp_df.to_csv(OUT_DIR / "timepoints_weissberg2025.csv", index=False)

        # Diagnostics
        log("", fh=fh)
        log("Selected experiments (compact):", fh=fh)
        for _, row in exp_df.iterrows():
            log(
                f"  {row['experiment_id']}\n"
                f"    organism={row['organism_name']:<32} omics={row['omics_type']:<10} "
                f"bg={row['background_factors']:<25} time_course={row['is_time_course']!s:<5} "
                f"distinct_genes={row['distinct_gene_count']}",
                fh=fh,
            )

        log("", fh=fh)
        log("Per-experiment timepoint structure (with growth phase):", fh=fh)
        for exp_id, group in tp_df.groupby("experiment_id"):
            n_tp = len(group)
            log(f"  {exp_id}: {n_tp} TP row(s)", fh=fh)
            for _, r in group.iterrows():
                log(
                    f"    tp={r['timepoint']!s:<15} "
                    f"phase={r['tp_growth_phase']!s:<18} "
                    f"hours={r['timepoint_hours']!s:<6} "
                    f"tp_gene_count={r['tp_gene_count']}",
                    fh=fh,
                )

        # QC: every selected experiment should have table_scope == all_detected_genes
        bad = exp_df[exp_df["table_scope"] != "all_detected_genes"]
        log("", fh=fh)
        if not bad.empty:
            log("WARNING: selected experiments NOT using all_detected_genes table_scope:", fh=fh)
            for _, r in bad.iterrows():
                log(f"  {r['experiment_id']}: table_scope={r['table_scope']}", fh=fh)
        else:
            log("QC: all 8 selected experiments use table_scope == 'all_detected_genes'", fh=fh)

        log("", fh=fh)
        log(f"Wrote {OUT_DIR/'experiments_weissberg2025.csv'}", fh=fh)
        log(f"Wrote {OUT_DIR/'experiments_weissberg2025_full.csv'}", fh=fh)
        log(f"Wrote {OUT_DIR/'timepoints_weissberg2025.csv'}", fh=fh)
        log(f"Wrote {LOG_PATH}", fh=fh)


if __name__ == "__main__":
    main()
