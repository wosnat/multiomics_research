"""Select Prochlorococcus MED4 experiments from Weissberg 2025.

Inputs:  none (queries the multiomics KG via the Python API)
Outputs:
  - data/experiments_pro_med4.csv          one row per selected experiment (compact)
  - data/timepoints_pro_med4.csv           one row per (experiment, timepoint), enriched with per-TP growth phase
  - data/experiments_pro_med4_full.csv     all metadata, flat (debugging / reference)

Filter funnel (printed to stdout and to data/01_select_experiments.log).

Note on filtering: the Python API's `organism` parameter does a partial match
on BOTH profiled organism AND coculture partner. Calling
`list_experiments(organism="Prochlorococcus MED4", ...)` therefore returns 6
experiments — 5 with MED4 as the profiled organism plus 1 HOT1A3 experiment
where MED4 is the coculture partner. We filter post-hoc on
`organism_name == "Prochlorococcus MED4"` to get the 5 we want.

Run from the multiomics_research repo root:
  .venv/bin/python analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/2_kg_selection/scripts/01_select_experiments.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

from multiomics_explorer import list_experiments
from multiomics_explorer.analysis import experiments_to_dataframe, to_dataframe

PUBLICATION_DOI = "10.1101/2025.11.24.690089"
ORGANISM_NAME = "Prochlorococcus MED4"

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

        # Pull all experiments for the publication; do not pass organism filter
        # (it does partial match on profiled organism OR coculture partner).
        result = list_experiments(
            publication_doi=[PUBLICATION_DOI], verbose=True, limit=None
        )
        n_pub = result["total_matching"]
        log(f"  publication_doi filter ........................ {n_pub} experiments", fh=fh)
        if result["truncated"]:
            sys.exit("ERROR: result truncated; rerun with limit=None")

        raw = result["results"]
        # Post-hoc filter on profiled organism only
        raw_med4 = [r for r in raw if r["organism_name"] == ORGANISM_NAME]
        n_med4 = len(raw_med4)
        log(f"  organism_name == {ORGANISM_NAME!r} (post-hoc) .. {n_med4} experiments", fh=fh)

        # Rebuild a result dict with the filtered subset for the converters
        med4 = {**result, "results": raw_med4, "total_matching": n_med4, "returned": n_med4}

        # Flat experiment table (debug)
        exp_df_full = to_dataframe(med4)
        exp_df_full.to_csv(OUT_DIR / "experiments_pro_med4_full.csv", index=False)

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
        exp_df.to_csv(OUT_DIR / "experiments_pro_med4.csv", index=False)

        # Per-TP table, enriched with per-TP growth phase
        tp_df = experiments_to_dataframe(med4)
        tp_df = attach_per_tp_phase(tp_df, raw_med4)
        # Reorder for legibility
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
        tp_df.to_csv(OUT_DIR / "timepoints_pro_med4.csv", index=False)

        # Diagnostics
        log("", fh=fh)
        log("Selected experiments (compact):", fh=fh)
        for _, row in exp_df.iterrows():
            log(
                f"  {row['experiment_id']}\n"
                f"    omics={row['omics_type']:<10} "
                f"bg={row['background_factors']:<25} "
                f"time_course={row['is_time_course']!s:<5} "
                f"genes={row['distinct_gene_count']}",
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

        # Cross-condition TP overlap, calendar-shared AND physiologically-shared
        log("", fh=fh)
        log("Cross-condition TP overlap (axenic vs coculture, per omics):", fh=fh)
        tp_df["_axenic"] = tp_df["background_factors"].astype(str).str.contains("axenic")
        tp_df["_coculture"] = tp_df["background_factors"].astype(str).str.contains("coculture")
        for omics in ["RNASEQ", "PROTEOMICS"]:
            sub = tp_df[tp_df["omics_type"] == omics]
            ax = sub[sub["_axenic"]][["timepoint", "tp_growth_phase"]].dropna(subset=["timepoint"])
            co = sub[sub["_coculture"]][["timepoint", "tp_growth_phase"]].dropna(subset=["timepoint"])
            ax_set = set(zip(ax["timepoint"], ax["tp_growth_phase"]))
            co_set = set(zip(co["timepoint"], co["tp_growth_phase"]))
            calendar_shared = sorted({t for (t, _) in ax_set} & {t for (t, _) in co_set})
            phys_shared = sorted(ax_set & co_set, key=lambda x: x[0])
            log(f"  [{omics}] axenic (TP, phase):    {sorted(ax_set, key=lambda x: x[0])}", fh=fh)
            log(f"  [{omics}] coculture (TP, phase): {sorted(co_set, key=lambda x: x[0])}", fh=fh)
            log(f"  [{omics}] calendar-shared TPs:    {calendar_shared}", fh=fh)
            log(f"  [{omics}] phys-shared (TP+phase): {phys_shared}", fh=fh)

        log("", fh=fh)
        log(f"Wrote {OUT_DIR/'experiments_pro_med4.csv'}", fh=fh)
        log(f"Wrote {OUT_DIR/'experiments_pro_med4_full.csv'}", fh=fh)
        log(f"Wrote {OUT_DIR/'timepoints_pro_med4.csv'}", fh=fh)
        log(f"Wrote {LOG_PATH}", fh=fh)


if __name__ == "__main__":
    main()
