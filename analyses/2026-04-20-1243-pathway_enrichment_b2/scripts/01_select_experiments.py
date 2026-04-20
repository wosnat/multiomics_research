"""Select and classify experiments for pathway_enrichment_b2 analysis.

Hard-codes the T/R/PC/NC/CTX classifications decided interactively in Step 1a.
Re-running this script reproduces experiments_classified.csv from the current KG.

Attribution note (per gaps_and_friction.md, Skill 2026-04-20 entry):
publication_doi → first_author + year is joined in via list_publications so the
paper attribution travels with the classified CSV. Do not name authors from
intrinsic memory.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd
from multiomics_explorer import list_experiments, list_publications
from multiomics_explorer.analysis import experiments_to_dataframe

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

# FILLED IN FROM STEP 1a INTERACTIVE CLASSIFICATION (researcher-approved 2026-04-20).
# class values: "T" (target), "R" (reference), "PC", "NC", "CTX".
CLASSIFICATIONS: list[dict[str, str]] = [
    # --- T (4): Weissberg 2025 MED4 N-starvation target conditions, all_detected_genes ---
    {"experiment_id": "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic",
     "class": "T",
     "rationale": "Weissberg 2025 MED4 axenic RNA-seq (single-timepoint N-limited snapshot)"},
    {"experiment_id": "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_coculture",
     "class": "T",
     "rationale": "Weissberg 2025 MED4 coculture RNA-seq time course (day 18/31/60/89)"},
    {"experiment_id": "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic",
     "class": "T",
     "rationale": "Weissberg 2025 MED4 axenic proteomics time course (day 14/31/89)"},
    {"experiment_id": "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_coculture",
     "class": "T",
     "rationale": "Weissberg 2025 MED4 coculture proteomics time course (day 18/31/60/89)"},

    # --- R (2): MED4 unambiguous N-limitation references ---
    {"experiment_id": "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray",
     "class": "R",
     "rationale": "Tolonen 2006 MED4 N-deprivation microarray, 6 timepoints (0/3/6/12/24/48h); all_detected_genes; canonical reference"},
    {"experiment_id": "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq",
     "class": "R",
     "rationale": "Read 2017 MED4 N-depleted RNA-seq, 3 timepoints (3/12/24h); filtered_subset (top 50% by expression) -> organism bg"},

    # --- PC (2): MED4 mild N-stress (alt N sources, growable) ---
    {"experiment_id": "10.1038/msb4100087_growth_medium_growth_on_cyanate_as_med4_microarray",
     "class": "PC",
     "rationale": "Tolonen 2006 MED4 on cyanate as sole N; exponential phase, can grow but under N-stress relative to NH4"},
    {"experiment_id": "10.1038/msb4100087_growth_medium_growth_on_urea_as_med4_microarray",
     "class": "PC",
     "rationale": "Tolonen 2006 MED4 on urea as sole N; exponential phase, mild N-stress"},

    # --- NC (5): MED4 experiments unrelated to N ---
    {"experiment_id": "10.1038/ismej.2016.70_coculture_alteromonas_hot1a3_med4_rnaseq",
     "class": "NC",
     "rationale": "Aharonovich 2016 MED4 coculture with HOT1A3 vs axenic; nutrient-replete RNA-seq, all_detected_genes"},
    {"experiment_id": "10.1101/2025.11.24.690089_coculture_alteromonas_hot1a3_med4_rnaseq",
     "class": "NC",
     "rationale": "Weissberg 2025 MED4 coculture vs axenic (nutrient-replete, exponential); matched-design RNA-seq NC for the coculture T"},
    {"experiment_id": "10.1128/spectrum.03275-22_light_high_glucose_med4_proteomics",
     "class": "NC",
     "rationale": "Moreno-Cabezuelo 2023 MED4 light+high glucose proteomics; all_detected_genes NC"},
    {"experiment_id": "10.1128/spectrum.03275-22_dark_high_glucose_med4_proteomics",
     "class": "NC",
     "rationale": "Moreno-Cabezuelo 2023 MED4 dark+high glucose proteomics; all_detected_genes NC"},
    {"experiment_id": "10.1128/JB.01097-06_light_stress_high_white_light_55_med4_microarray",
     "class": "NC",
     "rationale": "Steglich 2006 MED4 high white light (55 umol/m2/s); filtered_subset -> organism bg (provides organism-bg NC calibration)"},

    # --- CTX (2): non-MED4 N-limitation for figure context ---
    # Dropped 2026-04-20 during Step 1b ontology review:
    #   HOT1A3 Weissberg axenic N-starvation RNA-seq — Alteromonas has no
    #   cyanorak_role or tigr_role annotations (Gammaproteobacteria), so with
    #   the MED4-optimal cyanorak_role L1 pick it would appear as a blank
    #   column. Sister-Prochlorococcus conservation test (MIT9313 + SS120)
    #   is the right scope for the MED4 N-signature.
    {"experiment_id": "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray",
     "class": "CTX",
     "rationale": "Tolonen 2006 Prochlorococcus MIT9313 N-deprivation microarray; sister Prochlorococcus strain, same experimental design as R1"},
    {"experiment_id": "10.1128/mSystems.00008-17_azaserine_ss120_proteomics",
     "class": "CTX",
     "rationale": "SS120 azaserine (GS inhibitor, N-limitation proxy) proteomics; sister Prochlorococcus, significant_only -> organism bg"},
]


def _publications_dataframe() -> pd.DataFrame:
    """Fetch list_publications and reduce to doi + first_author + year + title."""
    pub_result = list_publications(verbose=True, limit=None)
    pub_rows = pub_result.get("results", []) if isinstance(pub_result, dict) else []
    if not pub_rows:
        return pd.DataFrame(columns=["publication_doi", "first_author", "publication_year", "publication_title_api"])
    rows = []
    for p in pub_rows:
        authors = p.get("authors") or []
        rows.append({
            "publication_doi": p.get("doi"),
            "first_author": authors[0] if authors else None,
            "publication_year": p.get("year"),
            "publication_title_api": p.get("title"),
        })
    return pd.DataFrame(rows)


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step1a.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step1a")
    log.info("Starting 01_select_experiments.py")

    if not CLASSIFICATIONS:
        log.error("CLASSIFICATIONS is empty — fill in from Step 1a.")
        return 1

    classified_df = pd.DataFrame(CLASSIFICATIONS)
    class_counts = classified_df["class"].value_counts().to_dict()
    log.info(f"Loaded {len(classified_df)} classifications: {class_counts}")

    # API gap workaround: list_experiments has no experiment_ids filter (see api_coverage.md A1).
    # Pull the full landscape and filter locally.
    result = list_experiments(verbose=True, limit=None)
    combined = experiments_to_dataframe(result)
    log.info(f"Pulled {len(combined)} experiment × timepoint rows from list_experiments(limit=None)")
    if combined.empty:
        log.error("list_experiments returned empty landscape; cannot proceed.")
        return 1

    # Filter to selected experiment_ids; join classification.
    selected = combined[combined["experiment_id"].isin(classified_df["experiment_id"])].copy()
    selected = selected.merge(classified_df, on="experiment_id", how="left")
    if selected.empty:
        log.error("No matching experiments found in the KG after filtering; check experiment_id spellings.")
        return 2

    # Fail loudly on any classified experiment_id that didn't land in the join.
    missing = sorted(set(classified_df["experiment_id"]) - set(selected["experiment_id"]))
    if missing:
        log.error(f"Missing experiments after metadata join ({len(missing)}): {missing}")
        print(
            f"ERROR: {len(missing)} classified experiment_ids not found in list_experiments: "
            f"{missing}",
            file=sys.stderr,
        )
        return 3

    # Publication attribution join (see gaps_and_friction.md, Skill 2026-04-20).
    pub_df = _publications_dataframe()
    log.info(f"Pulled {len(pub_df)} publication rows for attribution join")
    selected = selected.merge(pub_df, on="publication_doi", how="left")
    missing_pubs_df = selected.loc[selected["first_author"].isna(), ["experiment_id", "publication_doi"]].drop_duplicates()
    if len(missing_pubs_df):
        log.warning(f"Attribution missing for {len(missing_pubs_df)} experiment_id × doi pairs: "
                    f"{missing_pubs_df.to_dict(orient='records')}")

    # Persist.
    out_path = ANALYSIS_DIR / "data" / "experiments_classified.csv"
    selected.to_csv(out_path, index=False)
    log.info(f"Wrote {len(selected)} rows ({selected['experiment_id'].nunique()} unique experiments) to {out_path}")

    # Summaries to stdout for chat inspection.
    unique = selected.drop_duplicates("experiment_id")
    print("=== Class × organism × omics (unique experiments) ===")
    print(unique.groupby(["class", "organism_name", "omics_type"]).size().to_string())
    print("\n=== table_scope × class (unique experiments) ===")
    print(unique.groupby(["class", "table_scope"]).size().to_string())
    print("\n=== gene_count range per class (unique experiments) ===")
    print(unique.groupby("class")["gene_count"].agg(["min", "median", "max"]).to_string())
    print("\n=== Per-experiment attribution (unique experiments) ===")
    cols = ["class", "experiment_id", "first_author", "publication_year", "organism_name", "omics_type", "table_scope", "is_time_course"]
    print(unique[cols].to_string(index=False))

    # Cluster-count expectation (informational; for Step 3 signature-support arithmetic).
    # Rows-per-experiment in experiments_to_dataframe = # of timepoints (or 1 if non-time-course).
    rows_per_exp = selected.groupby("experiment_id").size()
    r_ids = set(classified_df.loc[classified_df["class"] == "R", "experiment_id"])
    r_rows = rows_per_exp.reindex(sorted(r_ids))
    log.info(f"R experiment timepoint counts: {r_rows.to_dict()}")
    print("\n=== R experiment timepoint counts (informs Step 3 support arithmetic) ===")
    print(r_rows.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
