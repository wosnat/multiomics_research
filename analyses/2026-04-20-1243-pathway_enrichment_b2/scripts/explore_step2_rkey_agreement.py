"""Step 2 explore: R clusters × key-pathway agreement with expected direction.

Produces a single CSV artifact + three derived summary CSVs so all agreement
claims in chat / notebook are sourced from files, not from memory.

Outputs:
- exploration/qc/step2_rkey_matrix.csv
    Full matrix: one row per (R cluster, key_pathway). Columns: class,
    first_author, experiment_short, timepoint, timepoint_hours, direction,
    term_id, term_name, expected_direction, signed_score, p_adjust, count,
    bg_count, fold_enrichment, is_significant, agrees_with_expected.

- exploration/qc/step2_rkey_summary_by_term.csv
    Per-term: n_significant, n_agree_sig, mean_abs_signed_score among sig,
    min_padj across all R clusters.

- exploration/qc/step2_rkey_summary_by_expTp.csv
    Per-(experiment, timepoint, direction): count of key pathways that are
    significant AND agree with expected direction. Lets the researcher see
    "at Tolonen 12h|down, how many expected-DOWN key pathways are enriched?"

- exploration/qc/step2_rkey_disagreements.csv
    All significant rows where direction disagrees with expected. Short list
    — should be near-zero.

- exploration/qc/step2_rkey_nc_enrichment.csv
    NC clusters × key pathways (same structure). Should be low p_adjust hits
    only, if any — lets the researcher see the noise floor.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
QC_DIR = ANALYSIS_DIR / "exploration" / "qc"


def _agreement_flag(direction: str, expected: str, is_sig: bool) -> bool | None:
    """True/False if row is significant; None if not significant (not counted)."""
    if not is_sig:
        return None
    if expected == "ambiguous":
        return True  # either direction counts; ambiguous pathway always "agrees"
    return direction == expected


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    log = logging.getLogger("explore_step2")

    enrichment = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    key_paths = pd.read_csv(ANALYSIS_DIR / "exploration" / "key_pathways.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    classified_unique = classified.drop_duplicates("experiment_id")
    class_map = classified_unique.set_index("experiment_id")["class"].to_dict()
    author_map = classified_unique.set_index("experiment_id")["first_author"].to_dict()

    # Join class + author + expected direction.
    enrichment["class"] = enrichment["experiment_id"].map(class_map)
    enrichment["first_author"] = enrichment["experiment_id"].map(author_map)

    # Filter to R and key pathways.
    key_term_ids = set(key_paths["term_id"])
    exp_dir_map = dict(zip(key_paths["term_id"], key_paths["expected_direction"]))

    def short_exp(s):
        # Pull publication-id fragment: e.g. msb4100087, ismej.2017.88, spectrum.03275-22
        for frag in ["msb4100087", "ismej.2017.88", "ismej.2016.70",
                     "spectrum.03275-22", "JB.01097-06", "mSystems.00008-17",
                     "2025.11.24.690089"]:
            if frag in s:
                return frag
        return s

    def build_matrix(class_filter: list[str]) -> pd.DataFrame:
        sub = enrichment[
            enrichment["class"].isin(class_filter)
            & enrichment["term_id"].isin(key_term_ids)
        ].copy()
        sub["expected_direction"] = sub["term_id"].map(exp_dir_map)
        sub["is_significant"] = sub["p_adjust"] < 0.05
        sub["agrees_with_expected"] = sub.apply(
            lambda r: _agreement_flag(r["direction"], r["expected_direction"],
                                      r["is_significant"]),
            axis=1,
        )
        sub["experiment_short"] = sub["experiment_id"].apply(short_exp)
        cols = [
            "class", "first_author", "experiment_short", "experiment_id",
            "timepoint", "timepoint_hours", "direction",
            "ontology", "term_id", "term_name",
            "expected_direction", "signed_score", "p_adjust",
            "count", "bg_count", "fold_enrichment",
            "is_significant", "agrees_with_expected",
        ]
        return sub[cols].sort_values(
            ["class", "first_author", "experiment_short", "timepoint_hours",
             "timepoint", "direction", "ontology", "term_id"]
        )

    # R matrix.
    r_matrix = build_matrix(["R"])
    out = QC_DIR / "step2_rkey_matrix.csv"
    r_matrix.to_csv(out, index=False)
    log.info(f"Wrote R matrix: {len(r_matrix)} rows → {out}")

    # NC matrix (noise floor).
    nc_matrix = build_matrix(["NC"])
    out = QC_DIR / "step2_rkey_nc_enrichment.csv"
    nc_matrix.to_csv(out, index=False)
    log.info(f"Wrote NC matrix: {len(nc_matrix)} rows → {out}")

    # Summary by term.
    sig_r = r_matrix[r_matrix["is_significant"]].copy()
    n_total_by_term = r_matrix.groupby("term_id").size().rename("n_tests_total")
    by_term = (
        sig_r.groupby(["term_id", "term_name", "expected_direction"])
        .agg(
            n_significant=("p_adjust", "count"),
            n_agree=("agrees_with_expected", "sum"),
            mean_abs_signed_score=("signed_score", lambda s: s.abs().mean()),
            min_padj=("p_adjust", "min"),
        )
        .reset_index()
        .join(n_total_by_term, on="term_id")
    )
    out = QC_DIR / "step2_rkey_summary_by_term.csv"
    by_term.to_csv(out, index=False)
    log.info(f"Wrote per-term summary → {out}")

    # Summary by (experiment, timepoint, direction) — count of expected-direction
    # significant hits out of 11 key pathways.
    by_expTp = (
        r_matrix.groupby(
            ["class", "experiment_short", "timepoint", "timepoint_hours", "direction"]
        )
        .apply(
            lambda g: pd.Series({
                "n_keypath_tested": len(g),
                "n_significant": int(g["is_significant"].sum()),
                "n_agree_expected": int((g["is_significant"] & (g["agrees_with_expected"] == True)).sum()),
                "n_disagree_expected": int((g["is_significant"] & (g["agrees_with_expected"] == False)).sum()),
            }),
            include_groups=False,
        )
        .reset_index()
    )
    out = QC_DIR / "step2_rkey_summary_by_expTp.csv"
    by_expTp.to_csv(out, index=False)
    log.info(f"Wrote per-experiment×TP summary → {out}")

    # Disagreements only.
    disag = r_matrix[
        r_matrix["is_significant"]
        & (r_matrix["agrees_with_expected"] == False)
    ].copy()
    out = QC_DIR / "step2_rkey_disagreements.csv"
    disag.to_csv(out, index=False)
    log.info(f"Wrote disagreements ({len(disag)} rows) → {out}")

    # Stdout summary.
    print("\n=== R-cluster × key-pathway agreement — per-term summary ===")
    print(by_term[[
        "term_id", "term_name", "expected_direction",
        "n_tests_total", "n_significant", "n_agree",
        "mean_abs_signed_score", "min_padj",
    ]].to_string(index=False))

    print("\n=== R-cluster × key-pathway agreement — per (exp, TP, dir) ===")
    print(by_expTp[[
        "class", "experiment_short", "timepoint", "direction",
        "n_keypath_tested", "n_significant", "n_agree_expected", "n_disagree_expected",
    ]].to_string(index=False))

    print("\n=== Disagreements (significant but direction mismatch with expected) ===")
    if disag.empty:
        print("  (none)")
    else:
        print(disag[[
            "experiment_short", "timepoint", "direction",
            "term_id", "term_name", "expected_direction",
            "signed_score", "p_adjust",
        ]].to_string(index=False))

    # NC noise floor.
    nc_sig = nc_matrix[nc_matrix["is_significant"]]
    print(f"\n=== NC noise floor (significant hits on key pathways) ===")
    if nc_sig.empty:
        print("  (none — clean noise floor)")
    else:
        print(nc_sig[[
            "first_author", "experiment_short", "timepoint", "direction",
            "term_id", "term_name", "signed_score", "p_adjust",
        ]].to_string(index=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
