"""Decompose T signature scores into per-term contributions + cross-ontology agreement.

Two outputs:

1. `exploration/qc/step4_t_top_contributors.csv` — for each T (exp × tp × ontology)
   row, enumerate per-term contributions (already applied: max-abs cluster pick,
   threshold, sign_ref × capped_signed_score). Ranked by |contribution|.
   Answers: "which specific signature pathways are driving this T row's score?"

2. `exploration/qc/step4_cross_ontology_agreement.csv` — per T (exp × tp), compare
   cyanorak and kegg scores. Agreement = both positive, both negative, or both
   near-zero. Disagreement = one positive + one negative in the same direction.
"""
from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path

import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ANALYSIS_DIR / "scripts"))
scoring = import_module("05_compute_scores")


def decompose_contributions(
    enrichment: pd.DataFrame, signature: pd.DataFrame, t_rows: pd.DataFrame
) -> pd.DataFrame:
    """For each T (exp × tp × ontology × bg), compute per-term max-abs
    contributions (same rule as compute_signature_score) and return long-form
    per-term rows.
    """
    rows = []
    for _, trow in t_rows.iterrows():
        exp_id = trow["experiment_id"]
        tp = trow["timepoint"]
        ont = trow["ontology"]
        bg = trow["background_used"]
        up_rows = enrichment[
            (enrichment["experiment_id"] == exp_id)
            & (enrichment["timepoint"] == tp)
            & (enrichment["ontology"] == ont)
            & (enrichment["direction"] == "up")
        ]
        down_rows = enrichment[
            (enrichment["experiment_id"] == exp_id)
            & (enrichment["timepoint"] == tp)
            & (enrichment["ontology"] == ont)
            & (enrichment["direction"] == "down")
        ]
        ont_sig = signature[signature["ontology"] == ont]
        for _, sigrow in ont_sig.iterrows():
            term_id = sigrow["term_id"]
            ref_dir = sigrow["direction"]
            sign_ref = 1.0 if ref_dir == "up" else -1.0
            up_match = up_rows[up_rows["term_id"] == term_id]
            down_match = down_rows[down_rows["term_id"] == term_id]
            up_signed = up_match.iloc[0]["signed_score"] if not up_match.empty else np.nan
            down_signed = down_match.iloc[0]["signed_score"] if not down_match.empty else np.nan
            up_c = scoring.thresholded_capped_contrib(up_signed, sign_ref)
            down_c = scoring.thresholded_capped_contrib(down_signed, sign_ref)
            max_abs_c = up_c if abs(up_c) >= abs(down_c) else down_c
            which = "-"
            if max_abs_c == up_c and up_c != 0:
                which = "UP"
            elif max_abs_c == down_c and down_c != 0:
                which = "DOWN"
            rows.append({
                "category": trow["category"],
                "omics": trow["omics"],
                "timepoint": tp,
                "ontology": ont,
                "experiment_id": exp_id,
                "term_id": term_id,
                "term_name": sigrow["term_name"],
                "expected_direction": ref_dir,
                "up_signed_score": up_signed if not pd.isna(up_signed) else None,
                "down_signed_score": down_signed if not pd.isna(down_signed) else None,
                "up_contribution": up_c,
                "down_contribution": down_c,
                "max_abs_contribution": max_abs_c,
                "contributing_cluster": which,
            })
    return pd.DataFrame(rows)


def top_contributors(decomp: pd.DataFrame, top_n: int = 3) -> pd.DataFrame:
    """Per (category × omics × tp × ontology), keep top-N terms by |max_abs_contribution|."""
    decomp = decomp.copy()
    decomp["abs_contrib"] = decomp["max_abs_contribution"].abs()
    decomp = decomp[decomp["abs_contrib"] > 0]  # drop zeros
    return (
        decomp.sort_values(
            ["category", "omics", "timepoint", "ontology", "abs_contrib"],
            ascending=[True, True, True, True, False],
        )
        .groupby(["category", "omics", "timepoint", "ontology"], sort=False)
        .head(top_n)
    )


def cross_ontology_agreement(scores_all: pd.DataFrame) -> pd.DataFrame:
    """Per T (exp × tp × bg), compare cyanorak vs kegg final_score / up / down.
    Agreement categories: both_positive, both_negative, both_zero (below 0.5),
    mixed, one_missing.
    """
    t = scores_all[scores_all["class"] == "T"].copy()
    rows = []
    for (exp_id, tp, bg), grp in t.groupby(["experiment_id", "timepoint", "background_used"]):
        cyan = grp[grp["ontology"] == "cyanorak_role"]
        kegg = grp[grp["ontology"] == "kegg"]
        if cyan.empty or kegg.empty:
            continue
        cyan = cyan.iloc[0]
        kegg = kegg.iloc[0]
        entry = {
            "category": cyan["category"],
            "omics": cyan["omics"],
            "timepoint": tp,
            "experiment_id": exp_id,
            "background_used": bg,
            "cyanorak_final": cyan["final_score"],
            "kegg_final": kegg["final_score"],
            "cyanorak_up": cyan["score_up"],
            "kegg_up": kegg["score_up"],
            "cyanorak_down": cyan["score_down"],
            "kegg_down": kegg["score_down"],
        }
        # Agreement on final_score. Categories:
        # - both_weak: both in [-tol, tol] — neither engaged
        # - same_sign: both same sign (at least one above tol)
        # - opposite_sign: one positive, one negative (actual disagreement)
        # - missing: NaN in either
        def categorize_pair(a, b, tol=0.5):
            if pd.isna(a) or pd.isna(b):
                return "missing"
            if abs(a) <= tol and abs(b) <= tol:
                return "both_weak"
            if (a > 0 and b >= 0) or (a < 0 and b <= 0) or (a >= 0 and b > 0) or (a <= 0 and b < 0):
                return "same_sign"
            return "opposite_sign"
        entry["agreement_final"] = categorize_pair(cyan["final_score"], kegg["final_score"])
        entry["agreement_up"] = categorize_pair(cyan["score_up"], kegg["score_up"])
        entry["agreement_down"] = categorize_pair(cyan["score_down"], kegg["score_down"])
        rows.append(entry)
    return pd.DataFrame(rows)


def main() -> None:
    enrichment = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    signature = pd.read_csv(ANALYSIS_DIR / "data" / "reference_signature.csv")
    scores_all = pd.read_csv(ANALYSIS_DIR / "results" / "scores_all.csv")
    t_rows = scores_all[scores_all["class"] == "T"].copy()

    # 1. Per-term contribution decomposition
    decomp = decompose_contributions(enrichment, signature, t_rows)
    decomp.to_csv(
        ANALYSIS_DIR / "exploration" / "qc" / "step4_t_contribution_decomposition.csv",
        index=False,
    )

    # Top-3 contributors per T (exp × tp × ontology)
    top3 = top_contributors(decomp, top_n=3)
    top3.to_csv(
        ANALYSIS_DIR / "exploration" / "qc" / "step4_t_top_contributors.csv",
        index=False,
    )

    # 2. Cross-ontology agreement
    agreement = cross_ontology_agreement(scores_all)
    agreement.to_csv(
        ANALYSIS_DIR / "exploration" / "qc" / "step4_cross_ontology_agreement.csv",
        index=False,
    )

    # Stdout summaries
    print("\n=== Top contributors per T (exp × tp × ontology) — sample ===")
    print(top3[[
        "category", "omics", "timepoint", "ontology",
        "term_id", "expected_direction", "max_abs_contribution", "contributing_cluster",
    ]].head(30).to_string(index=False))

    print("\n=== Cross-ontology agreement on final_score ===")
    print(agreement.groupby("agreement_final").size().to_string())

    print("\n=== Cross-ontology agreement on score_up ===")
    print(agreement.groupby("agreement_up").size().to_string())

    print("\n=== Cross-ontology agreement on score_down ===")
    print(agreement.groupby("agreement_down").size().to_string())

    print("\n=== Opposite-sign disagreements on final_score ===")
    opp = agreement[agreement["agreement_final"] == "opposite_sign"]
    if not opp.empty:
        print(opp[[
            "category", "omics", "timepoint",
            "cyanorak_final", "kegg_final",
            "cyanorak_up", "kegg_up", "cyanorak_down", "kegg_down",
        ]].to_string(index=False))
    else:
        print("(none)")


if __name__ == "__main__":
    main()
