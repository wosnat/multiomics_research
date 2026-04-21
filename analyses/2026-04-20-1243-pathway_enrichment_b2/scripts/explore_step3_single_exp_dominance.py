"""Single-R-experiment dominance check on the Step 3 signature.

For each signature pathway, expand per_experiment_breakdown (JSON) into a
long-form table of (term × experiment × direction → cluster count). Flag
pathways whose same-direction support comes from a SINGLE R experiment
(the other R experiment contributes zero in-direction clusters). These
are candidates for M4 stability weakness in Step 4 LOO-R.

Outputs:
    exploration/qc/step3_per_exp_support.csv       (long form)
    exploration/qc/step3_single_exp_dominated.csv  (flagged subset)
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


def main() -> None:
    sig = pd.read_csv(ANALYSIS_DIR / "data" / "reference_signature.csv")

    rows = []
    for _, row in sig.iterrows():
        direction = row["direction"]
        breakdown = json.loads(row["per_experiment_breakdown"])
        for exp_id, counts in breakdown.items():
            rows.append({
                "ontology": row["ontology"],
                "term_id": row["term_id"],
                "term_name": row["term_name"],
                "direction": direction,
                "experiment_id": exp_id,
                "n_up": counts.get("up", 0),
                "n_down": counts.get("down", 0),
                "n_in_direction": counts.get(direction, 0),
            })
    per_exp = pd.DataFrame(rows)
    per_exp.to_csv(
        ANALYSIS_DIR / "exploration" / "qc" / "step3_per_exp_support.csv",
        index=False,
    )

    # For each term, count how many R experiments contribute >=1 in-direction
    # cluster. Flag terms with exactly 1. Also compute share concentration.
    agg = (
        per_exp.groupby(["ontology", "term_id", "term_name", "direction"])
        .agg(
            n_exps_contributing=("n_in_direction", lambda s: (s > 0).sum()),
            n_exps_present=("experiment_id", "nunique"),
            total_in_direction=("n_in_direction", "sum"),
            max_per_exp=("n_in_direction", "max"),
            contributing_exps=(
                "experiment_id",
                lambda s: "|".join(
                    sorted(
                        per_exp.loc[s.index].loc[
                            per_exp.loc[s.index, "n_in_direction"] > 0,
                            "experiment_id",
                        ]
                    )
                ),
            ),
        )
        .reset_index()
    )
    agg["single_exp_dominated"] = agg["n_exps_contributing"] == 1
    agg["share_max_exp"] = (
        agg["max_per_exp"] / agg["total_in_direction"].clip(lower=1)
    )

    flagged = agg[agg["single_exp_dominated"]].copy()
    flagged.to_csv(
        ANALYSIS_DIR / "exploration" / "qc" / "step3_single_exp_dominated.csv",
        index=False,
    )

    print("\n=== Single-R-experiment dominance ===")
    print(f"Total signature terms: {len(agg)}")
    print(f"Dominated by 1 R experiment: {len(flagged)}")
    if len(flagged):
        print(flagged[[
            "ontology", "term_id", "term_name", "direction",
            "contributing_exps", "total_in_direction",
        ]].to_string(index=False))

    print("\n=== Per-experiment support summary ===")
    print(agg[[
        "ontology", "term_id", "direction",
        "n_exps_contributing", "total_in_direction",
        "max_per_exp", "share_max_exp",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
