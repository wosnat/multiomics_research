"""Per-cluster significance counts — full ontology, not just key pathways.

Produces exploration/qc/step2_cluster_totals.csv with:
- one row per (organism, ontology, bg, cluster)
- n_tested, n_sig (padj<0.05), n_sig_up, n_sig_down

Answers: "when the rkey summary says Tolonen 3h up has 2 sig, is that
2 of 11 key pathways or 2 of the whole ontology?" — 11 is key pathways
only; this script gives the whole-ontology picture.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


def main() -> None:
    df = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    df["is_sig"] = df["p_adjust"] < 0.05

    totals = (
        df.groupby(["organism", "ontology", "background_used", "experiment_id",
                    "timepoint", "direction"], dropna=False)
          .agg(n_tested=("p_adjust", "size"),
               n_sig=("is_sig", "sum"))
          .reset_index()
    )
    # Shorten experiment id to last DOI fragment for readability.
    totals["experiment_short"] = totals["experiment_id"].str.extract(
        r"(?:10\.\d+/)?([^/]+)")[0].str.split("_").str[0]

    out = ANALYSIS_DIR / "exploration" / "qc" / "step2_cluster_totals.csv"
    totals.to_csv(out, index=False)
    print(f"wrote {out} — {len(totals)} cluster-rows")

    # Focused view: R clusters at early TPs.
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    r_exps = set(classified.loc[classified["class"] == "R", "experiment_id"])
    r_view = totals[totals["experiment_id"].isin(r_exps)].copy()
    r_view = r_view.sort_values(["experiment_id", "timepoint", "ontology",
                                  "direction"])
    print("\n=== R clusters — total significant per (ontology, cluster) ===")
    cols = ["experiment_short", "timepoint", "direction", "ontology",
            "background_used", "n_tested", "n_sig"]
    print(r_view[cols].to_string(index=False))


if __name__ == "__main__":
    main()
