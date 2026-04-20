"""Main enriched pathways in R clusters with timepoint_hours > 3.

Writes exploration/qc/step2_r_top_pathways.csv — per (ontology, term, direction):
n_clusters (padj<0.05), mean/min p_adjust, mean |signed_score|,
contributing experiments.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


def main() -> None:
    df = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")

    r_exps = set(classified.loc[classified["class"] == "R", "experiment_id"])
    r = df[df["experiment_id"].isin(r_exps)].copy()
    # Apply resolved decision #1: exclude timepoint_hours <= 3.
    r = r[r["timepoint_hours"] > 3]
    sig = r[r["p_adjust"] < 0.05].copy()

    print(f"R clusters after filter: {sig['cluster'].nunique()} "
          f"({sig['experiment_id'].nunique()} experiments)")
    print(f"Total sig rows: {len(sig)}, unique terms: {sig['term_id'].nunique()}")

    sig["exp_short"] = sig["experiment_id"].str.split("_").str[0]
    agg = (
        sig.groupby(["ontology", "term_id", "term_name", "direction"])
           .agg(n_clusters=("cluster", "nunique"),
                min_padj=("p_adjust", "min"),
                mean_abs_signed=("signed_score", lambda s: s.abs().mean()),
                exps=("exp_short", lambda s: "|".join(sorted(set(s)))))
           .reset_index()
           .sort_values(["n_clusters", "mean_abs_signed"], ascending=[False, False])
    )

    out = ANALYSIS_DIR / "exploration" / "qc" / "step2_r_top_pathways.csv"
    agg.to_csv(out, index=False)
    print(f"\nwrote {out} — {len(agg)} (term × direction) rows\n")

    # Per ontology, show top 15 by n_clusters then |signed|.
    for ont in agg["ontology"].unique():
        sub = agg[agg["ontology"] == ont].head(20)
        print(f"\n=== {ont}: top {len(sub)} by n_clusters ===")
        cols = ["term_id", "term_name", "direction", "n_clusters",
                "mean_abs_signed", "min_padj", "exps"]
        with pd.option_context("display.max_colwidth", 60,
                                "display.width", 200):
            print(sub[cols].to_string(index=False))


if __name__ == "__main__":
    main()
