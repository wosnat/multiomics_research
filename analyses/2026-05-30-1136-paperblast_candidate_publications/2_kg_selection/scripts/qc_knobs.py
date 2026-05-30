#!/usr/bin/env python
"""
qc_knobs.py — compare per-experiment "top responder" set size under different knobs.

For ONE experiment: pull significant DE rows once, then show how the
deduplicated (per locus_tag) top-responder count changes under:
  - top-N threshold (5 vs 10)
  - direction (up only / down only / both)
  - min number of timepoints a gene must be top-N in (1 vs 2)
Also prints the per-row rank histogram.

Run from repo root:
  env -u VIRTUAL_ENV uv run python analyses/2026-05-30-1136-paperblast_candidate_publications/\
2_kg_selection/scripts/qc_knobs.py \
    --experiment 10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray
"""
from __future__ import annotations

import argparse

import pandas as pd

from multiomics_explorer import differential_expression_by_gene, to_dataframe


def load(experiment_id: str) -> pd.DataFrame:
    res = differential_expression_by_gene(
        experiment_ids=[experiment_id], significant_only=True,
        direction="both", verbose=True,
    )
    df = to_dataframe(res)
    df["row_rank"] = df[["rank_up", "rank_down"]].min(axis=1, skipna=True)
    df = df.dropna(subset=["row_rank"]).copy()
    df["dir"] = df.apply(lambda r: "up" if r["rank_up"] == r["row_rank"] else "down", axis=1)
    return df


def count(df: pd.DataFrame, top_n: int, direction: str, min_tps: int) -> int:
    d = df if direction == "both" else df[df["dir"] == direction]
    hits = d[d["row_rank"] <= top_n]
    if min_tps <= 1:
        return hits["locus_tag"].nunique()
    # genes that are top-N in >= min_tps distinct timepoints
    per = hits.groupby("locus_tag")["timepoint"].nunique()
    return int((per >= min_tps).sum())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", required=True)
    args = ap.parse_args()
    df = load(args.experiment)

    n_tps = df["timepoint"].nunique()
    print(f"# experiment: {args.experiment}")
    print(f"# significant rows: {len(df)} across {n_tps} timepoint(s), "
          f"{df['locus_tag'].nunique()} distinct genes\n")

    print("per-row best-rank histogram (rank bucket -> # rows):")
    buckets = pd.cut(df["row_rank"], bins=[0, 5, 10, 20, 50, 1e9],
                     labels=["1-5", "6-10", "11-20", "21-50", "50+"])
    print(buckets.value_counts().sort_index().to_string(), "\n")

    print(f"{'knob':38} {'count':>6}")
    for direction in ("both", "up", "down"):
        for top_n in (5, 10):
            for min_tps in (1, 2):
                c = count(df, top_n, direction, min_tps)
                label = f"dir={direction:4} top-{top_n} min_tps={min_tps}"
                print(f"{label:38} {c:>6}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
