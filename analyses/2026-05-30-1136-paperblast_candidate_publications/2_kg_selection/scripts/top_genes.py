#!/usr/bin/env python
"""
top_genes.py — identify the top responder genes of ONE experiment.

A gene is a "top responder" if its best (lowest) rank among significant rows is
<= TOP_N. Both directions: each significant row carries rank_up (sig-up) or
rank_down (sig-down), 1 = strongest |log2FC| within experiment x timepoint.

Conditional timepoint rule (locked 2026-05-30): for a time-course experiment
(>=2 distinct timepoints) a gene must be top-N in at least MIN_TPS (default 2)
timepoints to count; for a single-contrast experiment (1 timepoint) one hit
suffices. This drops transient single-timepoint spikes in time-courses without
zeroing out single-contrast experiments.

We collapse a gene's timepoints to its single best rank via pandas.

Run from repo root, e.g.:
  env -u VIRTUAL_ENV uv run python analyses/2026-05-30-1136-paperblast_candidate_publications/\
2_kg_selection/scripts/top_genes.py \
    --experiment 10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray --top-n 10
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from multiomics_explorer import differential_expression_by_gene, to_dataframe

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
pd.set_option("display.max_rows", None)
pd.set_option("display.width", 200)
pd.set_option("display.max_colwidth", 50)


def top_genes(experiment_id: str, top_n: int, min_tps: int = 2, conn=None) -> pd.DataFrame:
    res = differential_expression_by_gene(
        experiment_ids=[experiment_id],
        significant_only=True, direction="both", verbose=True, conn=conn,
    )
    organism = res.get("organism_name")
    df = to_dataframe(res)
    if df.empty:
        return df
    df["organism"] = organism

    # best rank across up/down per row
    df["row_rank"] = df[["rank_up", "rank_down"]].min(axis=1, skipna=True)
    df = df.dropna(subset=["row_rank"])
    df["direction"] = df.apply(
        lambda r: "up" if r["rank_up"] == r["row_rank"] else "down", axis=1)

    # conditional timepoint rule: time-course (>=2 tps) requires top-N hits in
    # >= min_tps distinct timepoints; single-contrast (1 tp) requires just 1.
    n_tps = df["timepoint"].nunique()
    hits = df[df["row_rank"] <= top_n]
    if n_tps >= 2 and min_tps > 1:
        tp_counts = hits.groupby("locus_tag")["timepoint"].nunique()
        keep = set(tp_counts[tp_counts >= min_tps].index)
        hits = hits[hits["locus_tag"].isin(keep)]

    # one row per locus_tag = the gene's best (lowest) rank across timepoints
    idx = hits.groupby("locus_tag")["row_rank"].idxmin()
    best = hits.loc[idx].copy()
    best = best.sort_values(["row_rank", "locus_tag"])
    best = best.drop_duplicates(subset="locus_tag", keep="first")  # explicit safety
    cols = ["organism", "locus_tag", "gene_name", "product", "gene_category",
            "row_rank", "direction", "log2fc", "padj", "timepoint"]
    return best[[c for c in cols if c in best.columns]].rename(columns={"row_rank": "best_rank"})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", required=True)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--min-tps", type=int, default=2,
                    help="for time-course experiments, min timepoints a gene must be top-N in")
    ap.add_argument("--out", action="store_true", help="also write CSV to ../data/")
    args = ap.parse_args()

    df = top_genes(args.experiment, args.top_n, args.min_tps)
    print(f"# experiment: {args.experiment}")
    print(f"# top responders (rank<={args.top_n}, both dirs, "
          f">={args.min_tps} tp if time-course): {len(df)}\n")
    print(df.to_string(index=False))

    if args.out:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        slug = args.experiment.replace("/", "_").replace(".", "_")
        path = DATA_DIR / f"top_genes_{slug}.csv"
        df.to_csv(path, index=False)
        print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
