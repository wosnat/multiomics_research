#!/usr/bin/env python
"""
02_all_experiments.py — pool top-responder genes across ALL experiments.

For every experiment, apply the locked per-experiment cut (top_genes: rank<=TOP_N,
both directions, >=MIN_TPS timepoints when time-course). Pool the results WITHOUT
deduplicating genes across experiments — a gene top in 5 experiments contributes
5 rows. This is the "before deduplication" pool, used to inspect the per-organism
and per-genus distribution before we decide dedup / quota.

Genus comes from list_organisms (KG as source), not string-splitting.
Experiments that error inside the package's diagnostic query are skipped + logged.

Outputs (../data/):
  02_pooled_top_responders.csv   one row per (experiment x gene), before dedup
  02_distribution_by_organism.csv
  02_distribution_by_genus.csv
  02_all_experiments.log

Run from repo root:
  env -u VIRTUAL_ENV uv run python analyses/2026-05-30-1136-paperblast_candidate_publications/\
2_kg_selection/scripts/02_all_experiments.py --top-n 10 --min-tps 2
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from multiomics_explorer import GraphConnection, list_experiments, list_organisms

from top_genes import top_genes  # same-dir import

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
pd.set_option("display.width", 200)
pd.set_option("display.max_rows", None)


def genus_map(conn) -> dict[str, str]:
    orgs = list_organisms(limit=None, conn=conn)["results"]
    return {o["organism_name"]: (o.get("genus") or "?") for o in orgs}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--min-tps", type=int, default=2)
    args = ap.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logf = open(DATA_DIR / "02_all_experiments.log", "w")

    def log(m: str) -> None:
        print(m)
        logf.write(m + "\n")

    log(f"== all_experiments: top_n={args.top_n} min_tps={args.min_tps} ==")

    with GraphConnection() as conn:
        g_map = genus_map(conn)
        exps = list_experiments(limit=None, conn=conn)["results"]
        log(f"experiments: {len(exps)} total; organisms in KG: {len(g_map)}")

        frames, n_ok, n_err, n_empty = [], 0, 0, 0
        for i, e in enumerate(exps, 1):
            eid = e.get("experiment_id")
            if not eid:
                continue
            try:
                df = top_genes(eid, args.top_n, args.min_tps, conn=conn)
            except Exception as ex:
                n_err += 1
                log(f"  ERR skip {e.get('omics_type')} {eid}: {type(ex).__name__}")
                continue
            if df.empty:
                n_empty += 1
                continue
            df = df.copy()
            df["experiment_id"] = eid
            df["omics_type"] = e.get("omics_type")
            frames.append(df)
            n_ok += 1
            if i % 25 == 0:
                log(f"  ...{i}/{len(exps)} processed")

    log(f"experiments: {n_ok} contributed, {n_err} errored (package bug), {n_empty} empty")

    pooled = pd.concat(frames, ignore_index=True)
    pooled["genus"] = pooled["organism"].map(g_map).fillna("?")
    # column order
    front = ["experiment_id", "organism", "genus", "omics_type", "locus_tag",
             "gene_name", "product", "gene_category", "best_rank", "direction",
             "log2fc", "padj", "timepoint"]
    pooled = pooled[[c for c in front if c in pooled.columns]]
    pooled.to_csv(DATA_DIR / "02_pooled_top_responders.csv", index=False)

    log(f"\npooled rows (before dedup): {len(pooled)}")
    log(f"distinct genes (locus_tag): {pooled['locus_tag'].nunique()}")
    log(f"distinct (organism, locus_tag): {pooled.groupby(['organism','locus_tag']).ngroups}")

    # distribution by organism (before dedup = row count; plus distinct genes)
    by_org = pooled.groupby("organism").agg(
        pooled_rows=("locus_tag", "size"),
        distinct_genes=("locus_tag", "nunique"),
        experiments=("experiment_id", "nunique"),
    ).sort_values("pooled_rows", ascending=False)
    by_org.to_csv(DATA_DIR / "02_distribution_by_organism.csv")

    by_genus = pooled.groupby("genus").agg(
        pooled_rows=("locus_tag", "size"),
        distinct_genes=("locus_tag", "nunique"),
        organisms=("organism", "nunique"),
        experiments=("experiment_id", "nunique"),
    ).sort_values("pooled_rows", ascending=False)
    by_genus.to_csv(DATA_DIR / "02_distribution_by_genus.csv")

    log("\n=== distribution by GENUS (before dedup) ===")
    log(by_genus.to_string())
    log("\n=== distribution by ORGANISM (before dedup) ===")
    log(by_org.to_string())

    logf.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
