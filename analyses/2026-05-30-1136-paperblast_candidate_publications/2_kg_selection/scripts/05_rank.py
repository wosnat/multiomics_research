#!/usr/bin/env python
"""
05_rank.py — score & rank the deduped seeds.

Two modes (the project has two goals; we rank for one at a time):
  --mode papers   (default, GOAL 1: more papers) score = response x conservation.
                  info_gain is DROPPED so the ranking is not pulled toward
                  hypotheticals; it favors strongly-responding, broadly-conserved
                  proteins = best PaperBLAST hit odds.
  --mode jackpot  (GOAL 2: info on low-info genes) score = response x
                  conservation x info_gain (boosts no_evidence/catch_all_only).

Multiplicative so a seed must score on every active factor to rank high.

Components (each normalized to [0,1], computed globally across all seeds):
  response     = RESP_RANK_W*rank_norm + RESP_FREQ_W*freq_norm
                   rank_norm = (11 - best_rank)/10   (best_rank 1->1.0, 10->0.1)
                   freq_norm = freq / max_freq        (cross-experiment reproducibility)
  conservation = 0.5*ogsize_norm + 0.5*genera_norm   (PaperBLAST hit-likelihood proxy)
                   ogsize_norm = log1p(ortholog_group_size)/log1p(max)
                   genera_norm = n_ortholog_genera / max_genera
  info_gain    = map(annotation_state): no_evidence 1.0, catch_all_only 0.75,
                   informative_single 0.35, informative_multi 0.15

Round-one selection (top N per pool) applies two soft caps while walking the
score-sorted pool, both skip-to-next when full:
  --cat-cap   max seeds per gene_category  (default 2; forces multi-pathway breadth)
  --org-cap   max seeds per organism       (default 3; stops one strain, e.g.
              Marinobacter, from filling the entire other_hetero pool)

Inputs:  ../data/04_pool_representatives.csv, ../data/04_genes_ortholog.csv
Outputs: ../data/05_ranked_seeds_<mode>.csv  (all seeds ranked)
         ../data/05_selected_<mode>.csv      (the round-one top-N per pool)
Run from scripts/ dir:
  env -u VIRTUAL_ENV uv run python 05_rank.py --mode papers --top 8 --cat-cap 2 --org-cap 3
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
pd.set_option("display.width", 240)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_colwidth", 40)

INFO_GAIN = {
    "no_evidence": 1.0,
    "catch_all_only": 0.75,
    "informative_single": 0.35,
    "informative_multi": 0.15,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["papers", "jackpot"], default="papers")
    ap.add_argument("--resp-rank-w", type=float, default=0.6)
    ap.add_argument("--resp-freq-w", type=float, default=0.4)
    ap.add_argument("--top", type=int, default=8, help="seeds to select per pool")
    ap.add_argument("--cat-cap", type=int, default=2,
                    help="max seeds per gene_category within a pool's selection (0 = no cap)")
    ap.add_argument("--org-cap", type=int, default=3,
                    help="max seeds per organism within a pool's selection (0 = no cap)")
    args = ap.parse_args()

    reps = pd.read_csv(DATA_DIR / "04_pool_representatives.csv")
    full = pd.read_csv(DATA_DIR / "04_genes_ortholog.csv").set_index("locus_tag")

    # attach conservation attributes from the full table (by representative gene)
    reps["ortholog_group_size"] = reps["rep_locus_tag"].map(full["ortholog_group_size"])
    reps["ortholog_genera"] = reps["rep_locus_tag"].map(full["ortholog_genera"]).fillna("")
    reps["annotation_state"] = reps["rep_annotation_state"]
    reps["n_genera"] = reps["ortholog_genera"].apply(
        lambda s: len([g for g in str(s).split(",") if g.strip()]) if s else 0)

    # ---- components ----
    reps["rank_norm"] = (11 - reps["best_rank"]).clip(lower=0) / 10.0
    reps["freq_norm"] = reps["freq"] / reps["freq"].max()
    reps["response"] = args.resp_rank_w * reps["rank_norm"] + args.resp_freq_w * reps["freq_norm"]

    ogs = reps["ortholog_group_size"].fillna(1).clip(lower=1)
    reps["ogsize_norm"] = ogs.apply(math.log1p) / math.log1p(ogs.max())
    reps["genera_norm"] = reps["n_genera"] / max(reps["n_genera"].max(), 1)
    reps["conservation"] = 0.5 * reps["ogsize_norm"] + 0.5 * reps["genera_norm"]

    reps["info_gain"] = reps["annotation_state"].map(INFO_GAIN).fillna(0.5)

    if args.mode == "papers":
        reps["score"] = reps["response"] * reps["conservation"]
    else:  # jackpot
        reps["score"] = reps["response"] * reps["conservation"] * reps["info_gain"]

    # gene_category for the diversity cap (from the full enriched table)
    reps["gene_category"] = reps["rep_locus_tag"].map(full["gene_category"]).fillna("Unknown")

    reps = reps.sort_values(["pool", "score"], ascending=[True, False])
    reps["pool_rank"] = reps.groupby("pool")["score"].rank(ascending=False, method="first").astype(int)

    # ---- per-pool selection with soft category + organism caps ----
    def select(sub: pd.DataFrame) -> pd.DataFrame:
        chosen, cat_count, org_count = [], {}, {}
        for _, row in sub.iterrows():  # sub is score-sorted desc
            if len(chosen) >= args.top:
                break
            c, o = row["gene_category"], row["rep_organism"]
            if args.cat_cap and cat_count.get(c, 0) >= args.cat_cap:
                continue
            if args.org_cap and org_count.get(o, 0) >= args.org_cap:
                continue
            chosen.append(row)
            cat_count[c] = cat_count.get(c, 0) + 1
            org_count[o] = org_count.get(o, 0) + 1
        return pd.DataFrame(chosen)

    selected = pd.concat(
        [select(reps[reps["pool"] == p]) for p in ["pro", "syn", "alt", "other_hetero"]],
        ignore_index=True,
    )
    selected["sel_rank"] = selected.groupby("pool").cumcount() + 1

    # full ranking (all seeds) + the selected round-one set
    full_cols = ["pool", "pool_rank", "rep_locus_tag", "rep_organism", "rep_gene_name",
                 "rep_product", "gene_category", "annotation_state", "freq", "best_rank",
                 "n_collapsed", "ortholog_group_size", "n_genera", "og_source", "og_level",
                 "response", "conservation", "info_gain", "score", "rep_protein_id", "dedup_key"]
    reps[full_cols].to_csv(DATA_DIR / f"05_ranked_seeds_{args.mode}.csv", index=False)
    sel_cols = ["pool", "sel_rank"] + [c for c in full_cols if c not in ("pool", "pool_rank")]
    selected[sel_cols].to_csv(DATA_DIR / f"05_selected_{args.mode}.csv", index=False)
    print(f"mode={args.mode}  top={args.top}/pool  cat_cap={args.cat_cap}  "
          f"org_cap={args.org_cap}  -> {len(selected)} selected "
          f"(full ranking: {len(reps)} seeds)\n")

    for pool in ["pro", "syn", "alt", "other_hetero"]:
        sub = selected[selected["pool"] == pool]
        n_org = sub["rep_organism"].nunique()
        print(f"\n===== SELECTED — {pool} ({(reps['pool']==pool).sum()} seeds in pool, "
              f"{n_org} organisms) =====")
        show = sub[["sel_rank", "rep_locus_tag", "rep_organism", "rep_gene_name",
                    "rep_product", "gene_category", "annotation_state", "freq",
                    "best_rank", "ortholog_group_size", "n_genera", "score"]].copy()
        show["rep_organism"] = show["rep_organism"].str.slice(0, 22)
        show["rep_gene_name"] = show["rep_gene_name"].fillna("-")
        show["rep_product"] = show["rep_product"].str.slice(0, 30)
        show["gene_category"] = show["gene_category"].str.slice(0, 20)
        show["score"] = show["score"].round(3)
        print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
