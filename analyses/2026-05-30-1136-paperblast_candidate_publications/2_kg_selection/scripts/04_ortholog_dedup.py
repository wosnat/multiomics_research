#!/usr/bin/env python
"""
04_ortholog_dedup.py — assign 4 pools, attach finest-available ortholog group,
and dedup genes within each pool at the ortholog-group level.

Input:  ../data/03_genes_dedup.csv   (1653 genes, one row per locus_tag)

Pools (by genus):
  pro          = Prochlorococcus
  syn          = Synechococcus + Parasynechococcus + Picosynechococcus (+ any *synechococcus)
  alt          = Alteromonas
  other_hetero = everything else (Marinobacter, Ruegeria, Pseudomonas, Shewanella, ...)

Dedup key (uniform "finest available", locked 2026-05-30):
  per gene, take the ortholog group with the SMALLEST specificity_rank
  (cyanorak curated = rank 0; eggNOG family=1, phylum=2, COG/Bacteria=3),
  tie-break cyanorak over eggnog. Genes with no group -> key = locus_tag (no collapse).

Dedup is WITHIN pool: genes sharing a dedup_key in the same pool collapse to one
representative. Representative = BEST BLAST ODDS (locked 2026-05-30): the member
with the largest ortholog_group_size, then most ortholog genera, then highest
freq, then best rank. Rep-choice only affects which sequence is BLASTed; the
group's response score still aggregates across ALL members (freq = max member
freq; best_rank = min member best_rank). Collapsed members are recorded.
Sequences are dropped per request (locus_tag is the PaperBLAST input).

Outputs (../data/):
  04_genes_ortholog.csv   all genes with pool + dedup_key + group metadata
  04_pool_representatives.csv   one row per (pool, dedup_key): the kept seed + members
  04_ortholog_dedup.log

Run from scripts/ dir:
  env -u VIRTUAL_ENV uv run python 04_ortholog_dedup.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from multiomics_explorer import GraphConnection, gene_homologs

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
pd.set_option("display.width", 220)
pd.set_option("display.max_rows", None)

SYN_GENERA = {"synechococcus", "parasynechococcus", "picosynechococcus",
              "cyanobium", "thermosynechococcus"}


def pool_of(genus: str) -> str:
    g = (genus or "").lower()
    if "prochlorococcus" in g:
        return "pro"
    if g in SYN_GENERA or "synechococcus" in g:
        return "syn"
    if "alteromonas" in g:
        return "alt"
    return "other_hetero"


def batched(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def finest_group(rows: list[dict]) -> dict | None:
    """Pick the finest ortholog group for one gene: smallest specificity_rank,
    cyanorak preferred over eggnog on ties."""
    if not rows:
        return None
    def key(r):
        src_pref = 0 if r.get("source") == "cyanorak" else 1
        return (r.get("specificity_rank", 99), src_pref)
    return sorted(rows, key=key)[0]


def main() -> int:
    genes = pd.read_csv(DATA_DIR / "03_genes_dedup.csv")
    genes["pool"] = genes["genus"].map(pool_of)

    logf = open(DATA_DIR / "04_ortholog_dedup.log", "w")

    def log(m: str) -> None:
        print(m)
        logf.write(m + "\n")

    log("== ortholog dedup: finest-available key, within-pool ==")
    log("\npool sizes (before ortholog dedup):")
    log(genes["pool"].value_counts().to_string())

    # ---- fetch ortholog groups (batched) ----
    loci = genes["locus_tag"].tolist()
    per_gene: dict[str, list[dict]] = {lt: [] for lt in loci}
    with GraphConnection() as conn:
        for chunk in batched(loci, 150):
            res = gene_homologs(locus_tags=chunk, verbose=True, limit=None, conn=conn)
            for r in res["results"]:
                per_gene.setdefault(r["locus_tag"], []).append(r)

    fg = {lt: finest_group(rows) for lt, rows in per_gene.items()}
    genes["dedup_key"] = genes["locus_tag"].map(
        lambda lt: (fg[lt]["group_id"] if fg.get(lt) else f"locus:{lt}"))
    genes["og_source"] = genes["locus_tag"].map(
        lambda lt: (fg[lt]["source"] if fg.get(lt) else None))
    genes["og_level"] = genes["locus_tag"].map(
        lambda lt: (fg[lt]["taxonomic_level"] if fg.get(lt) else None))
    genes["og_rank"] = genes["locus_tag"].map(
        lambda lt: (fg[lt]["specificity_rank"] if fg.get(lt) else None))
    genes["og_member_count"] = genes["locus_tag"].map(
        lambda lt: (fg[lt].get("member_count") if fg.get(lt) else None))

    n_nogroup = int((genes["dedup_key"].str.startswith("locus:")).sum())
    log(f"\ngenes with an ortholog group: {len(genes) - n_nogroup} / {len(genes)}  "
        f"({n_nogroup} keyed by locus_tag)")
    log("\nfinest-group source mix:")
    log(genes["og_source"].value_counts(dropna=False).to_string())

    genes.to_csv(DATA_DIR / "04_genes_ortholog.csv", index=False)

    # ---- within-pool ortholog dedup ----
    # Representative = BEST BLAST ODDS: largest ortholog_group_size, then most
    # genera, then highest freq, then best (lowest) rank. The GROUP's response
    # score aggregates across ALL members (freq=max, best_rank=min) so rep-choice
    # doesn't change ranking, only which sequence is BLASTed.
    genes["_n_genera"] = genes["ortholog_genera"].fillna("").apply(
        lambda s: len([g for g in str(s).split(",") if g.strip()]))
    reps = []
    for (pool, key), grp in genes.groupby(["pool", "dedup_key"], sort=False):
        rep = grp.sort_values(
            ["ortholog_group_size", "_n_genera", "freq", "best_rank"],
            ascending=[False, False, False, True],
        ).iloc[0]
        members = grp["locus_tag"].tolist()
        reps.append({
            "pool": pool,
            "dedup_key": key,
            "rep_locus_tag": rep["locus_tag"],
            "rep_organism": rep["organism"],
            "rep_gene_name": rep["gene_name"],
            "rep_product": rep["product"],
            "rep_annotation_state": rep["annotation_state"],
            "rep_protein_id": rep["protein_id"],
            # group-level response score (aggregated across all members)
            "freq": int(grp["freq"].max()),
            "best_rank": float(grp["best_rank"].min()),
            "og_source": rep["og_source"],
            "og_level": rep["og_level"],
            "n_collapsed": len(members),
            "collapsed_members": ";".join(m for m in members if m != rep["locus_tag"]),
        })
    reps_df = pd.DataFrame(reps).sort_values(
        ["pool", "freq", "best_rank"], ascending=[True, False, True])
    reps_df.to_csv(DATA_DIR / "04_pool_representatives.csv", index=False)

    log("\n=== seeds per pool: before -> after ortholog dedup ===")
    before = genes.groupby("pool").size()
    after = reps_df.groupby("pool").size()
    summary = pd.DataFrame({"genes_before": before, "seeds_after": after}).fillna(0).astype(int)
    summary["collapsed"] = summary["genes_before"] - summary["seeds_after"]
    log(summary.to_string())
    log(f"\nTOTAL seeds after dedup: {len(reps_df)}  (from {len(genes)} genes)")

    log("\n=== multi-gene collapses (n_collapsed > 1), top 15 ===")
    multi = reps_df[reps_df["n_collapsed"] > 1].head(15)
    for _, r in multi.iterrows():
        log(f"  [{r['pool']}] {r['rep_locus_tag']} ({r['rep_gene_name'] or r['rep_product'][:30]}) "
            f"<- {r['n_collapsed']} genes via {r['dedup_key']}")
    logf.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
