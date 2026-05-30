#!/usr/bin/env python
"""
qc_second_dedup.py — EXPLORATORY: what does a second, broader ortholog dedup do?

Current (step 4) dedup key = finest-available group (cyanorak rank0 / eggNOG
family rank1 / phylum rank2 / COG rank3), applied WITHIN pool. This explores
adding a SECOND dedup at a broader shared level. For each candidate broad level
we report seeds-before -> seeds-after, both within-pool and cross-pool (global).

Broad-level candidates per gene (from gene_homologs):
  family  = finest eggNOG with specificity_rank == 1   (e.g. Alteromonadaceae,
            Prochlorococcaceae) — genus/family scope, rarely cross-pool
  phylum  = eggNOG rank == 2
  cog     = eggNOG rank == 3 (Bacteria-level COG) — collapses across ALL pools
Genes lacking the chosen level keep their finest key (no over-merge).

Read-only: prints a comparison table + example cross-pool collapses. Writes nothing.

Run from scripts/ dir:
  env -u VIRTUAL_ENV uv run python qc_second_dedup.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from multiomics_explorer import GraphConnection, gene_homologs

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
pd.set_option("display.width", 200)


def batched(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def main() -> int:
    reps = pd.read_csv(DATA_DIR / "04_pool_representatives.csv")
    loci = reps["rep_locus_tag"].tolist()

    # fetch all groups for the representative genes
    per_gene: dict[str, list[dict]] = {lt: [] for lt in loci}
    with GraphConnection() as conn:
        for chunk in batched(loci, 150):
            for r in gene_homologs(locus_tags=chunk, verbose=True, limit=None, conn=conn)["results"]:
                per_gene.setdefault(r["locus_tag"], []).append(r)

    def group_at(lt, want_rank):
        """eggNOG group_id at the requested specificity_rank, else None."""
        for r in per_gene.get(lt, []):
            if r.get("source") == "eggnog" and r.get("specificity_rank") == want_rank:
                return r["group_id"]
        return None

    reps = reps.copy()
    reps["finest_key"] = reps["dedup_key"]
    for name, rank in [("family", 1), ("phylum", 2), ("cog", 3)]:
        reps[name + "_key"] = [
            group_at(lt, rank) or fk
            for lt, fk in zip(reps["rep_locus_tag"], reps["finest_key"])
        ]

    print(f"current seeds (finest, within-pool): {len(reps)}")
    print(f"  per pool: {reps['pool'].value_counts().to_dict()}\n")

    rows = []
    for level in ["finest", "family", "phylum", "cog"]:
        key = level + "_key" if level != "finest" else "finest_key"
        within = reps.groupby(["pool", key]).ngroups
        cross = reps[key].nunique()
        rows.append({"second_level": level,
                     "seeds_within_pool": within,
                     "seeds_cross_pool": cross})
    comp = pd.DataFrame(rows)
    print("=== seeds remaining after a SECOND dedup at each level ===")
    print(comp.to_string(index=False))

    # show cross-pool collapses at COG level (the aggressive one)
    print("\n=== example CROSS-POOL collapses at COG (Bacteria) level ===")
    g = reps.groupby("cog_key")
    shown = 0
    for key, grp in g:
        pools = grp["pool"].nunique()
        if pools > 1 and len(grp) > 1 and str(key).startswith("eggnog"):
            names = grp["rep_gene_name"].fillna(grp["rep_product"].str.slice(0, 20)).tolist()
            print(f"  {key}: {len(grp)} seeds across {pools} pools "
                  f"[{', '.join(sorted(grp['pool'].unique()))}] -> {names[:6]}")
            shown += 1
            if shown >= 12:
                break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
