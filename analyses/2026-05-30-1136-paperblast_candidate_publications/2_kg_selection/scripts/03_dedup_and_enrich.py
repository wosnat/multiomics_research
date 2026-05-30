#!/usr/bin/env python
"""
03_dedup_and_enrich.py — collapse the pooled top-responders to one row per gene,
add cross-experiment frequency, and enrich with annotation / ortholog signals.

Input:  ../data/02_pooled_top_responders.csv  (one row per experiment x gene)
Step 1 of the user's next-steps list: dedup per locus_tag (don't care about dups).

Per gene we keep:
  - freq            = # distinct experiments where it was a top responder
  - best_rank       = min best_rank across those experiments
  - n_genus, n_org  = breadth context (always 1 here; gene is org-specific)
  - organism, genus, gene_name, product, gene_category
Enriched from gene_overview (batched):
  - annotation_state          (no_evidence / catch_all_only / informative_single / informative_multi)
  - annotation_quality        (0..3)
  - closest_ortholog_group_size, closest_ortholog_genera  (BLAST-success proxies)
  - has_seq                   (PaperBLAST eligibility)

Output: ../data/03_genes_dedup.csv  (one row per locus_tag)

Run from scripts/ dir:
  env -u VIRTUAL_ENV uv run python 03_dedup_and_enrich.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from multiomics_explorer import GraphConnection, gene_aa_sequence, gene_overview

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
pd.set_option("display.width", 220)
pd.set_option("display.max_rows", None)


def batched(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def main() -> int:
    pooled = pd.read_csv(DATA_DIR / "02_pooled_top_responders.csv")

    # ---- dedup per locus_tag ----
    grp = pooled.groupby("locus_tag")
    dedup = grp.agg(
        organism=("organism", "first"),
        genus=("genus", "first"),
        gene_name=("gene_name", "first"),
        product=("product", "first"),
        gene_category=("gene_category", "first"),
        freq=("experiment_id", "nunique"),
        best_rank=("best_rank", "min"),
    ).reset_index()
    dedup = dedup.sort_values(["freq", "best_rank"], ascending=[False, True])

    print(f"pooled rows: {len(pooled)} -> dedup genes: {len(dedup)}")

    # ---- enrich via gene_overview + gene_aa_sequence (batched) ----
    loci = dedup["locus_tag"].tolist()
    ov_rows, seq_rows = {}, {}
    with GraphConnection() as conn:
        for chunk in batched(loci, 200):
            for r in gene_overview(locus_tags=chunk, limit=None, conn=conn)["results"]:
                ov_rows[r["locus_tag"]] = r
            sres = gene_aa_sequence(locus_tags=chunk, limit=None, conn=conn)["results"]
            for r in sres:
                seq_rows[r["locus_tag"]] = r

    def ov(lt, key):
        return ov_rows.get(lt, {}).get(key)

    dedup["annotation_state"] = dedup["locus_tag"].map(lambda lt: ov(lt, "annotation_state"))
    dedup["annotation_quality"] = dedup["locus_tag"].map(lambda lt: ov(lt, "annotation_quality"))
    dedup["ortholog_group_size"] = dedup["locus_tag"].map(
        lambda lt: ov(lt, "closest_ortholog_group_size"))
    dedup["ortholog_genera"] = dedup["locus_tag"].map(
        lambda lt: ", ".join(ov(lt, "closest_ortholog_genera") or []))
    dedup["protein_id"] = dedup["locus_tag"].map(lambda lt: seq_rows.get(lt, {}).get("protein_id"))
    dedup["has_seq"] = dedup["locus_tag"].map(
        lambda lt: bool(seq_rows.get(lt, {}).get("sequence")))

    out = DATA_DIR / "03_genes_dedup.csv"
    dedup.to_csv(out, index=False)
    print(f"wrote {out}  ({len(dedup)} genes)")

    # ---- quick distributions for review ----
    print("\n=== annotation_state distribution ===")
    print(dedup["annotation_state"].value_counts(dropna=False).to_string())
    print("\n=== has_seq ===")
    print(dedup["has_seq"].value_counts().to_string())
    print("\n=== freq distribution (genes appearing in N experiments) ===")
    print(dedup["freq"].value_counts().sort_index(ascending=False).to_string())
    print("\n=== by genus: genes, median ortholog_group_size, %informative ===")
    g = dedup.assign(
        informative=dedup["annotation_state"].isin(["informative_single", "informative_multi"]))
    summ = g.groupby("genus").agg(
        genes=("locus_tag", "size"),
        median_og_size=("ortholog_group_size", "median"),
        pct_informative=("informative", "mean"),
        pct_has_seq=("has_seq", "mean"),
    ).sort_values("genes", ascending=False)
    summ["pct_informative"] = (summ["pct_informative"] * 100).round(0)
    summ["pct_has_seq"] = (summ["pct_has_seq"] * 100).round(0)
    print(summ.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
