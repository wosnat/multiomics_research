"""Map each anchor locus_tag to its eggnog Cyanobacteria-level ortholog group(s).

This expands the per-strain Cyanorak-anchored gene set into a strain-agnostic
ortholog representation, which becomes the unit of cross-strain comparison.

Per step-1 method choice, the ortholog tier is eggnog at
`taxonomic_level='Cyanobacteria'` (specificity_rank=2). Lineage-specific
paralogs (e.g. Pro amt1 vs Syn amt1) sit in different Cyanobacteria-level
groups; both lineages are captured because the anchor set was seeded from
both Cyanorak-annotated Prochlorococcus and Synechococcus strains.

Inputs:
    data/01_anchor_genes.csv (375 unique anchor locus_tags)
Outputs:
    data/04_anchor_to_ortholog.csv         one row per (anchor locus_tag, ortholog group_id)
    data/04_anchors_without_group.csv      anchors that returned no Cyanobacteria-level group
    data/04_ortholog_groups.csv            one row per unique group (group_id, gene name, member counts)
    data/04_anchor_to_ortholog.log
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import gene_homologs

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
LOG_PATH = OUT_DIR / "04_anchor_to_ortholog.log"

BATCH_SIZE = 50  # gene_homologs accepts list input; batch to keep responses small.


def fetch_groups_for_batch(locus_tags: list[str]) -> tuple[pd.DataFrame, list[str]]:
    """Return (per-(gene,group) rows, list of anchor locus_tags with NO group)."""
    result = gene_homologs(
        locus_tags=locus_tags,
        source="eggnog",
        taxonomic_level="Cyanobacteria",
        limit=None,
        verbose=True,  # to capture group metadata
    )
    no_groups = result.get("no_groups", [])
    not_found = result.get("not_found", [])
    if not_found:
        logging.warning("Locus tags not found in KG: %s", not_found)
    rows = result.get("results", [])
    if not rows:
        return pd.DataFrame(), no_groups + not_found
    df = pd.DataFrame(rows)
    keep_cols = [
        "locus_tag", "organism_name", "group_id", "consensus_gene_name",
        "consensus_product", "taxonomic_level", "source", "specificity_rank",
        "member_count", "organism_count", "genera", "has_cross_genus_members",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df[keep_cols]
    # list column "genera" → joined string for CSV safety
    if "genera" in df.columns:
        df["genera"] = df["genera"].apply(
            lambda v: " | ".join(v) if isinstance(v, list) else (v if isinstance(v, str) else "")
        )
    return df, no_groups + not_found


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    anchor = pd.read_csv(OUT_DIR / "01_anchor_genes.csv")
    unique_loci = sorted(anchor["locus_tag"].dropna().unique().tolist())
    logging.info("Unique anchor locus_tags: %d (from %d anchor rows across %d strains)",
                 len(unique_loci), len(anchor), anchor["strain"].nunique())

    all_mappings: list[pd.DataFrame] = []
    all_no_group: list[str] = []
    for i in range(0, len(unique_loci), BATCH_SIZE):
        batch = unique_loci[i:i + BATCH_SIZE]
        df, no_group = fetch_groups_for_batch(batch)
        all_mappings.append(df)
        all_no_group.extend(no_group)
        logging.info("Batch %d/%d: %d locus_tags in, %d group rows out, %d anchors had no Cyanobacteria-level group",
                     i // BATCH_SIZE + 1,
                     (len(unique_loci) + BATCH_SIZE - 1) // BATCH_SIZE,
                     len(batch), len(df), len(no_group))

    mapping = pd.concat(all_mappings, ignore_index=True) if all_mappings else pd.DataFrame()
    mapping.to_csv(OUT_DIR / "04_anchor_to_ortholog.csv", index=False)

    # Anchors without a Cyanobacteria-level group: include their anchor metadata for review.
    anchor_meta = anchor.drop_duplicates("locus_tag")[
        ["locus_tag", "strain", "gene_name", "product", "gene_category"]
    ]
    no_group_df = anchor_meta[anchor_meta["locus_tag"].isin(set(all_no_group))].copy()
    no_group_df.to_csv(OUT_DIR / "04_anchors_without_group.csv", index=False)

    # Per-group summary (one row per unique group_id), with consensus gene name and member counts.
    if not mapping.empty:
        group_summary = (
            mapping.drop_duplicates("group_id")[
                ["group_id", "consensus_gene_name", "consensus_product",
                 "taxonomic_level", "source", "specificity_rank",
                 "member_count", "organism_count", "genera", "has_cross_genus_members"]
            ]
            .sort_values("group_id")
            .reset_index(drop=True)
        )
        # How many anchor locus_tags map to each group?
        anchor_per_group = mapping.groupby("group_id")["locus_tag"].nunique().rename("n_anchor_loci_in_group")
        group_summary = group_summary.merge(anchor_per_group, on="group_id", how="left")
        group_summary.to_csv(OUT_DIR / "04_ortholog_groups.csv", index=False)
    else:
        group_summary = pd.DataFrame()

    logging.info("---")
    logging.info("Anchor locus_tags mapped to a Cyanobacteria-level group: %d / %d",
                 mapping["locus_tag"].nunique() if not mapping.empty else 0, len(unique_loci))
    logging.info("Anchor locus_tags WITHOUT a Cyanobacteria-level group: %d", len(no_group_df))
    if len(no_group_df) > 0:
        for _, row in no_group_df.head(20).iterrows():
            logging.info("  %s (%s, %s, %s) — %s", row["locus_tag"], row["strain"],
                         row.get("gene_name") or "", row.get("gene_category") or "", row.get("product") or "")
        if len(no_group_df) > 20:
            logging.info("  ... %d more", len(no_group_df) - 20)
    logging.info("Unique ortholog groups: %d", len(group_summary))
    if not group_summary.empty:
        cross_genus_n = (group_summary["has_cross_genus_members"] == "cross_genus").sum()
        single_genus_n = (group_summary["has_cross_genus_members"] == "single_genus").sum()
        logging.info("  cross-genus groups: %d ; single-genus groups: %d", cross_genus_n, single_genus_n)
        logging.info("  median member_count: %.0f, max: %d",
                     group_summary["member_count"].median(),
                     int(group_summary["member_count"].max()))
    logging.info("Wrote %s, %s, %s",
                 OUT_DIR / "04_anchor_to_ortholog.csv",
                 OUT_DIR / "04_anchors_without_group.csv",
                 OUT_DIR / "04_ortholog_groups.csv")


if __name__ == "__main__":
    main()
