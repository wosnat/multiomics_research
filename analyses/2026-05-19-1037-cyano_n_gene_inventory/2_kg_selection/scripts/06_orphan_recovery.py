"""Recover orphan anchors into the inventory as synthetic singleton ortholog groups.

Background: 9 anchor locus_tags returned no eggnog Cyanobacteria-level group
(see 04_anchors_without_group.csv). Initial decision was to leave them out
for ortholog-tier purity; researcher overrode — don't throw out the orphans.

Recovery strategy: for each orphan, create a synthetic ortholog-group entry
with group_id pattern `anchor_singleton:<gene_name>:<strain_short>` (clearly
labeled as non-eggnog). When a strain has multiple orphan copies of the same
gene (e.g., WH8102's SYNW2462 + SYNW2463 both nrtP), collapse them into one
synthetic group with copy_count = 2.

Inputs (read; overwritten in place with augmented versions):
    data/04_ortholog_groups.csv           (54 groups → 61 after augmentation)
    data/05_inventory_members.csv         (549 rows → 558 after augmentation)
    data/05_inventory_matrix.csv          (19 × 54 → 19 × 61 after augmentation)
    data/04_anchors_without_group.csv     (9 orphan rows)
    data/01_anchor_genes.csv              (for product, gene_category metadata)
Outputs:
    Overwrites the four files above with the augmented versions.
    data/06_orphan_recovery.log
    data/06_synthetic_groups.csv          (audit: just the 7 added groups)
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
LOG_PATH = OUT_DIR / "06_orphan_recovery.log"


STRAIN_SHORT = {
    "Prochlorococcus AS9601": "AS9601",
    "Prochlorococcus MED4": "MED4",
    "Prochlorococcus MIT0801": "MIT0801",
    "Prochlorococcus MIT9301": "MIT9301",
    "Prochlorococcus MIT9303": "MIT9303",
    "Prochlorococcus MIT9312": "MIT9312",
    "Prochlorococcus MIT9313": "MIT9313",
    "Prochlorococcus NATL1A": "NATL1A",
    "Prochlorococcus NATL2A": "NATL2A",
    "Prochlorococcus RSP50": "RSP50",
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)": "SS120",
    "Synechococcus WH8102": "WH8102",
    "Synechococcus WH7803": "WH7803",
    "Synechococcus CC9311": "CC9311",
    "Synechococcus sp. BL107": "BL107",
    "Synechococcus PCC 7002": "PCC7002",
    "Synechococcus elongatus PCC 7942": "PCC7942",
    "Synechococcus elongatus UTEX 2973": "UTEX2973",
    "Thermosynechococcus vestitus BP-1": "BP-1",
}


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    orphans = pd.read_csv(OUT_DIR / "04_anchors_without_group.csv")
    anchor = pd.read_csv(OUT_DIR / "01_anchor_genes.csv")
    groups = pd.read_csv(OUT_DIR / "04_ortholog_groups.csv")
    members = pd.read_csv(OUT_DIR / "05_inventory_members.csv")
    matrix = pd.read_csv(OUT_DIR / "05_inventory_matrix.csv", index_col=0)

    logging.info("Inputs: %d orphan loci, %d existing groups, %d member rows, matrix %s",
                 len(orphans), len(groups), len(members), matrix.shape)

    # 04_anchors_without_group.csv already has product + gene_category for each orphan;
    # no merge needed. Just fill missing gene_name with locus_tag for group-id naming.
    orphan_meta = orphans.copy()
    orphan_meta["gene_name_filled"] = orphan_meta["gene_name"].fillna(orphan_meta["locus_tag"])

    synthetic_rows = []
    new_member_rows = []
    for (strain, gene_name_filled), block in orphan_meta.groupby(["strain", "gene_name_filled"]):
        loci = block["locus_tag"].tolist()
        # safe slug for group_id
        gene_slug = re.sub(r"[^A-Za-z0-9_.-]", "_", str(gene_name_filled))
        strain_slug = STRAIN_SHORT.get(strain, re.sub(r"[^A-Za-z0-9_.-]", "_", strain))
        group_id = f"anchor_singleton:{gene_slug}:{strain_slug}"
        product = block["product"].dropna().iloc[0] if block["product"].notna().any() else ""
        category = block["gene_category"].dropna().iloc[0] if block["gene_category"].notna().any() else ""
        synthetic_rows.append({
            "group_id": group_id,
            "consensus_gene_name": gene_name_filled,
            "consensus_product": product,
            "taxonomic_level": "anchor_singleton",
            "source": "anchor_singleton",
            "specificity_rank": pd.NA,
            "member_count": len(loci),
            "organism_count": 1,
            "genera": strain.split()[0],
            "has_cross_genus_members": "single_genus",
            "n_anchor_loci_in_group": len(loci),
        })
        for locus in loci:
            new_member_rows.append({
                "group_id": group_id,
                "organism_name": strain,
                "locus_tag": locus,
                "gene_name": gene_name_filled,
                "product": product,
                "gene_category": category,
            })

    synthetic_groups_df = pd.DataFrame(synthetic_rows)
    new_members_df = pd.DataFrame(new_member_rows)
    synthetic_groups_df.to_csv(OUT_DIR / "06_synthetic_groups.csv", index=False)

    logging.info("Synthetic singleton groups added: %d (one per orphan (strain, gene_name) pair)",
                 len(synthetic_groups_df))
    for _, row in synthetic_groups_df.iterrows():
        logging.info("  + %s — %s (%s) [member_count=%d]",
                     row["group_id"], row["consensus_gene_name"],
                     row["consensus_product"], int(row["member_count"]))

    # Augment 04_ortholog_groups.csv (append synthetic; preserve column order)
    augmented_groups = pd.concat([groups, synthetic_groups_df], ignore_index=True)
    augmented_groups.to_csv(OUT_DIR / "04_ortholog_groups.csv", index=False)

    # Augment 05_inventory_members.csv (append new member rows)
    augmented_members = pd.concat([members, new_members_df], ignore_index=True)
    augmented_members.to_csv(OUT_DIR / "05_inventory_members.csv", index=False)

    # Rebuild 05_inventory_matrix.csv from augmented members.
    copy_counts = augmented_members.groupby(["organism_name", "group_id"])["locus_tag"].nunique().unstack(fill_value=0)
    copy_counts = copy_counts.reindex(index=matrix.index.tolist(), columns=augmented_groups["group_id"], fill_value=0)
    copy_counts.index.name = "strain"
    copy_counts.to_csv(OUT_DIR / "05_inventory_matrix.csv")

    # Report deltas
    logging.info("---")
    logging.info("Augmented groups CSV: %d → %d rows", len(groups), len(augmented_groups))
    logging.info("Augmented members CSV: %d → %d rows", len(members), len(augmented_members))
    logging.info("Augmented matrix: %s → %s", matrix.shape, copy_counts.shape)
    # Per-strain inventory delta
    logging.info("Per-strain copy-count change:")
    for strain in matrix.index:
        before = int(matrix.loc[strain].sum())
        after = int(copy_counts.loc[strain].sum())
        delta = after - before
        marker = f"  +{delta}" if delta > 0 else f"  {delta}" if delta < 0 else ""
        logging.info("  %-55s before=%d  after=%d%s", strain, before, after, marker)


if __name__ == "__main__":
    main()
