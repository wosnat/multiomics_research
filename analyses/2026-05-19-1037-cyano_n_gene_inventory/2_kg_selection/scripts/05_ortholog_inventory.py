"""Build the 19-strain × N(ortholog-group) presence/copy-number inventory.

For each unique ortholog group from 04_ortholog_groups.csv, pull member
genes across all 19 cyano strains in scope (via genes_by_homolog_group).
Result is a long-format member table and a wide-format copy-number matrix.

Inputs:
    data/04_ortholog_groups.csv  (54 groups)
    data/02_strain_table.csv     (19 strains in scope)
Outputs:
    data/05_inventory_members.csv     long format: one row per (group, strain, locus_tag)
    data/05_inventory_matrix.csv      wide: rows = strain, cols = group_id, cells = copy count
    data/05_ortholog_inventory.log
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import genes_by_homolog_group

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
LOG_PATH = OUT_DIR / "05_ortholog_inventory.log"


def fetch_group_members(group_id: str) -> pd.DataFrame:
    result = genes_by_homolog_group(
        group_ids=[group_id],
        limit=None,
        verbose=False,
    )
    rows = result.get("results", [])
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    return df[["group_id", "organism_name", "locus_tag", "gene_name", "product", "gene_category"]]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    groups = pd.read_csv(OUT_DIR / "04_ortholog_groups.csv")
    strain_table = pd.read_csv(OUT_DIR / "02_strain_table.csv")
    scope_strains = strain_table["strain"].tolist()
    logging.info("Pulling members for %d ortholog groups", len(groups))

    member_frames: list[pd.DataFrame] = []
    for i, group_id in enumerate(groups["group_id"], start=1):
        df = fetch_group_members(group_id)
        member_frames.append(df)
        if i % 10 == 0 or i == len(groups):
            logging.info("  fetched %d/%d groups", i, len(groups))

    all_members = pd.concat(member_frames, ignore_index=True)
    logging.info("Total member rows across all groups: %d", len(all_members))
    logging.info("Distinct organisms across all group memberships: %d", all_members["organism_name"].nunique())

    # Restrict to in-scope strains (we don't care about Alteromonas/Marinobacter members for this analysis).
    in_scope = all_members[all_members["organism_name"].isin(scope_strains)].copy()
    out_of_scope_orgs = sorted(set(all_members["organism_name"]) - set(scope_strains))
    logging.info("In-scope (cyano-only) member rows: %d", len(in_scope))
    if out_of_scope_orgs:
        logging.info("Out-of-scope organisms present in groups (will be filtered out): %s", out_of_scope_orgs)
    in_scope.to_csv(OUT_DIR / "05_inventory_members.csv", index=False)

    # Wide matrix: rows = strain (in canonical order), cols = group_id, cells = copy count (n locus_tags in that strain × group).
    copy_counts = in_scope.groupby(["organism_name", "group_id"])["locus_tag"].nunique().unstack(fill_value=0)
    # Ensure every scope strain and every group appears as a row/column even if zero.
    copy_counts = copy_counts.reindex(index=scope_strains, columns=groups["group_id"], fill_value=0)
    copy_counts.index.name = "strain"
    copy_counts.to_csv(OUT_DIR / "05_inventory_matrix.csv")

    # Summary stats.
    logging.info("---")
    logging.info("Inventory matrix shape: %s strains × %s groups", *copy_counts.shape)
    per_strain_n_groups_present = (copy_counts > 0).sum(axis=1)
    per_strain_total_copies = copy_counts.sum(axis=1)
    logging.info("Per-strain N(groups present) summary:")
    for strain, n in per_strain_n_groups_present.items():
        total = int(per_strain_total_copies[strain])
        logging.info("  %-55s %2d groups present, %2d total gene copies",
                     strain, int(n), total)

    per_group_n_strains = (copy_counts > 0).sum(axis=0)
    logging.info("Group ubiquity distribution (n strains where present):")
    ubiquity_counts = per_group_n_strains.value_counts().sort_index()
    for n_strains, n_groups in ubiquity_counts.items():
        logging.info("  %2d strains: %d groups", int(n_strains), int(n_groups))

    logging.info("Wrote %s and %s", OUT_DIR / "05_inventory_members.csv", OUT_DIR / "05_inventory_matrix.csv")


if __name__ == "__main__":
    main()
