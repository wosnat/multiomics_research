"""Cluster the 19 × 61 N-gene inventory + render the primary clustered heatmap.

Reads from step 2's inventory matrix + strain table + group table, applies
the step-4 heatmap_clustering module, writes the clustered (reordered)
matrix + linkage CSVs + the heatmap figure.

Outputs:
    data/01_clustered_matrix.csv        reordered inventory matrix
    data/01_row_linkage.csv             scipy linkage matrix for strains
    data/01_col_linkage.csv             scipy linkage matrix for ortholog groups
    data/01_row_order.txt               final strain leaf order
    data/01_col_order.txt               final group leaf order
    figures/01_clustered_heatmap.png    primary figure
    figures/01_clustered_heatmap.pdf    publication-ready
    data/01_cluster_and_plot.log
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Local methods module (sibling step folder)
HERE = Path(__file__).resolve()
ANALYSIS_ROOT = HERE.parents[2]
sys.path.insert(0, str(ANALYSIS_ROOT / "4_methods"))
from heatmap_clustering import cluster_inventory, render_heatmap  # noqa: E402
from n_pathway_categories import categorize, CATEGORY_PALETTE  # noqa: E402

STEP_DIR = HERE.parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"
INVENTORY_DIR = ANALYSIS_ROOT / "2_kg_selection" / "data"
LOG_PATH = DATA_DIR / "01_cluster_and_plot.log"


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    matrix = pd.read_csv(INVENTORY_DIR / "05_inventory_matrix.csv", index_col=0)
    strain_table = pd.read_csv(INVENTORY_DIR / "02_strain_table.csv")
    groups = pd.read_csv(INVENTORY_DIR / "04_ortholog_groups.csv")
    logging.info("Loaded inventory matrix: %s strains × %s groups", *matrix.shape)
    logging.info("Loaded strain table: %d rows, group table: %d rows", len(strain_table), len(groups))

    # Short strain labels: genus is in annotation column; just strain name + clade
    SHORT = {
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
        "Synechococcus PCC 7002": "PCC 7002",
        "Synechococcus elongatus PCC 7942": "PCC 7942 (elong.)",
        "Synechococcus elongatus UTEX 2973": "UTEX 2973 (elong.)",
        "Thermosynechococcus vestitus BP-1": "BP-1 (Thermosyn.)",
    }
    matrix.index = matrix.index.map(lambda s: SHORT.get(s, s))

    # Cluster
    clustered = cluster_inventory(matrix)
    logging.info("Clustering complete. Strain leaf order:")
    for i, s in enumerate(clustered.row_order, 1):
        logging.info("  %2d. %s", i, s)
    logging.info("First 10 group leaves: %s", clustered.col_order[:10])
    logging.info("Last 10 group leaves: %s", clustered.col_order[-10:])

    # Persist reordered matrix and linkage data
    clustered.matrix.to_csv(DATA_DIR / "01_clustered_matrix.csv")
    pd.DataFrame(clustered.row_linkage, columns=["cluster_a", "cluster_b", "distance", "size"]).to_csv(
        DATA_DIR / "01_row_linkage.csv", index=False)
    pd.DataFrame(clustered.col_linkage, columns=["cluster_a", "cluster_b", "distance", "size"]).to_csv(
        DATA_DIR / "01_col_linkage.csv", index=False)
    (DATA_DIR / "01_row_order.txt").write_text("\n".join(clustered.row_order) + "\n")
    (DATA_DIR / "01_col_order.txt").write_text("\n".join(clustered.col_order) + "\n")

    # Prepare row annotations — drop pigment + n_genes; order genus, clade; abbreviate genus
    GENUS_SHORT = {
        "Prochlorococcus": "PRO",
        "Synechococcus": "SYN",
        "Synechococcus elongatus": "SYN-el",
        "Thermosynechococcus": "THERM",
        "Parasynechococcus": "PARA",
        "Picosynechococcus": "PICO",
    }
    strain_table_indexed = strain_table.set_index("strain")
    annot_cols = ["genus", "clade_canonical"]
    row_annot = strain_table_indexed[annot_cols].copy()
    row_annot["genus"] = row_annot["genus"].map(lambda g: GENUS_SHORT.get(g, g))
    row_annot.index = row_annot.index.map(lambda s: SHORT.get(s, s))
    row_annot = row_annot.rename(columns={"clade_canonical": "clade"})

    # Prepare column labels: prefer consensus_gene_name; fall back to group_id
    col_label_series = groups.set_index("group_id")["consensus_gene_name"].copy()
    # Where consensus_gene_name is NaN, use a short version of group_id
    for gid in col_label_series.index:
        v = col_label_series[gid]
        if pd.isna(v) or v == "":
            col_label_series[gid] = gid.split(":")[-1][:12]
    # Mark singleton groups
    for gid in groups[groups["source"] == "anchor_singleton"]["group_id"]:
        if gid in col_label_series.index:
            col_label_series[gid] = "*" + str(col_label_series[gid])

    # Build pathway category map (group_id → category)
    col_category = groups.set_index("group_id")["consensus_gene_name"].apply(categorize)
    cat_counts = col_category.value_counts()
    logging.info("Pathway category distribution across %d groups:", len(col_category))
    for cat, n in cat_counts.items():
        logging.info("  %-20s %d groups", cat, n)

    # Render
    out_png = FIG_DIR / "01_clustered_heatmap.png"
    out_pdf = FIG_DIR / "01_clustered_heatmap.pdf"
    title = "Cyano N-gene inventory — 19 strains × 61 ortholog groups (A3+B3+C2; * = anchor-singleton group)"
    render_heatmap(
        clustered, row_annot, col_label_series, out_png,
        col_category=col_category, col_category_palette=CATEGORY_PALETTE,
        title=title,
    )
    render_heatmap(
        clustered, row_annot, col_label_series, out_pdf,
        col_category=col_category, col_category_palette=CATEGORY_PALETTE,
        title=title,
    )
    logging.info("Wrote %s and %s", out_png, out_pdf)

    # Brief summary stats on the clustering
    n_total = matrix.size
    n_present = int((matrix > 0).sum().sum())
    n_paralog = int((matrix > 1).sum().sum())
    logging.info("Inventory totals: %d cells, %d present (%.1f%%), %d with paralogs",
                 n_total, n_present, 100*n_present/n_total, n_paralog)


if __name__ == "__main__":
    main()
