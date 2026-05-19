"""Pull the Cyanorak E.4 ∪ D.1.3 anchor gene set across the 13 Cyanorak-annotated cyano strains in the KG.

Inputs: none (queries KG directly).
Outputs:
    data/01_anchor_genes.csv          one row per (strain, locus_tag, term_id)
    data/01_anchor_summary.csv        one row per strain with anchor gene count
    data/01_anchor_genes.log          run log
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import genes_by_ontology
from multiomics_explorer.analysis import to_dataframe

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
LOG_PATH = OUT_DIR / "01_anchor_genes.log"

CYANORAK_ANNOTATED_STRAINS = [
    # exact list_organisms organism_name values for the 13 Cyanorak-annotated cyano in KG
    "Prochlorococcus AS9601",
    "Prochlorococcus MED4",
    "Prochlorococcus MIT9301",
    "Prochlorococcus MIT9303",
    "Prochlorococcus MIT9312",
    "Prochlorococcus MIT9313",
    "Prochlorococcus NATL1A",
    "Prochlorococcus NATL2A",
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)",
    "Synechococcus WH8102",
    "Synechococcus WH7803",
    "Synechococcus CC9311",
    "Synechococcus sp. BL107",
]

N_ROLE_TERMS = [
    "cyanorak.role:E.4",   # Central intermediary metabolism > Nitrogen metabolism
    "cyanorak.role:D.1.3", # Cellular processes > Adaptation > Nitrogen
]


def fetch_anchor_genes(strain: str) -> pd.DataFrame:
    result = genes_by_ontology(
        organism=strain,
        ontology="cyanorak_role",
        term_ids=N_ROLE_TERMS,
        min_gene_set_size=0,
        max_gene_set_size=1000,
        limit=None,
    )
    if result["total_matching"] == 0:
        return pd.DataFrame(columns=["strain", "organism_name_kg", "locus_tag", "gene_name", "product", "gene_category", "term_id"])
    df = to_dataframe(result)
    df.insert(0, "strain", strain)
    df.insert(1, "organism_name_kg", result["organism_name"])
    return df[["strain", "organism_name_kg", "locus_tag", "gene_name", "product", "gene_category", "term_id"]]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    per_strain_frames: list[pd.DataFrame] = []
    summary_rows: list[dict] = []

    for strain in CYANORAK_ANNOTATED_STRAINS:
        df = fetch_anchor_genes(strain)
        per_strain_frames.append(df)
        unique_loci = df["locus_tag"].nunique()
        unique_genes_named = df.dropna(subset=["gene_name"])["gene_name"].nunique()
        logging.info(
            "%s: %d rows, %d unique locus_tags, %d unique gene names",
            strain, len(df), unique_loci, unique_genes_named,
        )
        summary_rows.append({
            "strain": strain,
            "anchor_rows": len(df),
            "unique_locus_tags": unique_loci,
            "unique_gene_names": unique_genes_named,
        })

    anchor = pd.concat(per_strain_frames, ignore_index=True)
    anchor.to_csv(OUT_DIR / "01_anchor_genes.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "01_anchor_summary.csv", index=False)

    total_unique_loci = anchor["locus_tag"].nunique()
    logging.info("---")
    logging.info("Total anchor rows across 13 strains: %d", len(anchor))
    logging.info("Total unique locus_tags across 13 strains: %d", total_unique_loci)
    logging.info("Wrote %s and %s", OUT_DIR / "01_anchor_genes.csv", OUT_DIR / "01_anchor_summary.csv")


if __name__ == "__main__":
    main()
