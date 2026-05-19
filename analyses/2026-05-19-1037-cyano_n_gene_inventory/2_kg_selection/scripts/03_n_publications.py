"""Pull metadata for the 4 N-perturbation publications in the KG.

These publications are Background context for paper.md; their experiments
are not analyzed by this inventory analysis but the field-relevant N
literature represented in the KG should be cited.

Inputs: none (KG).
Outputs:
    data/03_n_publications.csv   one row per publication
    data/03_n_publications.log
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import list_publications

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
LOG_PATH = OUT_DIR / "03_n_publications.log"

N_TREATMENT_DOIS = [
    "10.1038/msb4100087",         # Tolonen 2006 — MED4, MIT9313, microarray
    "10.1038/ismej.2017.88",      # Read 2017 — MED4 N-deprived RNASEQ
    "10.1128/mSystems.00008-17",  # Domínguez-Martín 2017 — SS120 proteomics
    "10.1101/2025.11.24.690089",  # Weissberg 2025 — MED4 + HOT1A3 coculture N-recycling
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    result = list_publications(publication_dois=N_TREATMENT_DOIS, verbose=True, limit=None)
    not_found = result.get("not_found", [])
    if not_found:
        logging.warning("DOIs not found in KG: %s", not_found)

    rows = []
    for pub in result["results"]:
        rows.append({
            "doi": pub["doi"],
            "year": pub["year"],
            "title": pub["title"],
            "journal": pub["journal"],
            "first_author": pub["authors"][0] if pub["authors"] else None,
            "last_author": pub["authors"][-1] if pub["authors"] else None,
            "n_authors": len(pub.get("authors") or []),
            "study_type": pub["study_type"],
            "organisms": " | ".join(pub.get("organisms", [])),
            "experiment_count": pub["experiment_count"],
            "treatment_types": " | ".join(pub.get("treatment_types", [])),
            "background_factors": " | ".join(pub.get("background_factors", [])),
            "omics_types": " | ".join(pub.get("omics_types", [])),
            "growth_phases": " | ".join(pub.get("growth_phases", [])),
        })
    df = pd.DataFrame(rows).sort_values("year").reset_index(drop=True)
    df.to_csv(OUT_DIR / "03_n_publications.csv", index=False)

    logging.info("Pulled %d N-perturbation publications", len(df))
    for _, row in df.iterrows():
        logging.info("  %s — %s %s — %s (%s)",
                     row["year"], row["first_author"], row["doi"],
                     row["organisms"], row["omics_types"])
    logging.info("Wrote %s", OUT_DIR / "03_n_publications.csv")


if __name__ == "__main__":
    main()
