"""Build the enriched 19-strain cyano table:
 - KG facts from list_organisms (genus, species, gene_count, experiment_count, KG clade field)
 - Cyanorak source CSV enrichment for clade/subcluster/pigment (Prochlorococcus + marine Synechococcus)

KG `clade` values are known to be wrong for Synechococcus (null) and
NATL1A/NATL2A (LLII instead of LLI) — see gaps_and_friction.md. The
Cyanorak source CSV values are the authoritative reference.

Inputs:
    /home/oweissberg/work/multiomics/multiomics_biocypher_kg/data/Cyanorak  Organism Table  prochlorococcus.csv
    /home/oweissberg/work/multiomics/multiomics_biocypher_kg/data/Cyanorak  Organism Table  synechococcus.csv
    (queries list_organisms for the KG side)

Outputs:
    data/02_strain_table.csv  one row per cyano strain
    data/02_strain_table.log
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import list_organisms

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
LOG_PATH = OUT_DIR / "02_strain_table.log"

CYANORAK_DATA_DIR = Path("/home/oweissberg/work/multiomics/multiomics_biocypher_kg/data")
CYANORAK_PRO_CSV = CYANORAK_DATA_DIR / "Cyanorak  Organism Table  prochlorococcus.csv"
CYANORAK_SYN_CSV = CYANORAK_DATA_DIR / "Cyanorak  Organism Table  synechococcus.csv"

# Strain order — Prochlorococcus first by clade then alphabetical, Synechococcus marine, then non-Cyanorak.
SCOPE_STRAINS = [
    "Prochlorococcus AS9601",
    "Prochlorococcus MED4",
    "Prochlorococcus MIT0801",
    "Prochlorococcus MIT9301",
    "Prochlorococcus MIT9303",
    "Prochlorococcus MIT9312",
    "Prochlorococcus MIT9313",
    "Prochlorococcus NATL1A",
    "Prochlorococcus NATL2A",
    "Prochlorococcus RSP50",
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)",
    "Synechococcus WH8102",
    "Synechococcus WH7803",
    "Synechococcus CC9311",
    "Synechococcus sp. BL107",
    "Synechococcus PCC 7002",
    "Synechococcus elongatus PCC 7942",
    "Synechococcus elongatus UTEX 2973",
    "Thermosynechococcus vestitus BP-1",
]

# Map KG organism_name → Cyanorak CSV "Name" column.
KG_TO_CYANORAK_NAME = {
    "Prochlorococcus AS9601": "Pro_AS9601",
    "Prochlorococcus MED4": "Pro_MED4",
    "Prochlorococcus MIT0801": "Pro_MIT0801",
    "Prochlorococcus MIT9301": "Pro_MIT9301",
    "Prochlorococcus MIT9303": "Pro_MIT9303",
    "Prochlorococcus MIT9312": "Pro_MIT9312",
    "Prochlorococcus MIT9313": "Pro_MIT9313",
    "Prochlorococcus NATL1A": "Pro_NATL1A",
    "Prochlorococcus NATL2A": "Pro_NATL2A",
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)": "Pro_SS120",
    # RSP50 — not in Cyanorak Pro CSV
    "Synechococcus WH8102": "Syn_WH8102",
    "Synechococcus WH7803": "Syn_WH7803",
    "Synechococcus CC9311": "Syn_CC9311",
    "Synechococcus sp. BL107": "Syn_BL107",
    # PCC 7002, PCC 7942, UTEX 2973, T. vestitus — not in Cyanorak picocyano CSVs
}


def load_cyanorak_csv(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]
    return df


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    # KG side — pull all organisms (need to paginate; 37 total).
    kg_pages = []
    offset = 0
    while True:
        r = list_organisms(limit=50, offset=offset)
        kg_pages.extend(r["results"])
        if not r["truncated"]:
            break
        offset += 50
    kg_df = pd.DataFrame(kg_pages)
    logging.info("list_organisms returned %d organism rows total", len(kg_df))

    # Filter to our 19 cyano scope.
    scope_df = kg_df[kg_df["organism_name"].isin(SCOPE_STRAINS)].copy()
    missing = set(SCOPE_STRAINS) - set(scope_df["organism_name"])
    if missing:
        raise RuntimeError(f"Scope strains not found in list_organisms: {missing}")
    logging.info("Filtered to %d scope strains", len(scope_df))

    # Cyanorak source CSV side.
    pro = load_cyanorak_csv(CYANORAK_PRO_CSV)
    syn = load_cyanorak_csv(CYANORAK_SYN_CSV)
    cyanorak_combined = pd.concat([pro.assign(_source="pro"), syn.assign(_source="syn")], ignore_index=True)
    logging.info("Cyanorak Pro CSV: %d rows; Syn CSV: %d rows", len(pro), len(syn))

    # Build output.
    rows = []
    for strain in SCOPE_STRAINS:
        kg = scope_df[scope_df["organism_name"] == strain].iloc[0]
        cyanorak_name = KG_TO_CYANORAK_NAME.get(strain)
        if cyanorak_name:
            cy_match = cyanorak_combined[cyanorak_combined["Name"] == cyanorak_name]
            if cy_match.empty:
                logging.warning("KG_TO_CYANORAK_NAME mapping has '%s' but not in CSV", cyanorak_name)
                cy = None
            else:
                cy = cy_match.iloc[0]
        else:
            cy = None

        rows.append({
            "strain": strain,
            "genus": kg["genus"],
            "species_kg": kg["species"],
            "strain_short": kg["strain"],
            "ncbi_taxon_id": kg["ncbi_taxon_id"],
            "gene_count_kg": kg["gene_count"],
            "experiment_count_kg": kg["experiment_count"],
            "publication_count_kg": kg["publication_count"],
            "clade_kg": kg["clade"],
            # Cyanorak source CSV columns (NaN where strain not in CSV)
            "cyanorak_name": cyanorak_name if cy is not None else None,
            "cyanorak_species": cy["Species"] if cy is not None else None,
            "cyanorak_subcluster": cy["SubCluster"] if cy is not None else None,
            "cyanorak_clade": cy["Clade"] if cy is not None else None,
            "cyanorak_subclade": cy["SubClade"] if cy is not None else None,
            "cyanorak_pigment": cy["Pigment type"] if cy is not None else None,
            "in_cyanorak_picocyano_csv": cy is not None,
        })

    out = pd.DataFrame(rows)
    # Canonical clade = Cyanorak CSV value if available, else KG clade.
    out["clade_canonical"] = out["cyanorak_clade"].fillna(out["clade_kg"])

    out.to_csv(OUT_DIR / "02_strain_table.csv", index=False)

    # Report mismatches between KG and Cyanorak.
    mismatched = out[
        out["cyanorak_clade"].notna()
        & out["clade_kg"].notna()
        & (out["cyanorak_clade"] != out["clade_kg"])
    ]
    logging.info("---")
    logging.info("Strains with KG/Cyanorak clade mismatch: %d", len(mismatched))
    for _, row in mismatched.iterrows():
        logging.info(
            "  %s — KG=%s, Cyanorak=%s",
            row["strain"], row["clade_kg"], row["cyanorak_clade"],
        )
    logging.info("Strains with null KG clade (where Cyanorak has it): %d",
                 ((out["cyanorak_clade"].notna()) & (out["clade_kg"].isna())).sum())
    logging.info("Strains outside Cyanorak picocyano CSVs: %d",
                 (~out["in_cyanorak_picocyano_csv"]).sum())
    for _, row in out[~out["in_cyanorak_picocyano_csv"]].iterrows():
        logging.info("  %s (KG clade: %s)", row["strain"], row["clade_kg"])
    logging.info("Wrote %s", OUT_DIR / "02_strain_table.csv")


if __name__ == "__main__":
    main()
