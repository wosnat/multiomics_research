"""Run ontology_landscape per organism + search_ontology('nitrogen').

Produces CSVs for each organism's landscape and the nitrogen-ontology search.
Selection of 1-3 ontologies is done interactively by the researcher — this
script only extracts the data. Researcher-approved picks are locked in
ontology_selection.md and exploration/key_pathways.csv.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import ontology_landscape, search_ontology
from multiomics_explorer.analysis import to_dataframe

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


# Ontologies to survey via search_ontology("nitrogen"). Covers all ontologies
# supported by pathway_enrichment except BRITE (which needs a tree argument).
NITROGEN_SEARCH_ONTOLOGIES = [
    "cyanorak_role",
    "tigr_role",
    "go_bp",
    "go_mf",
    "kegg",
    "ec",
    "cog_category",
    "pfam",
]


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step1b.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step1b")
    log.info("Starting 02_ontology_landscape.py")

    # Load classifications.
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    log.info(f"Loaded {len(classified)} classified rows "
             f"({classified['experiment_id'].nunique()} unique experiments)")

    if classified.empty:
        log.error("experiments_classified.csv is empty; cannot produce landscape.")
        return 1

    # Dedupe to one row per experiment — landscape wants experiment_ids, not timepoints.
    unique = classified.drop_duplicates("experiment_id")

    # Per-organism landscape.
    # For MED4: T ∪ R ∪ PC ∪ NC (spec §5 Step 1b do).
    # For non-MED4 organisms: CTX experiments in that organism.
    per_org_rows: list[dict[str, object]] = []
    for org, group in unique.groupby("organism_name"):
        exp_ids = group["experiment_id"].tolist()
        if not exp_ids:
            log.warning(f"Organism {org!r}: no experiment IDs; skipping landscape")
            continue
        log.info(f"Organism {org!r}: {len(exp_ids)} experiments, requesting landscape")
        try:
            result = ontology_landscape(
                organism=org, experiment_ids=exp_ids, limit=100, verbose=True
            )
        except Exception as e:
            log.error(f"ontology_landscape failed for {org!r}: {e}")
            continue

        df = to_dataframe(result)
        if df.empty:
            log.warning(f"Organism {org!r}: landscape returned empty DataFrame")
            continue
        df["organism"] = org
        slug = org.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
        out = ANALYSIS_DIR / "data" / f"landscape_{slug}.csv"
        df.to_csv(out, index=False)
        log.info(f"Wrote {len(df)} landscape rows to {out}")

        # Top 5 by relevance_rank for stdout.
        per_org_rows.append({
            "organism": org,
            "n_exps": len(exp_ids),
            "landscape_rows": len(df),
            "output_path": str(out.relative_to(ANALYSIS_DIR)),
        })

        # Show top of the landscape.
        print(f"\n=== {org} — landscape top 10 by relevance_rank ===")
        top_cols = [
            "ontology_type", "level", "tree", "relevance_rank",
            "n_terms_with_genes", "genome_coverage", "median_genes_per_term",
        ]
        available_cols = [c for c in top_cols if c in df.columns]
        print(
            df.sort_values("relevance_rank").head(10)[available_cols].to_string(index=False)
        )

    # search_ontology('nitrogen') — one call per ontology (no ontology=None mode).
    nitrogen_rows: list[pd.DataFrame] = []
    for ont in NITROGEN_SEARCH_ONTOLOGIES:
        try:
            result = search_ontology(search_text="nitrogen", ontology=ont, limit=50)
        except Exception as e:
            log.error(f"search_ontology failed for {ont!r}: {e}")
            continue
        df = to_dataframe(result)
        if df.empty:
            log.info(f"search_ontology('nitrogen', {ont!r}) returned 0 rows")
            continue
        df["ontology"] = ont
        nitrogen_rows.append(df)
        log.info(f"search_ontology('nitrogen', {ont!r}): {len(df)} rows")

    if nitrogen_rows:
        nitrogen_df = pd.concat(nitrogen_rows, ignore_index=True)
        out = ANALYSIS_DIR / "data" / "nitrogen_ontology_search.csv"
        nitrogen_df.to_csv(out, index=False)
        log.info(f"Wrote {len(nitrogen_df)} nitrogen-term rows to {out}")
        print("\n=== nitrogen search — counts by ontology ===")
        print(nitrogen_df.groupby("ontology").size().to_string())
        print("\n=== nitrogen search — top 15 by score ===")
        show_cols = ["ontology", "id", "name", "score", "level"]
        available = [c for c in show_cols if c in nitrogen_df.columns]
        print(nitrogen_df.sort_values("score", ascending=False).head(15)[available].to_string(index=False))
    else:
        log.warning("All search_ontology calls returned empty; nitrogen_ontology_search.csv not written")

    # Per-organism landscape summary.
    print("\n=== Per-organism landscape outputs ===")
    for r in per_org_rows:
        print(f"  {r['organism']}: {r['n_exps']} exps, {r['landscape_rows']} landscape rows → {r['output_path']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
