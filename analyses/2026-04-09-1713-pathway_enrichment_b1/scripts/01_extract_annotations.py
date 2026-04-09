"""Step 1: Extract annotations for the full MED4 genome across all ontologies.

Extracts gene × term leaf annotations and hierarchies from the KG.
Saves per-ontology CSVs and logs coverage statistics.

Usage:
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/01_extract_annotations.py
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/01_extract_annotations.py --explore
"""

import argparse
import datetime
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from enrich_utils.extraction import ONTOLOGY_EDGES, extract_annotations, extract_hierarchy
from enrich_utils.io import load_gene_universes

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ANALYSIS_DIR / "data"
LOGS_DIR = ANALYSIS_DIR / "logs"
V2_DATA_DIR = ANALYSIS_DIR.parent / "2026-04-08-1038-n_limitation_signature_v2" / "data"

MARKER_GENES = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]
MARKER_NAMES = {
    "PMM0920": "glnA",
    "PMM0370": "cynA",
    "PMM0550": "rbcL",
    "PMM1452": "atpD",
    "PMM0030": "PMM0030",
    "PMM0346": "PMM0346",
}

# Ontologies that have hierarchies (non-empty hierarchy_rels)
HIERARCHICAL_ONTOLOGIES = {
    ont for ont, cfg in ONTOLOGY_EDGES.items() if cfg["hierarchy_rels"]
}


def get_full_med4_gene_list() -> list[str]:
    """Get all MED4 locus tags via genes_by_function wildcard query."""
    from multiomics_explorer import genes_by_function
    from multiomics_explorer.analysis import to_dataframe

    result = genes_by_function("*", organism="MED4", limit=None)
    if not result.get("results"):
        raise RuntimeError("genes_by_function returned no results for MED4")
    df = to_dataframe(result)
    if "locus_tag" not in df.columns:
        raise RuntimeError(f"Expected 'locus_tag' column, got: {list(df.columns)}")
    return sorted(df["locus_tag"].dropna().unique().tolist())


def main(explore: bool = False) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    log_path = LOGS_DIR / "01_extract_annotations.log"
    log_lines = []

    def log(msg: str) -> None:
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        log_lines.append(line)

    log("=== Step 1: Extract Annotations ===")

    # ------------------------------------------------------------------
    # 1. Load gene universes from v2 DE files
    # ------------------------------------------------------------------
    log(f"Loading gene universes from {V2_DATA_DIR}")
    universes = load_gene_universes(V2_DATA_DIR)
    log(f"  Loaded {len(universes)} experiment universes")
    for name, genes in universes.items():
        log(f"    {name}: {len(genes)} genes")

    # Save combined gene universe summary
    import pandas as pd
    rows = [{"experiment": name, "n_genes": len(genes)} for name, genes in universes.items()]
    pd.DataFrame(rows).to_csv(DATA_DIR / "gene_universes.csv", index=False)
    log(f"  Saved gene_universes.csv")

    # ------------------------------------------------------------------
    # 2. Get full MED4 gene list
    # ------------------------------------------------------------------
    log("Fetching full MED4 gene list via genes_by_function('*', organism='MED4')")
    all_locus_tags = get_full_med4_gene_list()
    log(f"  Full MED4 genome: {len(all_locus_tags)} genes")

    # ------------------------------------------------------------------
    # 3. Extract annotations for each ontology
    # ------------------------------------------------------------------
    log("Extracting annotations for all ontologies...")
    annotation_stats = {}

    for ontology in sorted(ONTOLOGY_EDGES):
        log(f"  Extracting {ontology}...")
        ann_df = extract_annotations(all_locus_tags, ontology)
        out_path = DATA_DIR / f"annotations_{ontology}.csv"
        ann_df.to_csv(out_path, index=False)

        n_genes = ann_df["locus_tag"].nunique() if len(ann_df) > 0 else 0
        n_terms = ann_df["term_id"].nunique() if len(ann_df) > 0 else 0
        n_unannotated = len(all_locus_tags) - n_genes
        annotation_stats[ontology] = {
            "n_genes": n_genes,
            "n_terms": n_terms,
            "n_unannotated": n_unannotated,
        }
        log(f"    genes={n_genes}, terms={n_terms}, unannotated={n_unannotated}")

        if explore:
            log(f"    --- Marker gene traces for {ontology} ---")
            if len(ann_df) > 0:
                for locus_tag in MARKER_GENES:
                    gene_name = MARKER_NAMES.get(locus_tag, locus_tag)
                    rows_for_gene = ann_df[ann_df["locus_tag"] == locus_tag]
                    if len(rows_for_gene) == 0:
                        log(f"      {locus_tag} ({gene_name}): no annotations")
                    else:
                        terms = rows_for_gene[["term_id", "term_name"]].values.tolist()
                        for term_id, term_name in terms:
                            log(f"      {locus_tag} ({gene_name}): {term_id} | {term_name}")
            else:
                log(f"      (no annotations extracted)")

    # ------------------------------------------------------------------
    # 4. Extract hierarchies for hierarchical ontologies
    # ------------------------------------------------------------------
    log("Extracting hierarchies for hierarchical ontologies...")

    for ontology in sorted(HIERARCHICAL_ONTOLOGIES):
        log(f"  Extracting hierarchy for {ontology}...")
        hier_df = extract_hierarchy(ontology)
        out_path = DATA_DIR / f"hierarchy_{ontology}.csv"
        hier_df.to_csv(out_path, index=False)

        n_edges = len(hier_df)
        if n_edges > 0:
            levels = sorted(hier_df["child_level"].unique().tolist())
            log(f"    edges={n_edges}, levels={levels}")
        else:
            log(f"    (no hierarchy edges found)")

    # ------------------------------------------------------------------
    # 5. Summary log
    # ------------------------------------------------------------------
    log("")
    log("=== Annotation Summary ===")
    for ontology, stats in sorted(annotation_stats.items()):
        log(
            f"  {ontology:20s}  genes={stats['n_genes']:4d}  "
            f"terms={stats['n_terms']:4d}  unannotated={stats['n_unannotated']:4d}"
        )

    # Write log file
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract MED4 annotations from KG")
    parser.add_argument(
        "--explore",
        action="store_true",
        help="Print marker gene traces through each ontology",
    )
    args = parser.parse_args()
    main(explore=args.explore)
