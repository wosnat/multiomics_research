"""Step 3: Build pathway definitions for a researcher-chosen ontology and level.

Takes --ontology and --level arguments (decided after reviewing step 2 output).
Builds pathway gene sets and scopes them to each experiment's universe.

Usage:
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/03_define_pathways.py \\
        --ontology go_bp --level 2
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/03_define_pathways.py \\
        --ontology go_bp --level 2 --explore
"""

import argparse
import datetime
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from enrich_utils.io import load_annotations, load_gene_universes, load_hierarchy
from enrich_utils.survey import build_pathway_definitions, scope_pathways_to_universe

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ANALYSIS_DIR / "data"
LOGS_DIR = ANALYSIS_DIR / "logs"
V2_DATA_DIR = Path(__file__).resolve().parents[3] / "2026-04-08-1038-n_limitation_signature_v2" / "data"

MARKER_GENES = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]
MARKER_NAMES = {
    "PMM0920": "glnA",
    "PMM0370": "cynA",
    "PMM0550": "rbcL",
    "PMM1452": "atpD",
    "PMM0030": "PMM0030",
    "PMM0346": "PMM0346",
}


def main(ontology: str, level: int, min_genes: int, explore: bool = False) -> None:
    import pandas as pd

    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    log_path = LOGS_DIR / "03_define_pathways.log"
    log_lines = []

    def log(msg: str) -> None:
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        log_lines.append(line)

    log(f"=== Step 3: Define Pathways (ontology={ontology}, level={level}, min_genes={min_genes}) ===")

    # ------------------------------------------------------------------
    # 1. Load annotations and hierarchy
    # ------------------------------------------------------------------
    ann_path = DATA_DIR / f"annotations_{ontology}.csv"
    if not ann_path.exists():
        raise FileNotFoundError(
            f"Annotations file not found: {ann_path}\n"
            f"Run 01_extract_annotations.py first."
        )
    ann_df = load_annotations(ann_path)
    log(f"Loaded annotations: {len(ann_df)} rows, {ann_df['locus_tag'].nunique()} genes")

    hier_path = DATA_DIR / f"hierarchy_{ontology}.csv"
    if hier_path.exists():
        hier_df = load_hierarchy(hier_path)
        log(f"Loaded hierarchy: {len(hier_df)} edges")
    else:
        hier_df = pd.DataFrame(columns=["child_id", "parent_id", "child_level", "parent_level"])
        log(f"No hierarchy file found — treating as flat ontology")

    # ------------------------------------------------------------------
    # 2. Build term_names dict from annotations
    # ------------------------------------------------------------------
    if "term_name" in ann_df.columns:
        term_names = dict(zip(ann_df["term_id"], ann_df["term_name"]))
        log(f"Built term_names dict: {len(term_names)} entries")
    else:
        term_names = None
        log("No term_name column — pathway names will default to term IDs")

    # ------------------------------------------------------------------
    # 3. Build pathway definitions
    # ------------------------------------------------------------------
    log(f"Building pathway definitions at level={level}, min_genes={min_genes}...")
    pathway_defs = build_pathway_definitions(
        ann_df, hier_df, level=level, min_genes=min_genes, term_names=term_names
    )
    log(f"  {len(pathway_defs)} pathways (>= {min_genes} genes)")

    if len(pathway_defs) == 0:
        log("WARNING: No pathways defined. Check ontology/level/min_genes arguments.")

    # Gene count distribution
    if len(pathway_defs) > 0:
        sizes = pathway_defs["gene_count"]
        log(
            f"  gene count distribution: min={sizes.min()}, "
            f"q25={sizes.quantile(0.25):.0f}, "
            f"median={sizes.median():.0f}, "
            f"q75={sizes.quantile(0.75):.0f}, "
            f"max={sizes.max()}"
        )

    if explore and len(pathway_defs) > 0:
        log("--- Marker gene pathway membership ---")
        for locus_tag in MARKER_GENES:
            gene_name = MARKER_NAMES.get(locus_tag, locus_tag)
            in_pathways = pathway_defs[
                pathway_defs["locus_tags"].apply(lambda tags: locus_tag in tags)
            ]
            if len(in_pathways) == 0:
                log(f"  {locus_tag} ({gene_name}): not in any pathway")
            else:
                for _, pw in in_pathways.iterrows():
                    log(f"  {locus_tag} ({gene_name}): {pw['pathway_id']} | {pw['pathway_name']} ({pw['gene_count']} genes)")

        log("--- Size distribution histogram ---")
        bins = [0, 5, 10, 20, 50, 100, 200, 500, 9999]
        bin_labels = ["<5", "5-9", "10-19", "20-49", "50-99", "100-199", "200-499", "500+"]
        for i, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
            count = ((pathway_defs["gene_count"] >= lo) & (pathway_defs["gene_count"] < hi)).sum()
            if i == 0:
                count = (pathway_defs["gene_count"] < hi).sum()
            log(f"  {bin_labels[i]:10s}: {count} pathways")

    # ------------------------------------------------------------------
    # 4. Scope pathways to each experiment universe
    # ------------------------------------------------------------------
    log("Loading gene universes and scoping pathways...")
    universes = load_gene_universes(V2_DATA_DIR)

    coverage_rows = []
    for exp_name, universe in sorted(universes.items()):
        scoped = scope_pathways_to_universe(pathway_defs, universe)
        # Summary for this experiment
        mean_cov = scoped["coverage"].mean() if len(scoped) > 0 else 0.0
        n_represented = (scoped["n_in_universe"] > 0).sum() if len(scoped) > 0 else 0
        log(
            f"  {exp_name}: {n_represented}/{len(pathway_defs)} pathways represented, "
            f"mean_coverage={mean_cov:.3f}"
        )
        coverage_rows.append({
            "experiment": exp_name,
            "n_pathways_total": len(pathway_defs),
            "n_pathways_represented": n_represented,
            "mean_coverage": mean_cov,
        })

    coverage_df = pd.DataFrame(coverage_rows)

    # ------------------------------------------------------------------
    # 5. Save outputs
    # ------------------------------------------------------------------
    # Pathway definitions — convert sets to pipe-delimited strings for CSV
    defs_out = pathway_defs.copy()
    defs_out["locus_tags"] = defs_out["locus_tags"].apply(
        lambda tags: "|".join(sorted(tags)) if isinstance(tags, set) else tags
    )
    defs_path = DATA_DIR / "pathway_definitions.csv"
    defs_out.to_csv(defs_path, index=False)
    log(f"Saved pathway_definitions.csv ({len(defs_out)} pathways)")

    cov_path = DATA_DIR / "pathway_coverage_per_experiment.csv"
    coverage_df.to_csv(cov_path, index=False)
    log(f"Saved pathway_coverage_per_experiment.csv")

    # Write log file
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Define pathway gene sets at a chosen ontology level"
    )
    parser.add_argument(
        "--ontology",
        required=True,
        help="Ontology name (e.g. go_bp, kegg, cyanorak_role). "
             "Must match an annotations_<ontology>.csv file in data/.",
    )
    parser.add_argument(
        "--level",
        type=int,
        required=True,
        help="Ontology level to roll up to (e.g. 2). "
             "Choose based on ontology_ranking.csv and ontology_profiles.csv from step 2.",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=5,
        dest="min_genes",
        help="Minimum genes per pathway (default: 5)",
    )
    parser.add_argument(
        "--explore",
        action="store_true",
        help="Trace marker genes through pathways and show size distribution",
    )
    args = parser.parse_args()
    main(ontology=args.ontology, level=args.level, min_genes=args.min_genes, explore=args.explore)
