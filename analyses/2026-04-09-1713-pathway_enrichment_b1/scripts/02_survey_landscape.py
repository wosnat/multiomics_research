"""Step 2: Profile each ontology and rank by coverage + term-size sweet spot.

Reads step 1 outputs, runs survey_ontology for each ontology, ranks them,
and saves profiles and ranking.

Usage:
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/02_survey_landscape.py
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/02_survey_landscape.py --explore
"""

import argparse
import datetime
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from enrich_utils.io import load_annotations, load_gene_universes, load_hierarchy
from enrich_utils.survey import rank_ontologies, survey_ontology

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ANALYSIS_DIR / "data"
RESULTS_DIR = ANALYSIS_DIR / "results"
LOGS_DIR = ANALYSIS_DIR / "logs"
V2_DATA_DIR = ANALYSIS_DIR.parent / "2026-04-08-1038-n_limitation_signature_v2" / "data"

MARKER_GENES = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]


def main(explore: bool = False) -> None:
    import pandas as pd

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    log_path = LOGS_DIR / "02_survey_landscape.log"
    log_lines = []

    def log(msg: str) -> None:
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        log_lines.append(line)

    log("=== Step 2: Survey Annotation Landscape ===")

    # ------------------------------------------------------------------
    # 1. Load gene universes
    # ------------------------------------------------------------------
    log(f"Loading gene universes from {V2_DATA_DIR}")
    universes = load_gene_universes(V2_DATA_DIR)
    log(f"  Loaded {len(universes)} experiment universes")

    # Use union of all universes as the combined background
    combined_universe = set()
    for genes in universes.values():
        combined_universe |= genes
    log(f"  Combined background universe: {len(combined_universe)} genes")

    # ------------------------------------------------------------------
    # 2. Discover available annotation files
    # ------------------------------------------------------------------
    annotation_files = sorted(DATA_DIR.glob("annotations_*.csv"))
    if not annotation_files:
        raise FileNotFoundError(
            f"No annotation files found in {DATA_DIR}. Run 01_extract_annotations.py first."
        )
    log(f"Found {len(annotation_files)} annotation files")

    # ------------------------------------------------------------------
    # 3. Survey each ontology
    # ------------------------------------------------------------------
    profiles = {}

    for ann_path in annotation_files:
        ontology = ann_path.stem.replace("annotations_", "")
        log(f"  Surveying {ontology}...")

        ann_df = load_annotations(ann_path)
        if len(ann_df) == 0:
            log(f"    (empty annotations — skipping)")
            continue

        hier_path = DATA_DIR / f"hierarchy_{ontology}.csv"
        if hier_path.exists():
            hier_df = load_hierarchy(hier_path)
        else:
            # Flat ontology — empty hierarchy
            hier_df = pd.DataFrame(columns=["child_id", "parent_id", "child_level", "parent_level"])

        profile = survey_ontology(ann_df, hier_df, combined_universe)
        profiles[ontology] = profile

        # Log per-ontology stats
        log(
            f"    genome_coverage={profile['genome_coverage']:.3f}  "
            f"annotated={profile['n_annotated']}/{profile['n_universe']}"
        )
        for lvl in profile.get("per_level", []):
            log(
                f"    level={lvl['level']}  "
                f"genes={lvl.get('n_genes_at_level', '?')}/{profile['n_universe']} ({lvl.get('genome_coverage', 0):.0%})  "
                f"n_terms={lvl['n_terms_with_genes']}  "
                f"median={lvl.get('median_genes')}  "
                f"[{lvl.get('min_genes')}, {lvl.get('max_genes')}]"
            )

        if explore:
            log(f"    --- Full profile for {ontology} ---")
            log(f"    {profile}")

    # ------------------------------------------------------------------
    # 4. Rank ontologies
    # ------------------------------------------------------------------
    log("Ranking ontologies...")
    ranking_df = rank_ontologies(profiles)
    log(f"  Ranked {len(ranking_df)} ontologies")
    for _, row in ranking_df.iterrows():
        lvl = row.get('best_level', None)
        lvl_str = f"level {int(lvl)}" if pd.notna(lvl) else "none"
        log(
            f"    rank={int(row['rank'])}  {row['ontology']:20s}  "
            f"genome_cov={row['genome_coverage']:.3f}  "
            f"best_level={lvl_str}  "
            f"level_cov={row.get('best_level_genome_coverage', 0):.0%}  "
            f"median={row.get('best_level_median_genes', 0):.0f}  "
            f"max={int(row.get('best_level_max_genes', 0))}  "
            f"n_terms={int(row.get('best_level_n_terms', 0))}  "
            f"score={row['score']:.3f}"
        )

    # ------------------------------------------------------------------
    # 5. Save outputs
    # ------------------------------------------------------------------
    # Flatten profiles to DataFrame
    profile_rows = []
    for ontology, profile in profiles.items():
        base = {
            "ontology": ontology,
            "genome_coverage": profile["genome_coverage"],
            "n_annotated": profile["n_annotated"],
            "n_universe": profile["n_universe"],
        }
        if not profile.get("per_level"):
            base["level"] = None
            base["n_terms_with_genes"] = None
            base["min_genes"] = None
            base["q1_genes"] = None
            base["median_genes"] = None
            base["q3_genes"] = None
            base["max_genes"] = None
            profile_rows.append(base)
        else:
            for lvl in profile["per_level"]:
                row = dict(base)
                row.update(lvl)
                profile_rows.append(row)

    profiles_df = pd.DataFrame(profile_rows)
    profiles_path = DATA_DIR / "ontology_profiles.csv"
    profiles_df.to_csv(profiles_path, index=False)
    log(f"Saved ontology_profiles.csv ({len(profiles_df)} rows)")

    ranking_path = RESULTS_DIR / "ontology_ranking.csv"
    ranking_df.to_csv(ranking_path, index=False)
    log(f"Saved ontology_ranking.csv ({len(ranking_df)} rows)")

    # Write log file
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Survey and rank annotation ontologies")
    parser.add_argument(
        "--explore",
        action="store_true",
        help="Print full profile comparison for each ontology",
    )
    args = parser.parse_args()
    main(explore=args.explore)
