"""Targeted check: Alteromonas HOT1A3 N-recycling genes in Weissberg coculture.

Step 7 of Approach A: query DE status of N-recycling genes in HOT1A3
coculture time course experiments.

Outputs: data/de_weissberg_hot1a3_nrecycling.csv, results/alteromonas_nrecycling.csv

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/explore_alteromonas.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from multiomics_explorer import genes_by_function, differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe
from sig_utils.io import save_data_csv, save_scores_csv

WEISSBERG_DOI = "10.1101/2025.11.24.690089"

HOT1A3_COCULTURE_EXPERIMENTS = [
    f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_hot1a3_rnaseq_coculture",
    f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_hot1a3_proteomics_coculture",
]

# N-recycling search terms to find relevant HOT1A3 genes
NRECYCLING_SEARCHES = [
    "ammonium AND transport",
    "amino acid AND degradation",
    "peptidase OR protease",
    "deaminase",
    "urease",
    "glutamate AND dehydrogenase",
    "organic nitrogen",
]


def find_nrecycling_genes() -> list[str]:
    """Search for N-recycling genes in HOT1A3."""
    all_tags = set()

    for search in NRECYCLING_SEARCHES:
        result = genes_by_function(
            search_text=search,
            organism="HOT1A3",
            limit=None,
        )
        tags = {r["locus_tag"] for r in result["results"]}
        print(f"  '{search}': {len(tags)} genes")
        all_tags.update(tags)

    print(f"  Total unique N-recycling locus tags: {len(all_tags)}")
    return sorted(all_tags)


def main():
    print("Finding HOT1A3 N-recycling genes...")
    tags = find_nrecycling_genes()

    if not tags:
        print("No N-recycling genes found! Check search terms.")
        return

    print(f"\nExtracting DE for {len(tags)} genes across coculture experiments...")
    all_dfs = []

    for exp_id in HOT1A3_COCULTURE_EXPERIMENTS:
        result = differential_expression_by_gene(
            locus_tags=tags,
            experiment_ids=[exp_id],
            verbose=True,
            limit=None,
        )
        if result["results"]:
            df = to_dataframe(result)
            all_dfs.append(df)
            print(f"  {exp_id.split('_')[-1]}: {len(df)} rows")

    if not all_dfs:
        print("No DE data found for N-recycling genes.")
        return

    combined = pd.concat(all_dfs, ignore_index=True)
    save_data_csv(combined, "de_weissberg_hot1a3_nrecycling.csv")

    # Summary: which genes are upregulated at which timepoints?
    sig = combined[combined["expression_status"].isin(["significant_up", "significant_down"])]
    summary = sig.groupby(["locus_tag", "gene_name", "product", "timepoint", "omics_type"]).agg(
        log2fc=("log2fc", "first"),
        expression_status=("expression_status", "first"),
        rank=("rank", "first"),
        rank_up=("rank_up", "first"),
        rank_down=("rank_down", "first"),
    ).reset_index()

    save_scores_csv(summary, "alteromonas_nrecycling.csv")
    print(f"\nSignificant N-recycling genes: {summary['locus_tag'].nunique()}")
    print(f"Summary saved to results/alteromonas_nrecycling.csv")
    print("\nSignificant genes per timepoint:")
    print(summary.groupby(["timepoint", "expression_status"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
