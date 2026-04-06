"""Extract DE data from Weissberg 2025 MED4 experiments.

Outputs:
    data/de_weissberg_med4_rnaseq_axenic.csv
    data/de_weissberg_med4_rnaseq_coculture.csv
    data/de_weissberg_med4_proteomics_axenic.csv
    data/de_weissberg_med4_proteomics_coculture.csv
    data/de_weissberg_med4_all.csv — combined for convenience

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_weissberg_de.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe
from sig_utils.io import save_data_csv

WEISSBERG_DOI = "10.1101/2025.11.24.690089"

EXPERIMENTS = {
    "de_weissberg_med4_rnaseq_axenic.csv": (
        f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic"
    ),
    "de_weissberg_med4_rnaseq_coculture.csv": (
        f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_coculture"
    ),
    "de_weissberg_med4_proteomics_axenic.csv": (
        f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic"
    ),
    "de_weissberg_med4_proteomics_coculture.csv": (
        f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_med4_proteomics_coculture"
    ),
}


def main():
    all_dfs = []

    for filename, experiment_id in EXPERIMENTS.items():
        print(f"Extracting {filename}...")

        result = differential_expression_by_gene(
            experiment_ids=[experiment_id],
            organism="MED4",
            verbose=True,
            limit=None,
        )
        assert not result["truncated"], f"Results truncated for {filename}!"

        df = to_dataframe(result)
        print(f"  {len(df)} rows, {df['locus_tag'].nunique()} genes")
        if "timepoint" in df.columns:
            print(f"  Timepoints: {sorted(df['timepoint'].dropna().unique())}")

        save_data_csv(df, filename)
        df["source_file"] = filename
        all_dfs.append(df)

    # Combined file
    combined = pd.concat(all_dfs, ignore_index=True)
    save_data_csv(combined, "de_weissberg_med4_all.csv")
    print(f"\nCombined: {len(combined)} rows")
    print("Done. Files saved to data/")


if __name__ == "__main__":
    main()
