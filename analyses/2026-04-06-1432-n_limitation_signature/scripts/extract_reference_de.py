"""Extract DE data from reference nitrogen experiments (Tolonen 2006, Read 2017).

Outputs:
    data/de_reference_tolonen_ndep.csv — all DE rows for Tolonen N-deprivation
    data/de_reference_read_ndep.csv — all DE rows for Read N-depleted
    data/de_reference_tolonen_cyanate.csv — Tolonen cyanate (supplementary)
    data/de_reference_tolonen_urea.csv — Tolonen urea (supplementary)

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_reference_de.py
"""

import sys
from pathlib import Path

# Add analysis dir to path for sig_utils
ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe
from sig_utils.io import save_data_csv

# Experiment IDs
TOLONEN_NDEP = "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray"
TOLONEN_CYANATE = "10.1038/msb4100087_growth_medium_growth_on_cyanate_as_med4_microarray"
TOLONEN_UREA = "10.1038/msb4100087_growth_medium_growth_on_urea_as_med4_microarray"
READ_NDEP = "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq"


def extract_experiment(experiment_id: str, label: str) -> None:
    """Extract all DE rows for an experiment and save to CSV."""
    print(f"Extracting {label}...")

    # First check total size
    result = differential_expression_by_gene(
        experiment_ids=[experiment_id],
        organism="MED4",
        summary=True,
    )
    total = result["total_matching"]
    print(f"  Total rows: {total}")

    # Extract all rows
    result = differential_expression_by_gene(
        experiment_ids=[experiment_id],
        organism="MED4",
        verbose=True,
        limit=None,
    )
    assert not result["truncated"], f"Results truncated for {label}!"

    df = to_dataframe(result)
    print(f"  Extracted {len(df)} rows, {df['locus_tag'].nunique()} genes")
    print(f"  Timepoints: {sorted(df['timepoint'].unique())}")
    print(f"  Status: {df['expression_status'].value_counts().to_dict()}")

    return df


def main():
    # Tolonen N-deprivation (main reference)
    df_tolonen = extract_experiment(TOLONEN_NDEP, "Tolonen 2006 N-deprivation")
    save_data_csv(df_tolonen, "de_reference_tolonen_ndep.csv")

    # Read N-depleted (second reference)
    df_read = extract_experiment(READ_NDEP, "Read 2017 N-depleted")
    save_data_csv(df_read, "de_reference_read_ndep.csv")

    # Supplementary: Tolonen cyanate and urea
    df_cyan = extract_experiment(TOLONEN_CYANATE, "Tolonen 2006 cyanate")
    save_data_csv(df_cyan, "de_reference_tolonen_cyanate.csv")

    df_urea = extract_experiment(TOLONEN_UREA, "Tolonen 2006 urea")
    save_data_csv(df_urea, "de_reference_tolonen_urea.csv")

    print("\nDone. Files saved to data/")


if __name__ == "__main__":
    main()
