"""Build core and extended N-limitation signatures from reference DE data.

Reads: data/de_reference_tolonen_ndep.csv, data/de_reference_read_ndep.csv
Outputs: data/core_signature_genes.csv, data/extended_signature_genes.csv

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/build_signature.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

from sig_utils.io import load_de_csv, save_data_csv
from sig_utils.signature import summarize_de_per_gene, intersect_de_lists


def main():
    # Load reference DE data
    tolonen = load_de_csv("de_reference_tolonen_ndep.csv")
    read = load_de_csv("de_reference_read_ndep.csv")

    print(f"Tolonen: {len(tolonen)} rows, {tolonen['locus_tag'].nunique()} genes")
    print(f"Read: {len(read)} rows, {read['locus_tag'].nunique()} genes")

    # Filter Tolonen to informative timepoints (skip 0h, 3h)
    tolonen_filtered = tolonen[~tolonen["timepoint"].isin(["0h", "3h"])]
    print(f"Tolonen after filtering 0h/3h: {len(tolonen_filtered)} rows")

    # Summarize to one row per gene (peak timepoint)
    tolonen_summary = summarize_de_per_gene(tolonen_filtered)
    read_summary = summarize_de_per_gene(read)

    print(f"\nTolonen significant genes: {len(tolonen_summary)}")
    print(f"  Up: {(tolonen_summary['direction'] == 'up').sum()}")
    print(f"  Down: {(tolonen_summary['direction'] == 'down').sum()}")

    print(f"Read significant genes: {len(read_summary)}")
    print(f"  Up: {(read_summary['direction'] == 'up').sum()}")
    print(f"  Down: {(read_summary['direction'] == 'down').sum()}")

    # Get all locus tags in Read dataset (not just significant) for tagging
    read_all_tags = set(read["locus_tag"].unique())
    print(f"Read total locus tags (including non-significant): {len(read_all_tags)}")

    # Build signatures
    core, extended = intersect_de_lists(
        tolonen_summary, read_summary,
        study_a_name="tolonen", study_b_name="read",
        study_b_all_locus_tags=read_all_tags,
    )

    print(f"\n=== CORE SIGNATURE ===")
    print(f"Total genes: {len(core)}")
    print(f"  Up: {(core['direction'] == 'up').sum()}")
    print(f"  Down: {(core['direction'] == 'down').sum()}")

    if len(core) < 30:
        print("  WARNING: Core signature < 30 genes. Permutation test will be underpowered.")
        print("  Extended signature results should be given more weight.")

    print(f"\n=== EXTENDED SIGNATURE ===")
    print(f"Total genes: {len(extended)}")
    print(f"  By type: {extended['signature_type'].value_counts().to_dict()}")

    # Save
    save_data_csv(core, "core_signature_genes.csv")
    save_data_csv(extended, "extended_signature_genes.csv")

    print("\nDone. Saved to data/core_signature_genes.csv and data/extended_signature_genes.csv")


if __name__ == "__main__":
    main()
