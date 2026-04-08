"""Step 2: Extract DE data for reference and control experiments.

Uses sig_utils/extraction.py for consistent schema.
Experiment IDs from Step 1 classification.

Outputs:
    data/de_*.csv — one per experiment
    logs/02_extract_reference_de.log

Run from multiomics_research root:
    uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/02_extract_reference_de.py [--explore]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

from sig_utils.extraction import extract_de, extraction_summary, check_marker_genes
from sig_utils.io import save_csv, write_log

# 6 marker genes (including PMM0346 for proteomics coverage)
MARKERS = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]

# Experiment registry from Step 1 classification
# (experiment_id, output_filename, label, role)
EXPERIMENTS = [
    # References
    ("10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray",
     "de_ref_tolonen_ndep.csv", "Tolonen 2006 N-deprivation", "reference"),
    ("10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq",
     "de_ref_read_ndep.csv", "Read 2017 N-depleted", "reference"),

    # Negative controls
    ("10.1038/msb4100087_growth_medium_growth_on_cyanate_as_med4_microarray",
     "de_ctrl_tolonen_cyanate.csv", "Tolonen 2006 cyanate", "negative_control"),
    ("10.1038/msb4100087_growth_medium_growth_on_urea_as_med4_microarray",
     "de_ctrl_tolonen_urea.csv", "Tolonen 2006 urea", "negative_control"),
    ("10.1038/ismej.2016.70_coculture_alteromonas_hot1a3_med4_rnaseq",
     "de_ctrl_aharonovich_coculture.csv", "Aharonovich 2016 coculture", "negative_control"),
    ("10.1128/JB.01097-06_light_stress_high_white_light_55_med4_microarray",
     "de_ctrl_steglich_high_white_light.csv", "Steglich 2006 high white light", "negative_control"),
]


def main():
    parser = argparse.ArgumentParser(description="Extract reference and control DE data")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces")
    args = parser.parse_args()

    log_lines = [f"Step 2: Reference DE Extraction — {datetime.now().isoformat()}", ""]

    for exp_id, filename, label, role in EXPERIMENTS:
        print(f"\nExtracting {label} ({role})...")
        df, envelope = extract_de(experiment_ids=[exp_id])

        summary = extraction_summary(df, envelope)
        print(summary)
        log_lines.append(f"=== {label} ({role}) ===")
        log_lines.append(summary)

        # Check for metadata issues
        if "timepoint" in df.columns:
            tps = df["timepoint"].unique()
            single_or_null = [t for t in tps if t in ("single", "", None) or str(t) == "nan"]
            if single_or_null:
                warning = f"  WARNING: contains single/null timepoints: {tps.tolist()}"
                print(warning)
                log_lines.append(warning)

        # Save
        save_csv(df, filename)
        log_lines.append(f"Saved: data/{filename}")

        # Marker gene exploration
        if args.explore:
            marker_df = check_marker_genes(df, MARKERS)
            if len(marker_df) > 0:
                print(f"\n  Marker genes in {label}:")
                cols = ["locus_tag", "gene_name", "timepoint", "log2fc",
                        "expression_status", "rank_up", "rank_down"]
                cols = [c for c in cols if c in marker_df.columns]
                print(marker_df[cols].to_string(index=False))
                log_lines.append(f"\nMarker genes:")
                log_lines.append(marker_df[cols].to_string(index=False))
            else:
                print(f"  No marker genes found in {label}")
                log_lines.append("No marker genes found")

        log_lines.append("")

    write_log("\n".join(log_lines), "02_extract_reference_de.log")
    print(f"\nDone. See logs/02_extract_reference_de.log")


if __name__ == "__main__":
    main()
