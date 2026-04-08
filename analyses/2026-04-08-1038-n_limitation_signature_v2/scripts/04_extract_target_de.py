"""Step 4: Extract Weissberg 2025 MED4 DE data for scoring.

Uses the same sig_utils/extraction.py utility as Step 2.

Outputs:
    data/de_weissberg_*.csv — one per experiment/platform/condition
    logs/04_extract_target_de.log

Run from multiomics_research root:
    uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/04_extract_target_de.py [--explore]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from sig_utils.extraction import extract_de, extraction_summary, check_marker_genes
from sig_utils.io import save_csv, write_log, load_signature_csv

MARKERS = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]

# Weissberg 2025 target experiments
WEISSBERG_EXPERIMENTS = [
    ("10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic",
     "de_weissberg_rnaseq_axenic.csv", "Weissberg RNA-seq axenic"),
    ("10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_coculture",
     "de_weissberg_rnaseq_coculture.csv", "Weissberg RNA-seq coculture"),
    ("10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic",
     "de_weissberg_proteomics_axenic.csv", "Weissberg proteomics axenic"),
    ("10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_coculture",
     "de_weissberg_proteomics_coculture.csv", "Weissberg proteomics coculture"),
]


def main():
    parser = argparse.ArgumentParser(description="Extract Weissberg 2025 DE data")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces and coverage")
    args = parser.parse_args()

    log_lines = [f"Step 4: Target DE Extraction — {datetime.now().isoformat()}", ""]

    # Load signature for coverage check
    try:
        core_sig = load_signature_csv("core_signature.csv")
        sig_tags = set(core_sig["locus_tag"].unique())
        log_lines.append(f"Core signature loaded: {len(sig_tags)} genes")
    except FileNotFoundError:
        sig_tags = set()
        log_lines.append("WARNING: core_signature.csv not found, skipping coverage check")

    for exp_id, filename, label in WEISSBERG_EXPERIMENTS:
        print(f"\nExtracting {label}...")
        df, envelope = extract_de(experiment_ids=[exp_id])

        summary = extraction_summary(df, envelope)
        print(summary)
        log_lines.append(f"=== {label} ===")
        log_lines.append(summary)

        # Signature coverage
        if sig_tags:
            exp_tags = set(df["locus_tag"].unique())
            overlap = sig_tags & exp_tags
            missing = sig_tags - exp_tags
            log_lines.append(f"Signature coverage: {len(overlap)}/{len(sig_tags)} genes detected")
            if missing:
                log_lines.append(f"Missing signature genes: {len(missing)}")
            print(f"  Signature coverage: {len(overlap)}/{len(sig_tags)} ({len(missing)} missing)")

        # Check for single timepoint issues
        if "timepoint" in df.columns:
            tps = df["timepoint"].unique()
            single_or_null = [str(t) for t in tps if t in ("single", "", None) or str(t) == "nan"]
            if single_or_null:
                warning = f"  WARNING: contains single/null timepoints: {tps.tolist()}"
                print(warning)
                log_lines.append(warning)

        save_csv(df, filename)
        log_lines.append(f"Saved: data/{filename}")

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

        log_lines.append("")

    write_log("\n".join(log_lines), "04_extract_target_de.log")
    print(f"\nDone. See logs/04_extract_target_de.log")


if __name__ == "__main__":
    main()
