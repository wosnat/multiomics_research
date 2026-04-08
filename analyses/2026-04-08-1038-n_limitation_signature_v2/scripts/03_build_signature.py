"""Step 3: Build N-limitation gene signature from reference studies.

Uses sig_utils/signature.py functions on real DE data from Step 2.

Outputs:
    data/core_signature.csv
    data/extended_signature.csv
    logs/03_build_signature.log

Run from multiomics_research root:
    uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/03_build_signature.py [--explore]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from sig_utils.signature import summarize_per_gene, intersect_references
from sig_utils.extraction import check_marker_genes
from sig_utils.io import load_de_csv, save_csv, write_log

MARKERS = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]

STUDY_A_FILE = "de_ref_tolonen_ndep.csv"
STUDY_B_FILE = "de_ref_read_ndep.csv"
STUDY_A_NAME = "tolonen"
STUDY_B_NAME = "read"
# Exclude early timepoints (transient responses, not sustained N-deprivation)
STUDY_A_EXCLUDE_TIMEPOINTS = ["0h", "3h"]  # Tolonen: near-zero DE at 0h/3h
STUDY_B_EXCLUDE_TIMEPOINTS = ["3h"]         # Read: transient up-bias at 3h reverses by 12h


def main():
    parser = argparse.ArgumentParser(description="Build N-limitation signature")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces")
    args = parser.parse_args()

    log_lines = [f"Step 3: Build Signature — {datetime.now().isoformat()}", ""]

    # Load reference DE data
    print("Loading reference DE data...")
    de_a = load_de_csv(STUDY_A_FILE)
    de_b = load_de_csv(STUDY_B_FILE)

    log_lines.append(f"Study A ({STUDY_A_NAME}): {len(de_a)} rows, {de_a['locus_tag'].nunique()} genes")
    log_lines.append(f"Study B ({STUDY_B_NAME}): {len(de_b)} rows, {de_b['locus_tag'].nunique()} genes")

    # Exclude early timepoints (transient responses)
    before_a = len(de_a)
    de_a = de_a[~de_a["timepoint"].isin(STUDY_A_EXCLUDE_TIMEPOINTS)]
    log_lines.append(f"Excluded {STUDY_A_NAME} timepoints {STUDY_A_EXCLUDE_TIMEPOINTS}: {before_a} → {len(de_a)} rows")
    print(f"Excluded Tolonen {STUDY_A_EXCLUDE_TIMEPOINTS}: {before_a} → {len(de_a)} rows")

    before_b = len(de_b)
    de_b = de_b[~de_b["timepoint"].isin(STUDY_B_EXCLUDE_TIMEPOINTS)]
    log_lines.append(f"Excluded {STUDY_B_NAME} timepoints {STUDY_B_EXCLUDE_TIMEPOINTS}: {before_b} → {len(de_b)} rows")
    print(f"Excluded Read {STUDY_B_EXCLUDE_TIMEPOINTS}: {before_b} → {len(de_b)} rows")

    # Filter funnel
    print("\n=== Filter funnel ===")

    # Summarize per gene
    summary_a = summarize_per_gene(de_a)
    summary_b = summarize_per_gene(de_b)

    a_up = (summary_a["direction"] == "up").sum()
    a_down = (summary_a["direction"] == "down").sum()
    b_up = (summary_b["direction"] == "up").sum()
    b_down = (summary_b["direction"] == "down").sum()

    print(f"{STUDY_A_NAME}: {len(summary_a)} significant genes (up={a_up}, down={a_down})")
    print(f"{STUDY_B_NAME}: {len(summary_b)} significant genes (up={b_up}, down={b_down})")
    log_lines.append(f"\n--- Filter funnel ---")
    log_lines.append(f"{STUDY_A_NAME} significant genes: {len(summary_a)} (up={a_up}, down={a_down})")
    log_lines.append(f"{STUDY_B_NAME} significant genes: {len(summary_b)} (up={b_up}, down={b_down})")

    # All locus tags in each study (for classifying one-study-only genes)
    study_a_all_tags = set(de_a["locus_tag"].unique())
    study_b_all_tags = set(de_b["locus_tag"].unique())
    log_lines.append(f"{STUDY_A_NAME} total locus tags in dataset: {len(study_a_all_tags)}")
    log_lines.append(f"{STUDY_B_NAME} total locus tags in dataset: {len(study_b_all_tags)}")

    # Intersect
    core, extended, discordant = intersect_references(
        summary_a, summary_b,
        study_a_name=STUDY_A_NAME,
        study_b_name=STUDY_B_NAME,
        study_a_all_locus_tags=study_a_all_tags,
        study_b_all_locus_tags=study_b_all_tags,
    )

    core_up = (core["direction"] == "up").sum()
    core_down = (core["direction"] == "down").sum()
    print(f"\nCore signature: {len(core)} genes (up={core_up}, down={core_down})")
    print(f"Extended: {len(extended)} genes")
    print(f"Discordant: {len(discordant)} genes")
    log_lines.append(f"\nCore: {len(core)} genes (up={core_up}, down={core_down})")
    log_lines.append(f"Extended: {len(extended)} genes")
    log_lines.append(f"Discordant: {len(discordant)} genes")

    if "signature_type" in extended.columns:
        ext_types = extended["signature_type"].value_counts().to_dict()
        print(f"  Extended by type: {ext_types}")
        log_lines.append(f"  {ext_types}")

    # Save
    save_csv(core, "core_signature.csv")
    save_csv(extended, "extended_signature.csv")
    save_csv(discordant, "discordant_genes.csv")
    print(f"\nSaved: data/core_signature.csv, data/extended_signature.csv, data/discordant_genes.csv")

    # Top 10 core genes
    print("\n=== Top 10 core genes by cross-study rank ===")
    top_cols = ["locus_tag", "gene_name", "direction", "cross_study_best_dir_rank"]
    if "product" in core.columns:
        top_cols.append("product")
    print(core[top_cols].head(10).to_string(index=False))
    log_lines.append("\nTop 10 core genes:")
    log_lines.append(core[top_cols].head(10).to_string(index=False))

    # Marker gene traces
    if args.explore:
        print("\n=== Marker gene traces ===")
        log_lines.append("\n=== Marker gene traces ===")
        for tag in MARKERS:
            in_core = core[core["locus_tag"] == tag]
            in_ext = extended[extended["locus_tag"] == tag]
            if len(in_core) > 0:
                row = in_core.iloc[0]
                line = (f"  {tag} ({row.get('gene_name', '?')}): CORE, direction={row['direction']}, "
                        f"cross_rank={row['cross_study_best_dir_rank']}, "
                        f"{STUDY_A_NAME}_rank={row.get(f'{STUDY_A_NAME}_best_dir_rank', '?')}, "
                        f"{STUDY_B_NAME}_rank={row.get(f'{STUDY_B_NAME}_best_dir_rank', '?')}")
            elif len(in_ext) > 0:
                row = in_ext.iloc[0]
                line = (f"  {tag} ({row.get('gene_name', '?')}): EXTENDED ({row['signature_type']}), "
                        f"direction={row['direction']}")
            else:
                in_a = summary_a[summary_a["locus_tag"] == tag]
                in_b = summary_b[summary_b["locus_tag"] == tag]
                if len(in_a) > 0 and len(in_b) > 0:
                    line = (f"  {tag}: DISCORDANT — {STUDY_A_NAME}={in_a.iloc[0]['direction']}, "
                            f"{STUDY_B_NAME}={in_b.iloc[0]['direction']}")
                elif len(in_a) > 0:
                    line = f"  {tag}: in {STUDY_A_NAME} only (not in {STUDY_B_NAME} summary)"
                elif len(in_b) > 0:
                    line = f"  {tag}: in {STUDY_B_NAME} only"
                else:
                    line = f"  {tag}: NOT IN EITHER STUDY (not significant anywhere)"
            print(line)
            log_lines.append(line)

    write_log("\n".join(log_lines), "03_build_signature.log")
    print(f"\nDone. See logs/03_build_signature.log")


if __name__ == "__main__":
    main()
