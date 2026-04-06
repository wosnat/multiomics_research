"""Score Weissberg and reference experiments against N-limitation signatures.

Reads: data/core_signature_genes.csv, data/extended_signature_genes.csv,
       data/de_weissberg_med4_*.csv, data/de_reference_*.csv
Outputs: results/signature_scores_core.csv, results/signature_scores_extended.csv,
         results/permutation_tests.csv, results/core_vs_extended_comparison.csv

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/score_signature.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
import numpy as np
from sig_utils.io import load_de_csv, load_signature_csv, save_scores_csv
from sig_utils.metrics import (
    hit_rate,
    mean_signed_log2fc,
    mean_signed_normalized_rank,
    permutation_test_mean_signed_log2fc,
)

# Experiment files and their metadata
WEISSBERG_EXPERIMENTS = [
    ("de_weissberg_med4_rnaseq_axenic.csv", "rnaseq", "axenic"),
    ("de_weissberg_med4_rnaseq_coculture.csv", "rnaseq", "coculture"),
    ("de_weissberg_med4_proteomics_axenic.csv", "proteomics", "axenic"),
    ("de_weissberg_med4_proteomics_coculture.csv", "proteomics", "coculture"),
]

REFERENCE_EXPERIMENTS = [
    ("de_reference_tolonen_ndep.csv", "microarray", "tolonen_ndep"),
    ("de_reference_read_ndep.csv", "rnaseq", "read_ndep"),
]

# Timepoints to skip
SKIP_TIMEPOINTS_TOLONEN = {"0h", "3h"}


def score_one_timepoint(
    signature_df: pd.DataFrame,
    de_timepoint_df: pd.DataFrame,
    total_genes: int,
) -> dict:
    """Compute all metrics for one signature x one timepoint."""
    hr = hit_rate(signature_df, de_timepoint_df)
    ms = mean_signed_log2fc(signature_df, de_timepoint_df)
    mr = mean_signed_normalized_rank(signature_df, de_timepoint_df, total_genes)

    return {
        "hit_rate_concordant": hr["hit_rate_concordant"],
        "hit_rate_reversed": hr["hit_rate_reversed"],
        "concordant_hits": hr["concordant_hits"],
        "reversed_hits": hr["reversed_hits"],
        "not_significant": hr["not_significant"],
        "mean_signed_log2fc": ms["score"],
        "n_detected": ms["n_detected"],
        "rank_score": mr["score"],
        "n_concordant_rank": mr["n_concordant"],
        "n_reversed_rank": mr["n_reversed"],
    }


def score_experiment(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
    platform: str,
    condition: str,
    study: str,
    skip_timepoints: set | None = None,
) -> list[dict]:
    """Score all timepoints in an experiment."""
    rows = []

    if "timepoint" in de_df.columns and de_df["timepoint"].nunique() > 1:
        for tp, tp_df in de_df.groupby("timepoint"):
            if skip_timepoints and tp in skip_timepoints:
                continue
            total_genes = tp_df["locus_tag"].nunique()
            scores = score_one_timepoint(signature_df, tp_df, total_genes)
            scores.update({
                "study": study,
                "platform": platform,
                "condition": condition,
                "timepoint": tp,
                "timepoint_hours": tp_df["timepoint_hours"].iloc[0] if "timepoint_hours" in tp_df.columns else np.nan,
                "total_genes_in_experiment": total_genes,
                "signature_genes_total": len(signature_df),
            })
            rows.append(scores)
    else:
        total_genes = de_df["locus_tag"].nunique()
        scores = score_one_timepoint(signature_df, de_df, total_genes)
        scores.update({
            "study": study,
            "platform": platform,
            "condition": condition,
            "timepoint": "single",
            "timepoint_hours": np.nan,
            "total_genes_in_experiment": total_genes,
            "signature_genes_total": len(signature_df),
        })
        rows.append(scores)

    return rows


def main():
    # Load signatures
    core = load_signature_csv("core_signature_genes.csv")
    extended = load_signature_csv("extended_signature_genes.csv")
    print(f"Core signature: {len(core)} genes")
    print(f"Extended signature: {len(extended)} genes")

    for sig_name, sig_df, output_file in [
        ("core", core, "signature_scores_core.csv"),
        ("extended", extended, "signature_scores_extended.csv"),
    ]:
        print(f"\n=== Scoring {sig_name} signature ===")
        all_rows = []

        # Score Weissberg experiments
        for filename, platform, condition in WEISSBERG_EXPERIMENTS:
            print(f"  {filename}...")
            de = load_de_csv(filename)
            rows = score_experiment(sig_df, de, platform, condition, "weissberg_2025")
            all_rows.extend(rows)

        # Score reference experiments (baselines)
        for filename, platform, study in REFERENCE_EXPERIMENTS:
            print(f"  {filename} (baseline)...")
            de = load_de_csv(filename)
            skip = SKIP_TIMEPOINTS_TOLONEN if "tolonen" in filename else None
            rows = score_experiment(sig_df, de, platform, "reference", study, skip)
            all_rows.extend(rows)

        scores_df = pd.DataFrame(all_rows)
        save_scores_csv(scores_df, output_file)
        print(f"  Saved {len(scores_df)} rows to results/{output_file}")

    # Permutation tests for core signature (Weissberg only)
    print("\n=== Permutation tests (core, mean_signed_log2fc) ===")
    perm_rows = []
    for filename, platform, condition in WEISSBERG_EXPERIMENTS:
        de = load_de_csv(filename)

        if "timepoint" in de.columns and de["timepoint"].nunique() > 1:
            for tp, tp_df in de.groupby("timepoint"):
                tp_fc = tp_df.set_index("locus_tag")["log2fc"]
                result = permutation_test_mean_signed_log2fc(core, tp_df, tp_fc)
                result.update({"platform": platform, "condition": condition, "timepoint": tp})
                perm_rows.append(result)
        else:
            all_fc = de.set_index("locus_tag")["log2fc"]
            result = permutation_test_mean_signed_log2fc(core, de, all_fc)
            result.update({"platform": platform, "condition": condition, "timepoint": "single"})
            perm_rows.append(result)

    perm_df = pd.DataFrame(perm_rows)
    save_scores_csv(perm_df, "permutation_tests.csv")
    print(f"  Saved {len(perm_df)} permutation tests")

    # Core vs extended comparison
    print("\n=== Core vs extended comparison ===")
    core_scores = pd.read_csv(ANALYSIS_DIR / "results" / "signature_scores_core.csv")
    ext_scores = pd.read_csv(ANALYSIS_DIR / "results" / "signature_scores_extended.csv")

    core_scores = core_scores.rename(columns={
        "hit_rate_concordant": "core_hit_rate",
        "mean_signed_log2fc": "core_mean_signed_log2fc",
        "rank_score": "core_rank_score",
    })
    ext_scores = ext_scores.rename(columns={
        "hit_rate_concordant": "extended_hit_rate",
        "mean_signed_log2fc": "extended_mean_signed_log2fc",
        "rank_score": "extended_rank_score",
    })

    merge_keys = ["study", "platform", "condition", "timepoint"]
    comparison = core_scores[merge_keys + ["core_hit_rate", "core_mean_signed_log2fc", "core_rank_score"]].merge(
        ext_scores[merge_keys + ["extended_hit_rate", "extended_mean_signed_log2fc", "extended_rank_score"]],
        on=merge_keys,
    )
    save_scores_csv(comparison, "core_vs_extended_comparison.csv")
    print(f"  Saved {len(comparison)} comparison rows")

    print("\nDone.")


if __name__ == "__main__":
    main()
