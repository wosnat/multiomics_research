"""Step 5: Score all experiments against N-limitation signature (3 tiers).

Scores references (positive controls), negative controls, and targets
using top/core/extended signature tiers.

Outputs:
    data/applied_*.csv — signature applied to each experiment
    results/scores_all.csv — all scores with tier, role, breakdown
    logs/05_score_experiments.log

Run from multiomics_research root:
    uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/05_score_experiments.py [--explore]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
import numpy as np
from sig_utils.scoring import apply_signature, rank_score, score_with_significance
from sig_utils.extraction import check_marker_genes
from sig_utils.io import load_de_csv, load_signature_csv, save_csv, write_log

MARKERS = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]
TOP_RANK_THRESHOLD = 10
N_PERMS = 1000
SEED = 42

# All experiments to score: (de_filename, label, role, skip_timepoints)
ALL_EXPERIMENTS = [
    # References (positive controls — self-scoring, circularity flagged)
    ("de_ref_tolonen_ndep.csv", "Tolonen N-dep", "reference", {"0h", "3h"}),
    ("de_ref_read_ndep.csv", "Read N-dep", "reference", {"3h"}),

    # Negative controls
    ("de_ctrl_tolonen_cyanate.csv", "Tolonen cyanate", "negative_control", set()),
    ("de_ctrl_tolonen_urea.csv", "Tolonen urea", "negative_control", set()),
    ("de_ctrl_aharonovich_coculture.csv", "Aharonovich coculture", "negative_control", set()),
    ("de_ctrl_steglich_high_white_light.csv", "Steglich high white light", "negative_control", set()),

    # Targets
    ("de_weissberg_rnaseq_axenic.csv", "Weissberg RNA-seq axenic", "target", set()),
    ("de_weissberg_rnaseq_coculture.csv", "Weissberg RNA-seq coculture", "target", set()),
    ("de_weissberg_proteomics_axenic.csv", "Weissberg proteomics axenic", "target", set()),
    ("de_weissberg_proteomics_coculture.csv", "Weissberg proteomics coculture", "target", set()),
]


def build_tiers():
    """Load signature and build top/core/extended tiers."""
    core = load_signature_csv("core_signature.csv")
    extended = load_signature_csv("extended_signature.csv")

    top = core[core["cross_study_best_dir_rank"] <= TOP_RANK_THRESHOLD].copy()

    # Extended = core + extended combined
    full_extended = pd.concat([core, extended], ignore_index=True)

    return {
        "top": top[["locus_tag", "gene_name", "direction"]],
        "core": core[["locus_tag", "gene_name", "direction"]],
        "extended": full_extended[["locus_tag", "gene_name", "direction"]],
    }


def score_experiment_all_tiers(
    tiers: dict,
    de_df: pd.DataFrame,
    label: str,
    role: str,
    skip_timepoints: set,
) -> tuple[list[dict], dict[str, list[pd.DataFrame]]]:
    """Score all timepoints for one experiment across all tiers.

    Returns:
        (score_rows, applied_by_tier) — score rows for results table,
        and {tier_name: [applied_dfs]} for saving.
    """
    score_rows = []
    applied_by_tier = {t: [] for t in tiers}

    # Determine timepoints
    if "timepoint" in de_df.columns and de_df["timepoint"].nunique() > 1:
        timepoints = [
            (tp, tp_df)
            for tp, tp_df in de_df.groupby("timepoint")
            if tp not in skip_timepoints
        ]
    else:
        tp_label = de_df["timepoint"].iloc[0] if "timepoint" in de_df.columns else "single"
        if tp_label in skip_timepoints:
            return score_rows, applied_by_tier
        timepoints = [(tp_label, de_df)]

    for tp, tp_df in timepoints:
        tp_hours = (
            tp_df["timepoint_hours"].iloc[0]
            if "timepoint_hours" in tp_df.columns and pd.notna(tp_df["timepoint_hours"].iloc[0])
            else np.nan
        )
        total_genes = tp_df["locus_tag"].nunique()

        for tier_name, tier_sig in tiers.items():
            result = score_with_significance(
                tier_sig, tp_df, n_perms=N_PERMS, seed=SEED
            )
            result.update({
                "tier": tier_name,
                "label": label,
                "role": role,
                "timepoint": tp,
                "timepoint_hours": tp_hours,
                "total_genes_in_experiment": total_genes,
            })
            score_rows.append(result)

            applied = apply_signature(tier_sig, tp_df)
            applied_by_tier[tier_name].append(
                applied.assign(timepoint=tp, label=label)
            )

    return score_rows, applied_by_tier


def main():
    parser = argparse.ArgumentParser(description="Score all experiments")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces")
    args = parser.parse_args()

    log_lines = [f"Step 5: Score All Experiments — {datetime.now().isoformat()}", ""]

    tiers = build_tiers()
    for name, sig in tiers.items():
        n_up = (sig["direction"] == "up").sum()
        n_down = (sig["direction"] == "down").sum()
        log_lines.append(f"Tier '{name}': {len(sig)} genes ({n_up} up, {n_down} down)")
        print(f"Tier '{name}': {len(sig)} genes ({n_up} up, {n_down} down)")
    print()

    all_score_rows = []

    for de_filename, label, role, skip_tps in ALL_EXPERIMENTS:
        print(f"Scoring {label} ({role})...")
        de = load_de_csv(de_filename)

        score_rows, applied_by_tier = score_experiment_all_tiers(
            tiers, de, label, role, skip_tps
        )
        all_score_rows.extend(score_rows)

        # Save applied subsets (core tier only — most useful for inspection)
        core_applied = applied_by_tier.get("core", [])
        if core_applied:
            combined = pd.concat(core_applied, ignore_index=True)
            save_csv(combined, f"applied_{de_filename}")

        # Print core tier summary
        core_rows = [r for r in score_rows if r["tier"] == "core"]
        for row in core_rows:
            p_str = f"{row['empirical_p']:.3f}" if not np.isnan(row['empirical_p']) else "NaN"
            print(f"  {row['timepoint']}: score={row['rank_score']:.3f}, p={p_str}, "
                  f"conc={row['n_concordant']}, rev={row['n_reversed']}, "
                  f"ns={row['n_not_significant']}, absent={row['n_absent']}, "
                  f"hit={row['hit_rate']:.2f}")

        if args.explore and core_applied:
            combined = pd.concat(core_applied, ignore_index=True)
            marker_rows = combined[combined["locus_tag"].isin(MARKERS)]
            if len(marker_rows) > 0:
                cols = ["locus_tag", "direction", "timepoint", "log2fc",
                        "concordance", "dir_rank"]
                cols = [c for c in cols if c in marker_rows.columns]
                log_lines.append(f"\nMarker genes in {label}:")
                log_lines.append(marker_rows[cols].to_string(index=False))

        log_lines.append("")

    # Save all scores
    scores_df = pd.DataFrame(all_score_rows)
    save_csv(scores_df, "scores_all.csv", subdir="results")
    print(f"\nSaved: results/scores_all.csv ({len(scores_df)} rows)")

    # Control separation summary
    print("\n=== Control separation (core tier) ===")
    log_lines.append("\n=== Control separation (core tier) ===")
    core_scores = scores_df[scores_df["tier"] == "core"]
    for role_name in ["reference", "negative_control", "target"]:
        role_df = core_scores[core_scores["role"] == role_name]
        if role_df.empty:
            continue
        line = (f"  {role_name}: n={len(role_df)}, "
                f"mean_score={role_df['rank_score'].mean():.3f}, "
                f"range=[{role_df['rank_score'].min():.3f}, {role_df['rank_score'].max():.3f}], "
                f"mean_hit_rate={role_df['hit_rate'].mean():.2f}")
        print(line)
        log_lines.append(line)

    # Tier comparison
    print("\n=== Tier comparison (targets only) ===")
    log_lines.append("\n=== Tier comparison (targets only) ===")
    target_scores = scores_df[scores_df["role"] == "target"]
    for tier_name in ["top", "core", "extended"]:
        tier_df = target_scores[target_scores["tier"] == tier_name]
        if tier_df.empty:
            continue
        line = (f"  {tier_name}: mean_score={tier_df['rank_score'].mean():.3f}, "
                f"range=[{tier_df['rank_score'].min():.3f}, {tier_df['rank_score'].max():.3f}], "
                f"mean_hit_rate={tier_df['hit_rate'].mean():.2f}")
        print(line)
        log_lines.append(line)

    write_log("\n".join(log_lines), "05_score_experiments.log")
    print(f"\nDone. See logs/05_score_experiments.log")


if __name__ == "__main__":
    main()
