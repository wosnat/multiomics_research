"""Step 6: Visualize scoring results.

Produces:
    results/trajectory_rnaseq.png — rank score over time, RNA-seq
    results/trajectory_proteomics.png — rank score over time, proteomics
    results/control_separation.png — all experiments by role
    results/tier_comparison.png — top vs core vs extended for targets

Run from multiomics_research root:
    uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/06_plot_results.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sig_utils.io import RESULTS_DIR

RESULTS_DIR.mkdir(exist_ok=True)


def load_scores():
    return pd.read_csv(RESULTS_DIR / "scores_all.csv")


def plot_trajectories(scores_df: pd.DataFrame):
    """Plot rank score trajectories for RNA-seq and proteomics (core tier)."""
    core = scores_df[scores_df["tier"] == "core"].copy()

    # Reference range for shading
    ref = core[core["role"] == "reference"]
    ref_min = ref["rank_score"].min() if len(ref) > 0 else 0
    ref_max = ref["rank_score"].max() if len(ref) > 0 else 0

    for platform, experiments in [
        ("RNA-seq", ["Weissberg RNA-seq axenic", "Weissberg RNA-seq coculture"]),
        ("Proteomics", ["Weissberg proteomics axenic", "Weissberg proteomics coculture"]),
    ]:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

        for ax_idx, (ax, metric, ylabel) in enumerate(zip(
            axes, ["rank_score", "hit_rate"],
            ["Rank Score", "Hit Rate (concordant fraction)"]
        )):
            # Reference band
            if metric == "rank_score":
                ax.axhspan(ref_min, ref_max, alpha=0.08, color="green",
                           label="Reference range")

            ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)

            colors = {"axenic": "#d62728", "coculture": "#1f77b4"}
            for exp_label in experiments:
                subset = core[core["label"] == exp_label].copy()
                if subset.empty:
                    continue

                condition = "axenic" if "axenic" in exp_label.lower() else "coculture"
                color = colors[condition]

                # Time-course points
                tc = subset[subset["timepoint_hours"].notna()].sort_values("timepoint_hours")
                # Single timepoint
                single = subset[subset["timepoint_hours"].isna()]

                if len(tc) > 0:
                    # Skip "days 60+89" combined timepoint
                    tc = tc[~tc["timepoint"].str.contains(r"\+", na=False)]
                    ax.plot(tc["timepoint_hours"], tc[metric],
                            marker="o", color=color, linewidth=2,
                            label=f"{condition}")

                if len(single) > 0:
                    # Plot single timepoint as a star at a reasonable x position
                    # Use first coculture timepoint x for visual alignment
                    x_pos = 14 * 24  # day 14 in hours
                    for _, row in single.iterrows():
                        ax.scatter([x_pos], [row[metric]],
                                   marker="*", s=200, color=color, zorder=5,
                                   label=f"{condition} (single tp)")

            ax.set_xlabel("Time (hours)")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{platform} — {ylabel}")
            ax.legend(loc="best", fontsize=9)

        fig.suptitle(f"N-Limitation Signature Score — {platform}", fontsize=14, y=1.02)
        fig.tight_layout()

        outpath = RESULTS_DIR / f"trajectory_{platform.lower().replace('-', '')}.png"
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {outpath.name}")


def plot_control_separation(scores_df: pd.DataFrame):
    """Plot all experiments by role (core tier), with breakdown annotation."""
    core = scores_df[scores_df["tier"] == "core"].copy()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: rank score
    # Right: hit rate
    for ax, metric, title in zip(
        axes,
        ["rank_score", "hit_rate"],
        ["Rank Score", "Hit Rate"],
    ):
        role_colors = {
            "reference": "#2ca02c",
            "negative_control": "#d62728",
            "target": "#1f77b4",
        }

        for role in ["reference", "negative_control", "target"]:
            subset = core[core["role"] == role]
            if subset.empty:
                continue

            # Jitter x position
            x_base = list(role_colors.keys()).index(role)
            x_jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(subset))

            ax.scatter(
                x_base + x_jitter, subset[metric],
                c=role_colors[role], s=60, alpha=0.7,
                label=role.replace("_", " "),
            )

            # Annotate outliers
            for _, row in subset.iterrows():
                if (role == "negative_control" and row[metric] > 0.3) or \
                   (role == "target" and row[metric] > 0.5):
                    ax.annotate(
                        f"{row['label']}\n{row['timepoint']}",
                        (x_base + x_jitter[subset.index.get_loc(row.name)], row[metric]),
                        fontsize=7, ha="center", va="bottom",
                    )

        ax.set_xticks(range(3))
        ax.set_xticklabels(["Reference", "Negative\ncontrol", "Target"])
        ax.set_ylabel(metric.replace("_", " ").title())
        ax.set_title(title)
        ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)

    fig.suptitle("Control Separation — Core Signature (189 genes)", fontsize=13)
    fig.tight_layout()

    outpath = RESULTS_DIR / "control_separation.png"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath.name}")


def plot_tier_comparison(scores_df: pd.DataFrame):
    """Plot top vs core vs extended for target experiments."""
    targets = scores_df[scores_df["role"] == "target"].copy()

    fig, ax = plt.subplots(figsize=(12, 6))

    tier_colors = {"top": "#ff7f0e", "core": "#1f77b4", "extended": "#2ca02c"}
    tier_offsets = {"top": -0.2, "core": 0.0, "extended": 0.2}

    # Create x positions from unique label+timepoint combos
    targets["x_label"] = targets["label"].str.replace("Weissberg ", "") + "\n" + targets["timepoint"].astype(str)
    x_labels = targets.drop_duplicates("x_label")["x_label"].tolist()
    x_map = {label: i for i, label in enumerate(x_labels)}
    targets["x_pos"] = targets["x_label"].map(x_map)

    for tier in ["top", "core", "extended"]:
        tier_data = targets[targets["tier"] == tier]
        ax.scatter(
            tier_data["x_pos"] + tier_offsets[tier],
            tier_data["rank_score"],
            c=tier_colors[tier], s=80, alpha=0.8,
            label=f"{tier} ({len(tier_data.iloc[0:1])})" if len(tier_data) > 0 else tier,
            zorder=3,
        )

    # Add tier sizes to legend
    ax.legend(
        [f"top (38 genes)", f"core (189 genes)", f"extended (531 genes)"],
        loc="upper right", fontsize=9,
    )

    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Rank Score")
    ax.set_title("Tier Comparison — Target Experiments")
    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)

    fig.tight_layout()
    outpath = RESULTS_DIR / "tier_comparison.png"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath.name}")


def main():
    print("Loading scores...")
    scores = load_scores()
    print(f"  {len(scores)} score rows loaded")

    print("\nPlotting trajectories...")
    plot_trajectories(scores)

    print("\nPlotting control separation...")
    plot_control_separation(scores)

    print("\nPlotting tier comparison...")
    plot_tier_comparison(scores)

    print("\nDone.")


if __name__ == "__main__":
    main()
