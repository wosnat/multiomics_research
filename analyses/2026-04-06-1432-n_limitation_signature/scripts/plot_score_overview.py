"""Overview strip plot of signature scores across all conditions.

Shows all timepoints as dots, grouped by condition, colored by platform,
split into core vs extended signature panels.

Reads: results/signature_scores_core.csv, results/signature_scores_extended.csv
Outputs: results/score_overview.png

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/plot_score_overview.py
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

METRICS = [
    ("hit_rate_concordant", "Hit rate (concordant)"),
    ("mean_signed_log2fc", "Mean signed log2FC"),
    ("rank_score", "Rank score"),
]

PLATFORM_COLORS = {
    "microarray": "#888888",
    "rnaseq": "#2171b5",
    "proteomics": "#e6550d",
}

# Map study/condition to display groups
def assign_group(row):
    if row["study"] != "weissberg_2025":
        if "tolonen" in row["study"]:
            return "Tolonen ref"
        else:
            return "Read ref"
    return row["condition"].capitalize()

GROUP_ORDER = ["Tolonen ref", "Read ref", "Axenic", "Coculture"]


def load_and_prepare(filename: str) -> pd.DataFrame:
    df = pd.read_csv(RESULTS_DIR / filename)
    # Exclude combined timepoints
    df = df[~df["timepoint"].str.contains(r"\+", na=False)]
    df["group"] = df.apply(assign_group, axis=1)
    return df


def plot_overview(core_df, ext_df, output_path):
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    fig.suptitle("N-Limitation Signature Scores — All Conditions", fontsize=14, y=0.98)

    for row_idx, (metric, label) in enumerate(METRICS):
        for col_idx, (sig_name, sig_df) in enumerate([("Core (198 genes)", core_df),
                                                        ("Extended (367 genes)", ext_df)]):
            ax = axes[row_idx, col_idx]

            for gi, group in enumerate(GROUP_ORDER):
                group_data = sig_df[sig_df["group"] == group]
                if group_data.empty:
                    continue

                for platform, pdf in group_data.groupby("platform"):
                    color = PLATFORM_COLORS.get(platform, "black")
                    # Jitter x position
                    jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(pdf))
                    ax.scatter(
                        gi + jitter,
                        pdf[metric],
                        c=color, s=50, alpha=0.7, edgecolors="white", linewidth=0.5,
                        label=platform if row_idx == 0 and gi == 0 else None,
                        zorder=3,
                    )

                # Box whisker for groups with enough points
                vals = group_data[metric].dropna()
                if len(vals) >= 3:
                    bp = ax.boxplot(
                        [vals], positions=[gi], widths=0.4,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor="none", edgecolor="gray", linewidth=0.8),
                        whiskerprops=dict(color="gray", linewidth=0.8),
                        medianprops=dict(color="gray", linewidth=1),
                        capprops=dict(color="gray", linewidth=0.8),
                    )

            ax.set_xticks(range(len(GROUP_ORDER)))
            ax.set_xticklabels(GROUP_ORDER, rotation=30, ha="right", fontsize=9)
            ax.set_ylabel(label, fontsize=10)
            ax.axhline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.3)
            ax.grid(True, axis="y", alpha=0.2)

            if row_idx == 0:
                ax.set_title(sig_name, fontsize=11, fontweight="bold")

    # Legend for platforms
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=8, label=p)
        for p, c in PLATFORM_COLORS.items()
    ]
    fig.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(0.98, 0.96),
              fontsize=9, title="Platform")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def main():
    core = load_and_prepare("signature_scores_core.csv")
    ext = load_and_prepare("signature_scores_extended.csv")

    output = RESULTS_DIR / "score_overview.png"
    plot_overview(core, ext, output)
    print("Done.")


if __name__ == "__main__":
    main()
