"""Plot activation score trajectories for axenic vs coculture.

Reads: results/signature_scores_core.csv, results/signature_scores_extended.csv
Outputs: results/trajectory_rnaseq.png, results/trajectory_proteomics.png

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/plot_trajectories.py
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


def load_weissberg_scores(filename: str) -> pd.DataFrame:
    """Load scores and filter to Weissberg experiments with numeric timepoints."""
    df = pd.read_csv(RESULTS_DIR / filename)
    df = df[df["study"] == "weissberg_2025"].copy()
    # Exclude combined timepoints (d60+89)
    df = df[~df["timepoint"].str.contains(r"\+", na=False)]
    # Parse day number from timepoint strings like "day 18", "day 31"
    df["day"] = df["timepoint"].str.extract(r"(\d+)").astype(float)
    return df


def load_reference_scores(filename: str) -> pd.DataFrame:
    """Load reference baseline scores."""
    df = pd.read_csv(RESULTS_DIR / filename)
    df = df[df["study"] != "weissberg_2025"].copy()
    return df


METRICS = [
    ("hit_rate_concordant", "Hit rate (concordant)", 0, 1),
    ("mean_signed_log2fc", "Mean signed log2FC", None, None),
    ("rank_score", "Rank score", None, None),
]


def plot_trajectories_for_platform(
    platform: str,
    core_df: pd.DataFrame,
    extended_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    output_path: Path,
):
    """Plot 3-panel figure (one per metric) for a platform."""
    core_plat = core_df[core_df["platform"] == platform]
    ext_plat = extended_df[extended_df["platform"] == platform]
    ref_plat = ref_df[ref_df["study"].str.contains("tolonen" if platform == "microarray" else "")]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"N-Limitation Signature Activation — {platform.upper()}", fontsize=14)

    for ax, (metric, label, ymin, ymax) in zip(axes, METRICS):
        # Weissberg core
        for condition, style, color in [("axenic", "--", "red"), ("coculture", "-", "blue")]:
            subset = core_plat[core_plat["condition"] == condition].sort_values("day")
            if not subset.empty:
                ax.plot(subset["day"], subset[metric], style, color=color,
                        marker="o", label=f"{condition} (core)", linewidth=2)

        # Weissberg extended (lighter)
        for condition, style, color in [("axenic", "--", "salmon"), ("coculture", "-", "lightblue")]:
            subset = ext_plat[ext_plat["condition"] == condition].sort_values("day")
            if not subset.empty:
                ax.plot(subset["day"], subset[metric], style, color=color,
                        marker="s", alpha=0.6, label=f"{condition} (extended)", linewidth=1)

        # Reference baselines as horizontal bands
        for study_label, ref_color in [("tolonen_ndep", "gray"), ("read_ndep", "darkgray")]:
            ref_study = ref_plat[ref_plat["study"] == study_label]
            if not ref_study.empty:
                ref_min = ref_study[metric].min()
                ref_max = ref_study[metric].max()
                ax.axhspan(ref_min, ref_max, alpha=0.15, color=ref_color,
                          label=f"{study_label} range")

        ax.set_xlabel("Day")
        ax.set_ylabel(label)
        if ymin is not None:
            ax.set_ylim(ymin, ymax)
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def main():
    core_w = load_weissberg_scores("signature_scores_core.csv")
    ext_w = load_weissberg_scores("signature_scores_extended.csv")
    core_all = pd.read_csv(RESULTS_DIR / "signature_scores_core.csv")
    ref = core_all[core_all["study"] != "weissberg_2025"]

    for platform in ["rnaseq", "proteomics"]:
        output = RESULTS_DIR / f"trajectory_{platform}.png"
        plot_trajectories_for_platform(platform, core_w, ext_w, ref, output)

    print("Done.")


if __name__ == "__main__":
    main()
