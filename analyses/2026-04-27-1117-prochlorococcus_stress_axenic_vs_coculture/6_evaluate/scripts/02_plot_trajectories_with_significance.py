"""Step 6 figure: trajectory plots with permutation-significance annotations.

Re-renders the step-5 headline trajectory layout, overlaying:
  - filled markers for cells significant after BH-FDR correction (p_bh < 0.05)
  - open markers for non-significant cells
  - * = p_bh < 0.05; ** = p_bh < 0.01; *** = p_bh < 0.001

Inputs:
  - 5_analyze/data/all_axes_scores.csv
  - data/permutation_pvalues.csv

Outputs:
  - figures/trajectories_positive_panel_sig.png + .pdf
  - figures/trajectories_cyanorak_panel_sig.png + .pdf

Run:
  uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/6_evaluate/scripts/02_plot_trajectories_with_significance.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parents[2]
DATA_5 = ANALYSIS_DIR / "5_analyze" / "data"
DATA_6 = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

AXIS_ORDER = ["n_stress", "photo", "proteotoxic", "oxidative", "cell_death"]
OMICS_ORDER = ["RNASEQ", "PROTEOMICS"]
COND_COLOR = {"axenic": "#d04545", "coculture": "#5b9bd5"}


def sig_stars(p_bh: float) -> str:
    if pd.isna(p_bh):
        return ""
    if p_bh < 0.001:
        return "***"
    if p_bh < 0.01:
        return "**"
    if p_bh < 0.05:
        return "*"
    return ""


def trajectories_grid_with_sig(
    scores: pd.DataFrame, pvals: pd.DataFrame, panel_kind: str,
    out_path: Path, title: str,
    *, axes_to_plot: list[str] | None = None,
) -> None:
    axes_to_plot = axes_to_plot or AXIS_ORDER
    n_rows = len(axes_to_plot)
    n_cols = len(OMICS_ORDER)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11, 2.5 * n_rows + 0.6),
                              sharex="col", sharey="row")
    if n_rows == 1:
        axes = np.array([axes])

    sub_pvals = pvals[pvals["panel_kind"] == panel_kind]
    sub_scores = scores[scores["panel_kind"] == panel_kind]

    for i, ax_name in enumerate(axes_to_plot):
        for j, omics in enumerate(OMICS_ORDER):
            ax = axes[i, j]
            cell = sub_scores[(sub_scores["axis"] == ax_name) & (sub_scores["omics"] == omics)]
            cell_p = sub_pvals[(sub_pvals["axis"] == ax_name) & (sub_pvals["omics"] == omics)]

            for cond in ("axenic", "coculture"):
                cs = cell[cell["condition"] == cond].copy()
                cps = cell_p[cell_p["condition"] == cond].copy()
                if cs.empty:
                    continue
                cs = cs.sort_values("timepoint_hours")
                cps = cps.sort_values("timepoint_hours")
                # Merge p-values onto scores on (timepoint)
                merged = cs.merge(cps[["timepoint", "permutation_p_bh", "sig_05"]], on="timepoint", how="left")

                x_days = (merged["timepoint_hours"].astype(float) / 24.0).to_numpy()
                y = merged["axis_score"].to_numpy()
                sig = merged["sig_05"].fillna(False).to_numpy().astype(bool)
                pbh = merged["permutation_p_bh"].to_numpy()

                # Single-point experiment (axenic-RNA)
                if len(merged) == 1 and pd.isna(merged["timepoint_hours"].iloc[0]):
                    ax.axhline(y[0], color=COND_COLOR[cond], linestyle=":",
                               linewidth=1.5, label=f"{cond} (single point)", alpha=0.6)
                    s = sig_stars(pbh[0])
                    if s:
                        ax.text(40, y[0], s, ha="center", va="bottom",
                                color="black", fontsize=10, fontweight="bold")
                else:
                    # Connect with line
                    ax.plot(x_days, y, "-", color=COND_COLOR[cond], linewidth=1.2, alpha=0.6)
                    # Significant cells: filled marker; non-sig: open marker
                    for xi, yi, si, pi in zip(x_days, y, sig, pbh):
                        if si:
                            ax.plot(xi, yi, "o", color=COND_COLOR[cond],
                                    markeredgecolor="black", markersize=8, markeredgewidth=0.8)
                        else:
                            ax.plot(xi, yi, "o", markerfacecolor="white",
                                    markeredgecolor=COND_COLOR[cond], markersize=8, markeredgewidth=1.2)
                        s = sig_stars(pi)
                        if s:
                            ax.text(xi, yi + 0.18, s, ha="center", va="bottom",
                                    color="black", fontsize=9, fontweight="bold")
                    # Add a single legend entry per condition
                    ax.plot([], [], "o-", color=COND_COLOR[cond],
                            markeredgecolor="black", markersize=7, label=cond)

            ax.axhline(0, color="black", lw=0.6)
            ax.axhline(2, color="gray", lw=0.5, linestyle="--", alpha=0.5)
            ax.axhline(-2, color="gray", lw=0.5, linestyle="--", alpha=0.5)
            if i == 0:
                ax.set_title(omics, fontsize=10)
            if j == 0:
                ax.set_ylabel(ax_name, fontsize=10, fontweight="bold")
            if i == n_rows - 1:
                ax.set_xlabel("timepoint (days)", fontsize=9)
            ax.tick_params(axis="both", labelsize=8)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3,
                   bbox_to_anchor=(0.5, 1.0), fontsize=10)

    fig.suptitle(
        f"{title}\n"
        "filled marker = sig at BH-FDR p<0.05; open = not sig; "
        "* p_bh<0.05; ** p_bh<0.01; *** p_bh<0.001",
        fontsize=10, y=1.02,
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Wrote {out_path}")


def main() -> None:
    scores = pd.read_csv(DATA_5 / "all_axes_scores.csv")
    pvals = pd.read_csv(DATA_6 / "permutation_pvalues.csv")
    print(f"Loaded {len(scores)} scores and {len(pvals)} p-values")

    trajectories_grid_with_sig(
        scores, pvals, panel_kind="positive",
        out_path=FIG_DIR / "trajectories_positive_panel_sig.png",
        title="Stress-axis trajectories with permutation significance (positive panel)\n"
              "Pro MED4, Weissberg 2025",
    )
    trajectories_grid_with_sig(
        scores, pvals, panel_kind="cyanorak",
        out_path=FIG_DIR / "trajectories_cyanorak_panel_sig.png",
        title="Stress-axis trajectories with permutation significance (cyanorak panel)\n"
              "Pro MED4, Weissberg 2025 — broader gene sets per axis",
        axes_to_plot=[a for a in AXIS_ORDER if a != "cell_death"],
    )


if __name__ == "__main__":
    main()
