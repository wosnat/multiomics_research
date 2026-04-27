"""Trajectory figures for step 5.

Outputs:
  - figures/trajectories_positive_panel.png        5 axes × 2 omics, axenic vs coculture overlaid (signed-Z)
  - figures/trajectories_positive_panel_raw.png    same as above but axis_mean signed log2FC (raw)
  - figures/trajectories_cyanorak_panel.png        cyanorak (broader gene set) signed-Z, where applicable
  - figures/panel_comparison.png                   positive-panel vs cyanorak-panel, scatter per (axis, condition, omics, TP)

Run:
  uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/5_analyze/scripts/02_plot_trajectories.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

AXIS_ORDER = ["n_stress", "photo", "proteotoxic", "oxidative", "cell_death"]
OMICS_ORDER = ["RNASEQ", "PROTEOMICS"]
COND_COLOR = {"axenic": "#d04545", "coculture": "#5b9bd5"}


def trajectories_grid(
    scores: pd.DataFrame, value_col: str, title: str, out_path: Path,
    *,
    axes_to_plot: list[str] | None = None,
    z2_line: bool = True,
) -> None:
    axes_to_plot = axes_to_plot or AXIS_ORDER
    n_rows = len(axes_to_plot)
    n_cols = len(OMICS_ORDER)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 2.4 * n_rows + 0.6),
                              sharex="col", sharey="row")
    if n_rows == 1:
        axes = np.array([axes])

    for i, ax_name in enumerate(axes_to_plot):
        for j, omics in enumerate(OMICS_ORDER):
            ax = axes[i, j]
            sub = scores[(scores["axis"] == ax_name) & (scores["omics"] == omics)]
            for cond in ("axenic", "coculture"):
                cs = sub[sub["condition"] == cond].copy()
                if cs.empty:
                    continue
                cs = cs.sort_values("timepoint_hours")
                x = (cs["timepoint_hours"].astype(float) / 24.0).to_numpy()
                y = cs[value_col].to_numpy()
                # Single-point experiments → big marker, no line
                if len(cs) == 1 and pd.isna(cs["timepoint_hours"].iloc[0]):
                    # No timepoint hours — place at x=0 for visibility, draw differently
                    ax.axhline(y[0], color=COND_COLOR[cond], linestyle=":",
                               linewidth=1.5, label=f"{cond} (single point)",
                               alpha=0.6)
                else:
                    ax.plot(x, y, marker="o", linewidth=1.5,
                            color=COND_COLOR[cond], label=cond,
                            markeredgecolor="black", markersize=6)
            ax.axhline(0, color="black", lw=0.6)
            if z2_line and value_col == "axis_score":
                ax.axhline(2, color="gray", lw=0.5, linestyle="--")
                ax.axhline(-2, color="gray", lw=0.5, linestyle="--")
            if i == 0:
                ax.set_title(omics, fontsize=10)
            if j == 0:
                ax.set_ylabel(ax_name, fontsize=10, fontweight="bold")
            if i == n_rows - 1:
                ax.set_xlabel("timepoint (days)", fontsize=9)
            ax.tick_params(axis="both", labelsize=8)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0), fontsize=10)

    fig.suptitle(title, fontsize=11, y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Wrote {out_path}")


def panel_comparison(scores: pd.DataFrame, out_path: Path) -> None:
    """Scatter of positive-panel score vs cyanorak-panel score, per (axis, cell, TP).

    Where they agree (close to y=x), the result is robust to gene-set choice.
    Where they diverge, the wider cyanorak set is washing out (or surfacing) signal
    that the focused positive panel doesn't see.
    """
    pivot = scores.pivot_table(
        index=["axis", "condition", "omics", "timepoint", "growth_phase"],
        columns="panel_kind",
        values="axis_score",
        aggfunc="first",
    ).reset_index()
    pivot = pivot.dropna(subset=["positive", "cyanorak"])  # only axes with both panels

    fig, ax = plt.subplots(figsize=(7, 6))
    cmap = {"n_stress": "#d04545", "photo": "#1f7a3e", "proteotoxic": "#d27000",
            "oxidative": "#7b3fa0", "cell_death": "#5b9bd5"}
    for axis_name, group in pivot.groupby("axis"):
        ax.scatter(group["cyanorak"], group["positive"],
                   s=70, color=cmap.get(axis_name, "gray"),
                   edgecolor="black", linewidth=0.5, label=axis_name, alpha=0.85)

    # y = x reference
    lo = min(pivot["positive"].min(), pivot["cyanorak"].min()) - 0.5
    hi = max(pivot["positive"].max(), pivot["cyanorak"].max()) + 0.5
    ax.plot([lo, hi], [lo, hi], "--", color="gray", lw=0.8, label="y = x")
    ax.axhline(0, color="black", lw=0.4)
    ax.axvline(0, color="black", lw=0.4)

    ax.set_xlabel("axis_score (cyanorak panel — broader, ~25-30 genes)")
    ax.set_ylabel("axis_score (positive panel — 3-5 validated genes)")
    ax.set_title("Sensitivity check: axis_score by panel kind\n"
                 "(each point = one (axis, condition, omics, TP) cell)")
    ax.legend(fontsize=9, loc="lower right")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Wrote {out_path}")


def main() -> None:
    scores = pd.read_csv(DATA_DIR / "all_axes_scores.csv")
    print(f"Loaded {len(scores)} score rows")

    pos = scores[scores["panel_kind"] == "positive"]
    cy = scores[scores["panel_kind"] == "cyanorak"]

    trajectories_grid(
        pos, value_col="axis_score",
        title="Stress-axis trajectories (positive panel, signed-Z vs genome background)\n"
              "Pro MED4, Weissberg 2025 — within-condition trajectories per omics",
        out_path=FIG_DIR / "trajectories_positive_panel.png",
    )
    trajectories_grid(
        pos, value_col="axis_mean_signed_lfc",
        title="Stress-axis raw response (positive panel, mean signed log2FC over axis genes)\n"
              "Pro MED4, Weissberg 2025 — within-condition trajectories per omics",
        out_path=FIG_DIR / "trajectories_positive_panel_raw.png",
        z2_line=False,
    )
    trajectories_grid(
        cy, value_col="axis_score",
        title="Stress-axis trajectories (cyanorak panel, signed-Z vs genome background)\n"
              "Pro MED4, Weissberg 2025 — broader gene sets per axis (no cell_death panel)",
        out_path=FIG_DIR / "trajectories_cyanorak_panel.png",
        axes_to_plot=[a for a in AXIS_ORDER if a != "cell_death"],
    )
    panel_comparison(scores, FIG_DIR / "panel_comparison.png")


if __name__ == "__main__":
    main()
