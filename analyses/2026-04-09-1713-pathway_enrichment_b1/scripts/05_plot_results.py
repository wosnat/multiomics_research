"""Step 5: Plot enrichment heatmaps.

Loads enrichment_all.csv and generates heatmaps of -log10(padj) for
up- and down-regulated pathways across all experiments and timepoints.

Only pathways significant (padj < 0.05) in at least one condition are shown.

Usage:
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/05_plot_results.py
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/05_plot_results.py --explore
"""

import argparse
import datetime
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = ANALYSIS_DIR / "results"
LOGS_DIR = ANALYSIS_DIR / "logs"

PADJ_THRESHOLD = 0.05
NEG_LOG10_THRESHOLD = -np.log10(PADJ_THRESHOLD)  # ~1.30

# Colormap for -log10(padj): white (not significant) → deep color (highly significant)
CMAP_UP = "Reds"
CMAP_DOWN = "Blues"


def make_heatmap(
    matrix: pd.DataFrame,
    title: str,
    out_path: Path,
    cmap: str,
    threshold: float,
    log_fn,
) -> None:
    """Plot a pathway × condition heatmap of -log10(padj).

    Parameters
    ----------
    matrix:
        DataFrame with pathways as rows and condition labels as columns.
        Values are -log10(padj); 0 = not significant.
    title:
        Plot title.
    out_path:
        Output PNG path.
    cmap:
        Matplotlib colormap name.
    threshold:
        -log10(padj) at significance threshold (used as vmin).
    log_fn:
        Logging function.
    """
    n_pathways, n_conditions = matrix.shape
    log_fn(f"  Plotting {n_pathways} pathways × {n_conditions} conditions → {out_path.name}")

    # Figure size: scale with dimensions
    fig_width = max(8, n_conditions * 0.6 + 3)
    fig_height = max(4, n_pathways * 0.3 + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    im = ax.imshow(
        matrix.values,
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=max(matrix.values.max(), threshold + 1),
        interpolation="nearest",
    )

    # Threshold line (significance boundary)
    plt.colorbar(im, ax=ax, label="-log10(padj)")

    ax.set_xticks(range(n_conditions))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_pathways))
    ax.set_yticklabels(matrix.index, fontsize=7)

    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xlabel("Experiment · Timepoint", fontsize=9)
    ax.set_ylabel("Pathway", fontsize=9)

    # Draw significance threshold gridline on colorbar (visual guide)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log_fn(f"  Saved {out_path}")


def build_matrix(
    enrich_df: pd.DataFrame,
    direction: str,
    sig_pathways: set,
) -> pd.DataFrame:
    """Build pathway × condition matrix of -log10(padj).

    Parameters
    ----------
    enrich_df:
        Full enrichment results DataFrame.
    direction:
        'up' or 'down'.
    sig_pathways:
        Set of pathway_id values that are significant in at least one condition.

    Returns
    -------
    DataFrame with pathway_name as index and condition labels as columns.
    """
    subset = enrich_df[
        (enrich_df["direction"] == direction)
        & (enrich_df["pathway_id"].isin(sig_pathways))
    ].copy()

    # Build condition label: experiment + timepoint
    if "timepoint" in subset.columns and "experiment" in subset.columns:
        subset["condition"] = subset["experiment"].str.replace("de_", "", regex=False) + "·" + subset["timepoint"].astype(str)
    elif "experiment" in subset.columns:
        subset["condition"] = subset["experiment"].str.replace("de_", "", regex=False)
    else:
        subset["condition"] = "unknown"

    # -log10(padj), clamp NaN to 0
    subset["neg_log10_padj"] = subset["padj"].apply(
        lambda p: -np.log10(p) if pd.notna(p) and p > 0 else 0.0
    )

    # Use pathway_name as row label (fall back to pathway_id if missing)
    if "pathway_name" in subset.columns:
        subset["label"] = subset["pathway_name"].fillna(subset["pathway_id"])
    else:
        subset["label"] = subset["pathway_id"]

    if len(subset) == 0:
        return pd.DataFrame()

    matrix = subset.pivot_table(
        index="label",
        columns="condition",
        values="neg_log10_padj",
        aggfunc="max",
        fill_value=0.0,
    )

    # Sort pathways by max significance across conditions (descending)
    matrix = matrix.loc[matrix.max(axis=1).sort_values(ascending=False).index]

    return matrix


def main(explore: bool = False) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    log_path = LOGS_DIR / "05_plot_results.log"
    log_lines = []

    def log(msg: str) -> None:
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        log_lines.append(line)

    log("=== Step 5: Plot Enrichment Heatmaps ===")

    # ------------------------------------------------------------------
    # 1. Load enrichment results
    # ------------------------------------------------------------------
    enrich_path = RESULTS_DIR / "enrichment_all.csv"
    if not enrich_path.exists():
        raise FileNotFoundError(
            f"Enrichment results not found: {enrich_path}\n"
            f"Run 04_run_enrichment.py first."
        )
    enrich_df = pd.read_csv(enrich_path)
    log(f"Loaded enrichment_all.csv: {len(enrich_df)} rows")

    # ------------------------------------------------------------------
    # 2. Find pathways significant in at least one condition
    # ------------------------------------------------------------------
    sig_mask = enrich_df["padj"] < PADJ_THRESHOLD
    sig_pathway_ids = set(enrich_df.loc[sig_mask, "pathway_id"].dropna().unique())
    log(f"Pathways significant in >= 1 condition: {len(sig_pathway_ids)}")

    if len(sig_pathway_ids) == 0:
        log("WARNING: No significant pathways found. Heatmaps will be empty.")

    if explore:
        log("--- Significant pathway IDs ---")
        for pid in sorted(sig_pathway_ids):
            n_sig = sig_mask[enrich_df["pathway_id"] == pid].sum()
            log(f"  {pid}: significant in {n_sig} conditions")

    # ------------------------------------------------------------------
    # 3. Plot up direction
    # ------------------------------------------------------------------
    up_matrix = build_matrix(enrich_df, direction="up", sig_pathways=sig_pathway_ids)
    if len(up_matrix) > 0:
        make_heatmap(
            matrix=up_matrix,
            title="Pathway Enrichment — Up-regulated",
            out_path=RESULTS_DIR / "heatmap_enrichment_up.png",
            cmap=CMAP_UP,
            threshold=NEG_LOG10_THRESHOLD,
            log_fn=log,
        )
    else:
        log("No significant up-regulated pathways to plot")

    # ------------------------------------------------------------------
    # 4. Plot down direction
    # ------------------------------------------------------------------
    down_matrix = build_matrix(enrich_df, direction="down", sig_pathways=sig_pathway_ids)
    if len(down_matrix) > 0:
        make_heatmap(
            matrix=down_matrix,
            title="Pathway Enrichment — Down-regulated",
            out_path=RESULTS_DIR / "heatmap_enrichment_down.png",
            cmap=CMAP_DOWN,
            threshold=NEG_LOG10_THRESHOLD,
            log_fn=log,
        )
    else:
        log("No significant down-regulated pathways to plot")

    log("Done.")

    # Write log file
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot pathway enrichment heatmaps"
    )
    parser.add_argument(
        "--explore",
        action="store_true",
        help="Print list of significant pathway IDs before plotting",
    )
    args = parser.parse_args()
    main(explore=args.explore)
