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

# Short labels for experiments (for x-axis)
EXPERIMENT_SHORT = {
    "de_ref_tolonen_ndep": "Tolonen N-dep",
    "de_ref_read_ndep": "Read N-dep",
    "de_ctrl_tolonen_cyanate": "Tolonen cyanate",
    "de_ctrl_tolonen_urea": "Tolonen urea",
    "de_ctrl_aharonovich_coculture": "Aharonovich cocult",
    "de_ctrl_steglich_high_white_light": "Steglich light",
    "de_weissberg_rnaseq_axenic": "Wb RNA ax",
    "de_weissberg_rnaseq_coculture": "Wb RNA cocult",
    "de_weissberg_proteomics_axenic": "Wb prot ax",
    "de_weissberg_proteomics_coculture": "Wb prot cocult",
}

# Order: references → negative controls → targets (RNA then prot, axenic then cocult)
EXPERIMENT_ORDER = [
    "de_ref_tolonen_ndep",
    "de_ref_read_ndep",
    "de_ctrl_tolonen_cyanate",
    "de_ctrl_tolonen_urea",
    "de_ctrl_aharonovich_coculture",
    "de_ctrl_steglich_high_white_light",
    "de_weissberg_rnaseq_axenic",
    "de_weissberg_rnaseq_coculture",
    "de_weissberg_proteomics_axenic",
    "de_weissberg_proteomics_coculture",
]

# Timepoint sort order
TIMEPOINT_ORDER = {
    "0h": 0, "3h": 1, "6h": 2, "12h": 3, "24h": 4, "48h": 5,
    "20h": 3,
    "day 14": 10, "day 18": 11, "day 31": 12, "day 60": 13, "day 89": 14,
    "days 60+89": 15,
}

# Short pathway names: strip the "Category > " prefix, keep the specific part
def _short_pathway_name(name: str) -> str:
    if " > " in name:
        return name.split(" > ", 1)[1]
    return name


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
    ax.set_xticklabels(matrix.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(n_pathways))
    ax.set_yticklabels(matrix.index, fontsize=8)

    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)

    # Light gridlines
    ax.set_xticks([x - 0.5 for x in range(1, n_conditions)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n_pathways)], minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", size=0)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
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

    # Build condition label: short experiment name + timepoint
    def _condition_label(row):
        exp_short = EXPERIMENT_SHORT.get(row["experiment"], row["experiment"])
        tp = row.get("timepoint")
        if pd.notna(tp) and str(tp) != "None":
            return f"{exp_short} {tp}"
        return exp_short

    # Sort key for conditions: experiment order, then timepoint order
    def _condition_sort_key(row):
        exp_idx = EXPERIMENT_ORDER.index(row["experiment"]) if row["experiment"] in EXPERIMENT_ORDER else 99
        tp = row.get("timepoint")
        tp_idx = TIMEPOINT_ORDER.get(str(tp), 50) if pd.notna(tp) else 0
        return (exp_idx, tp_idx)

    subset["condition"] = subset.apply(_condition_label, axis=1)
    subset["_sort_key"] = subset.apply(_condition_sort_key, axis=1)
    subset = subset.sort_values("_sort_key")

    # -log10(padj), clamp NaN to 0
    subset["neg_log10_padj"] = subset["padj"].apply(
        lambda p: -np.log10(p) if pd.notna(p) and p > 0 else 0.0
    )

    # Short pathway names
    if "pathway_name" in subset.columns:
        subset["label"] = subset["pathway_name"].fillna(subset["pathway_id"]).apply(_short_pathway_name)
    else:
        subset["label"] = subset["pathway_id"]

    if len(subset) == 0:
        return pd.DataFrame()

    # Preserve condition order from sorted subset
    condition_order = list(dict.fromkeys(subset["condition"]))

    matrix = subset.pivot_table(
        index="label",
        columns="condition",
        values="neg_log10_padj",
        aggfunc="max",
        fill_value=0.0,
    )

    # Reorder columns to match experiment/timepoint order
    matrix = matrix[[c for c in condition_order if c in matrix.columns]]

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
