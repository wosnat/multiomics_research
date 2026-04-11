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

from enrich_utils.enrichment import signed_enrichment_score

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

# Experiment groups for visual separation (vertical lines between groups)
EXPERIMENT_GROUPS = [
    # (group_label, [experiment_ids])
    ("Reference", [
        "de_ref_tolonen_ndep",
        "de_ref_read_ndep",
    ]),
    ("Neg. control", [
        "de_ctrl_tolonen_cyanate",
        "de_ctrl_tolonen_urea",
        "de_ctrl_aharonovich_coculture",
        "de_ctrl_steglich_high_white_light",
    ]),
    ("RNA-seq axenic", ["de_weissberg_rnaseq_axenic"]),
    ("RNA-seq coculture", ["de_weissberg_rnaseq_coculture"]),
    ("Proteomics axenic", ["de_weissberg_proteomics_axenic"]),
    ("Proteomics coculture", ["de_weissberg_proteomics_coculture"]),
]

# Flat experiment order derived from groups
EXPERIMENT_ORDER = [exp for _, exps in EXPERIMENT_GROUPS for exp in exps]

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


def _find_group_boundaries(labels: list[str], group_key_fn) -> list[int]:
    """Find indices where the group key changes. Returns positions for separator lines."""
    boundaries = []
    prev_key = None
    for i, label in enumerate(labels):
        key = group_key_fn(label)
        if prev_key is not None and key != prev_key:
            boundaries.append(i)
        prev_key = key
    return boundaries


def make_heatmap(
    matrix: pd.DataFrame,
    title: str,
    out_path: Path,
    cmap: str,
    threshold: float,
    log_fn,
    col_group_boundaries: list[int] | None = None,
    row_group_boundaries: list[int] | None = None,
    col_group_labels: list[tuple[str, int, int]] | None = None,
    row_group_labels: list[tuple[str, int, int]] | None = None,
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
    col_group_boundaries:
        Column indices where group separators should be drawn.
    row_group_boundaries:
        Row indices where group separators should be drawn.
    """
    n_pathways, n_conditions = matrix.shape
    log_fn(f"  Plotting {n_pathways} pathways × {n_conditions} conditions → {out_path.name}")

    # Figure size: scale with dimensions, extra left margin for row group labels
    fig_width = max(10, n_conditions * 0.35 + 6)
    fig_height = max(5, n_pathways * 0.35 + 2.5)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    # Leave space on the left for row group labels
    fig.subplots_adjust(left=0.40)

    # Cap color scale at vmax for readability (extreme values compress the range)
    vmax = min(matrix.values.max(), 10.0) if matrix.values.max() > 0 else 1.0

    im = ax.imshow(
        np.clip(matrix.values, 0, vmax),
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        interpolation="nearest",
    )

    plt.colorbar(im, ax=ax, label="-log10(padj)", shrink=0.8)

    # Mark significant cells with *
    for i in range(n_pathways):
        for j in range(n_conditions):
            val = matrix.values[i, j]
            if val >= threshold:
                ax.text(j, i, "*", ha="center", va="center",
                        fontsize=8, fontweight="bold", color="black" if val < vmax * 0.7 else "white")

    ax.set_xticks(range(n_conditions))
    ax.set_xticklabels(matrix.columns, rotation=60, ha="right", fontsize=7.5)
    ax.set_yticks(range(n_pathways))
    ax.set_yticklabels(matrix.index, fontsize=8)

    fig.suptitle(title, fontsize=12, fontweight="bold", y=1.02)

    # Light cell gridlines
    ax.set_xticks([x - 0.5 for x in range(1, n_conditions)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n_pathways)], minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", size=0)

    # Group separator lines (thicker, dark)
    if col_group_boundaries:
        for b in col_group_boundaries:
            ax.axvline(x=b - 0.5, color="black", linewidth=1.5)
    if row_group_boundaries:
        for b in row_group_boundaries:
            ax.axhline(y=b - 0.5, color="black", linewidth=1.5)

    # Column group labels (above the heatmap)
    if col_group_labels:
        for label, start, end in col_group_labels:
            mid = (start + end - 1) / 2
            ax.text(mid, -1.5, label, ha="center", va="bottom", fontsize=7,
                    fontweight="bold", clip_on=False)

    # Row group labels with bracket lines (left of the heatmap)
    if row_group_labels:
        for label, start, end in row_group_labels:
            if end - start < 1:
                continue
            mid = (start + end - 1) / 2
            # Bracket line in axes fraction x, data y
            bracket_x = -0.08  # left of the y-axis labels
            ax.plot([bracket_x, bracket_x], [start - 0.4, end - 1 + 0.4],
                    transform=ax.get_yaxis_transform(),
                    color="0.4", linewidth=1.2, clip_on=False)
            # Group label text further left
            ax.text(bracket_x - 0.015, mid, label,
                    transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=6, fontweight="bold",
                    clip_on=False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    # Also save SVG for easy editing
    svg_path = out_path.with_suffix(".svg")
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
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

    # Track pathway_id per label for sorting
    label_to_id = dict(zip(subset["label"], subset["pathway_id"]))

    matrix = subset.pivot_table(
        index="label",
        columns="condition",
        values="neg_log10_padj",
        aggfunc="max",
        fill_value=0.0,
    )

    # Reorder columns to match experiment/timepoint order
    matrix = matrix[[c for c in condition_order if c in matrix.columns]]

    # Sort pathways by pathway_id (alphabetical = grouped by level-0 parent)
    pathway_ids = [label_to_id.get(label, label) for label in matrix.index]
    sort_order = sorted(range(len(pathway_ids)), key=lambda i: pathway_ids[i])
    matrix = matrix.iloc[sort_order]

    return matrix


def plot_signed_heatmap(
    enrich_df: pd.DataFrame,
    sig_pathway_ids: set,
    results_dir: Path,
    log_fn,
    col_groups_fn,
    row_groups_fn,
) -> None:
    """Combined signed heatmap: red = up-enriched, blue = down-enriched."""
    scores_df = signed_enrichment_score(enrich_df)

    # Filter to significant pathways
    scores_df = scores_df[scores_df["pathway_id"].isin(sig_pathway_ids)]

    if len(scores_df) == 0:
        log_fn("  No significant pathways for signed heatmap")
        return

    # Build condition label and sort
    def _cond_label(row):
        exp_short = EXPERIMENT_SHORT.get(row.get("experiment", ""), row.get("experiment", ""))
        tp = row.get("timepoint")
        if pd.notna(tp) and str(tp) != "None":
            return f"{exp_short} {tp}"
        return exp_short

    def _cond_sort(row):
        exp = row.get("experiment", "")
        exp_idx = EXPERIMENT_ORDER.index(exp) if exp in EXPERIMENT_ORDER else 99
        tp = row.get("timepoint")
        tp_idx = TIMEPOINT_ORDER.get(str(tp), 50) if pd.notna(tp) else 0
        return (exp_idx, tp_idx)

    scores_df["condition"] = scores_df.apply(_cond_label, axis=1)
    scores_df["_sort"] = scores_df.apply(_cond_sort, axis=1)
    scores_df = scores_df.sort_values("_sort")

    scores_df["label"] = scores_df["pathway_name"].apply(_short_pathway_name)

    condition_order = list(dict.fromkeys(scores_df["condition"]))
    label_to_id = dict(zip(scores_df["label"], scores_df["pathway_id"]))

    matrix = scores_df.pivot_table(
        index="label", columns="condition", values="score",
        aggfunc="first", fill_value=0.0,
    )
    matrix = matrix[[c for c in condition_order if c in matrix.columns]]

    # Sort rows by pathway_id
    pathway_ids = [label_to_id.get(l, l) for l in matrix.index]
    sort_order = sorted(range(len(pathway_ids)), key=lambda i: pathway_ids[i])
    matrix = matrix.iloc[sort_order]

    n_pathways, n_conditions = matrix.shape
    log_fn(f"  Plotting signed heatmap: {n_pathways} pathways × {n_conditions} conditions")

    fig_width = max(10, n_conditions * 0.35 + 6)
    fig_height = max(5, n_pathways * 0.35 + 2.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.40)

    vmax = 10.0
    im = ax.imshow(
        np.clip(matrix.values, -vmax, vmax),
        aspect="auto",
        cmap="RdBu_r",  # red = positive (up), blue = negative (down)
        vmin=-vmax, vmax=vmax,
        interpolation="nearest",
    )

    plt.colorbar(im, ax=ax, label="Signed score: +up / −down", shrink=0.8)

    # Mark significant cells
    for i in range(n_pathways):
        for j in range(n_conditions):
            val = matrix.values[i, j]
            if abs(val) >= NEG_LOG10_THRESHOLD:
                text_color = "white" if abs(val) > vmax * 0.6 else "black"
                ax.text(j, i, "*", ha="center", va="center",
                        fontsize=8, fontweight="bold", color=text_color)

    ax.set_xticks(range(n_conditions))
    ax.set_xticklabels(matrix.columns, rotation=60, ha="right", fontsize=7.5)
    ax.set_yticks(range(n_pathways))
    ax.set_yticklabels(matrix.index, fontsize=8)

    fig.suptitle("Pathway Enrichment — Signed Score", fontsize=12, fontweight="bold", y=1.02)

    # Cell gridlines
    ax.set_xticks([x - 0.5 for x in range(1, n_conditions)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n_pathways)], minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", size=0)

    # Group separators
    col_bounds, col_labels = col_groups_fn(matrix)
    row_bounds, row_labels = row_groups_fn(matrix, enrich_df)

    for b in col_bounds:
        ax.axvline(x=b - 0.5, color="black", linewidth=1.5)
    for b in row_bounds:
        ax.axhline(y=b - 0.5, color="black", linewidth=1.5)

    # Column group labels
    for label, start, end in col_labels:
        mid = (start + end - 1) / 2
        ax.text(mid, -1.5, label, ha="center", va="bottom", fontsize=7,
                fontweight="bold", clip_on=False)

    # Row group bracket labels
    for label, start, end in row_labels:
        if end - start < 1:
            continue
        mid = (start + end - 1) / 2
        bracket_x = -0.08
        ax.plot([bracket_x, bracket_x], [start - 0.4, end - 1 + 0.4],
                transform=ax.get_yaxis_transform(),
                color="0.4", linewidth=1.2, clip_on=False)
        ax.text(bracket_x - 0.015, mid, label,
                transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=6, fontweight="bold",
                clip_on=False)

    fig.tight_layout()
    out = results_dir / "heatmap_enrichment_signed.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), format="svg", bbox_inches="tight")
    plt.close(fig)
    log_fn(f"  Saved {out.name} + .svg")


def plot_discordance_scatter(enrich_df: pd.DataFrame, results_dir: Path, log_fn) -> None:
    """Scatter: RNA-seq vs Proteomics signed enrichment score per pathway.

    Each point = one pathway. ● axenic, ▲ coculture.
    Signed score: positive = up-enriched, negative = down-enriched.
    Points on the diagonal = concordant; off-diagonal = discordance.
    """
    scores = signed_enrichment_score(enrich_df)

    pairs = [
        ("de_weissberg_rnaseq_axenic", "de_weissberg_proteomics_axenic", "Axenic", "o"),
        ("de_weissberg_rnaseq_coculture", "de_weissberg_proteomics_coculture", "Coculture", "^"),
    ]

    fig, ax = plt.subplots(figsize=(7, 6))

    for rna_exp, prot_exp, cond_label, marker in pairs:
        # Max absolute score across timepoints, preserving sign
        rna = scores[scores["experiment"] == rna_exp].copy()
        rna_best = rna.loc[rna.groupby("pathway_id")["score"].apply(lambda s: s.abs().idxmax())]
        rna_best = rna_best.set_index("pathway_id")

        prot = scores[scores["experiment"] == prot_exp].copy()
        prot_best = prot.loc[prot.groupby("pathway_id")["score"].apply(lambda s: s.abs().idxmax())]
        prot_best = prot_best.set_index("pathway_id")

        common = set(rna_best.index) & set(prot_best.index)
        if not common:
            continue

        x_vals, y_vals, labels = [], [], []
        for pid in sorted(common):
            x_vals.append(rna_best.loc[pid, "score"])
            y_vals.append(prot_best.loc[pid, "score"])
            labels.append(_short_pathway_name(str(rna_best.loc[pid, "pathway_name"])))

        alpha = 0.8 if cond_label == "Coculture" else 0.5
        ax.scatter(x_vals, y_vals, alpha=alpha, s=50,
                   marker=marker, edgecolors="gray", linewidths=0.5,
                   label=cond_label, c=["firebrick" if y > 0 else "steelblue" if y < 0 else "gray"
                                        for y in y_vals])

        # Label points with significant signal on either axis
        for x, y, label in zip(x_vals, y_vals, labels):
            if abs(x) > NEG_LOG10_THRESHOLD or abs(y) > NEG_LOG10_THRESHOLD:
                ax.annotate(label, (x, y), fontsize=5.5, ha="left", va="bottom",
                            xytext=(4, 4), textcoords="offset points")

    # Diagonal and axis lines
    lim = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]),
              abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1]), 5) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", alpha=0.3, linewidth=0.8)
    ax.axhline(0, color="gray", linewidth=0.5, alpha=0.3)
    ax.axvline(0, color="gray", linewidth=0.5, alpha=0.3)

    ax.set_xlabel("RNA-seq signed score", fontsize=10)
    ax.set_ylabel("Proteomics signed score", fontsize=10)
    ax.set_title("RNA/Protein Discordance — Signed Enrichment Score",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="lower right")

    fig.tight_layout()
    out = results_dir / "discordance_scatter.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), format="svg", bbox_inches="tight")
    plt.close(fig)
    log_fn(f"  Saved {out.name} + .svg")


def plot_proteomics_trajectory(enrich_df: pd.DataFrame, results_dir: Path, log_fn) -> None:
    """Line plot: signed enrichment score over time for top pathways in proteomics coculture."""
    scores = signed_enrichment_score(enrich_df)

    prot_co = scores[
        (scores["experiment"] == "de_weissberg_proteomics_coculture") &
        (scores["timepoint"].notna()) &
        (scores["timepoint"] != "days 60+89")
    ].copy()

    if len(prot_co) == 0:
        log_fn("  No proteomics coculture data for trajectory plot")
        return

    # Top pathways: non-zero score in at least 2 timepoints
    nonzero_counts = prot_co[prot_co["score"] != 0].groupby("pathway_id").size()
    top_pathways = nonzero_counts[nonzero_counts >= 2].index.tolist()

    if not top_pathways:
        log_fn("  No pathways with signal in >=2 timepoints")
        return

    tp_order = ["day 18", "day 31", "day 60", "day 89"]
    tp_x = {tp: i for i, tp in enumerate(tp_order)}

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = plt.cm.tab10(np.linspace(0, 1, len(top_pathways)))

    for j, pid in enumerate(sorted(top_pathways)):
        pw_data = prot_co[
            (prot_co["pathway_id"] == pid) &
            (prot_co["timepoint"].isin(tp_order))
        ].sort_values("timepoint", key=lambda s: s.map(tp_x))

        if len(pw_data) == 0:
            continue

        xs = [tp_x[tp] for tp in pw_data["timepoint"]]
        ys = pw_data["score"].values
        name = _short_pathway_name(str(pw_data["pathway_name"].iloc[0]))

        ax.plot(xs, ys, marker="o", markersize=5, linewidth=1.8,
                color=colors[j], label=name, alpha=0.85)

    ax.axhline(0, color="gray", linewidth=0.8, alpha=0.5)
    ax.axhline(NEG_LOG10_THRESHOLD, color="gray", linestyle=":", alpha=0.3, linewidth=0.8)
    ax.axhline(-NEG_LOG10_THRESHOLD, color="gray", linestyle=":", alpha=0.3, linewidth=0.8)

    ax.set_xticks(range(len(tp_order)))
    ax.set_xticklabels(tp_order, fontsize=9)
    ax.set_ylabel("Signed enrichment score (+up / −down)", fontsize=9)
    ax.set_xlabel("Timepoint", fontsize=9)
    ax.set_title("Proteomics Coculture — Pathway Enrichment Trajectory",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, loc="best", framealpha=0.8)

    fig.tight_layout()
    out = results_dir / "trajectory_proteomics_coculture.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), format="svg", bbox_inches="tight")
    plt.close(fig)
    log_fn(f"  Saved {out.name} + .svg")


def plot_control_separation(enrich_df: pd.DataFrame, results_dir: Path, log_fn) -> None:
    """Bar chart: number of significant pathways per experiment."""
    sig = enrich_df[enrich_df["padj"] < PADJ_THRESHOLD].copy()

    # Count unique significant pathways per experiment (any direction)
    counts = sig.groupby("experiment")["pathway_id"].nunique().reset_index()
    counts.columns = ["experiment", "n_significant_pathways"]

    # Add experiments with zero significant
    all_exps = enrich_df["experiment"].unique()
    for exp in all_exps:
        if exp not in counts["experiment"].values:
            counts = pd.concat([counts, pd.DataFrame([{
                "experiment": exp, "n_significant_pathways": 0
            }])], ignore_index=True)

    # Short names and order
    counts["short_name"] = counts["experiment"].map(EXPERIMENT_SHORT).fillna(counts["experiment"])
    exp_order_map = {e: i for i, e in enumerate(EXPERIMENT_ORDER)}
    counts["sort_key"] = counts["experiment"].map(exp_order_map).fillna(99)
    counts = counts.sort_values("sort_key").reset_index(drop=True)

    # Color by role
    def _role_color(exp):
        if "ref_" in exp:
            return "forestgreen"
        elif "ctrl_" in exp:
            return "salmon"
        elif "axenic" in exp:
            return "darkorange"
        else:
            return "royalblue"

    counts["color"] = counts["experiment"].apply(_role_color)

    fig, ax = plt.subplots(figsize=(10, 4))
    bars = ax.bar(range(len(counts)), counts["n_significant_pathways"],
                  color=counts["color"], edgecolor="white", linewidth=0.5)

    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts["short_name"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Significant pathways (padj < 0.05)", fontsize=9)
    ax.set_title("Control Separation — Pathway Enrichment",
                 fontsize=12, fontweight="bold")

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="forestgreen", label="Reference (+)"),
        Patch(facecolor="salmon", label="Negative control"),
        Patch(facecolor="darkorange", label="Target (axenic)"),
        Patch(facecolor="royalblue", label="Target (coculture)"),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc="upper right")

    # Group separators
    boundaries = []
    prev_group = None
    for i, exp in enumerate(counts["experiment"]):
        group = _col_group_for_exp(exp)
        if prev_group is not None and group != prev_group:
            boundaries.append(i)
        prev_group = group
    for b in boundaries:
        ax.axvline(x=b - 0.5, color="black", linewidth=1, alpha=0.5)

    fig.tight_layout()
    out = results_dir / "control_separation.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), format="svg", bbox_inches="tight")
    plt.close(fig)
    log_fn(f"  Saved {out.name} + .svg")


def _col_group_for_exp(exp):
    """Map experiment ID to group label."""
    for g_label, g_exps in EXPERIMENT_GROUPS:
        if exp in g_exps:
            return g_label
    return "?"


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
    # 3. Compute group boundaries for column separators
    # ------------------------------------------------------------------
    # Build the full condition order to find where experiment groups change
    # CyanoRak level-0 parent names for row group labels
    LEVEL0_NAMES = {
        "A": "Amino acid",
        "B": "Cofactors",
        "C": "Cell envelope",
        "D": "Cell processes",
        "E": "Central metab.",
        "F": "DNA metab.",
        "G": "Energy metab.",
        "H": "Fatty acid",
        "J": "Photosynthesis",
        "K": "Protein synth.",
        "L": "Protein fate",
        "M": "Nucleotides",
        "N": "Regulatory",
        "O": "Signal transd.",
        "P": "Transcription",
        "Q": "Transport",
        "R": "Other",
    }

    def _col_group_for(col):
        """Map a condition column label to its experiment group label."""
        for g_label, g_exps in EXPERIMENT_GROUPS:
            for e in g_exps:
                if col.startswith(EXPERIMENT_SHORT.get(e, e)):
                    return g_label
        return "?"

    def _col_groups(matrix):
        """Compute column group boundaries and labels."""
        boundaries = []
        labels = []  # (label, start_idx, end_idx)
        prev_group = None
        group_start = 0
        for i, col in enumerate(matrix.columns):
            group = _col_group_for(col)
            if prev_group is not None and group != prev_group:
                boundaries.append(i)
                labels.append((prev_group, group_start, i))
                group_start = i
            prev_group = group
        if prev_group is not None:
            labels.append((prev_group, group_start, len(matrix.columns)))
        return boundaries, labels

    def _row_parent(label, enrich_df):
        """Get CyanoRak level-0 parent letter for a pathway label."""
        name_to_id = dict(zip(
            enrich_df["pathway_name"].apply(_short_pathway_name),
            enrich_df["pathway_id"]
        ))
        pid = name_to_id.get(label, label)
        code = pid.replace("cyanorak.role:", "")
        return code.split(".")[0]

    def _row_groups(matrix, enrich_df):
        """Compute row group boundaries and labels by level-0 parent."""
        boundaries = []
        labels = []
        prev_parent = None
        group_start = 0
        for i, label in enumerate(matrix.index):
            parent = _row_parent(label, enrich_df)
            if prev_parent is not None and parent != prev_parent:
                boundaries.append(i)
                parent_name = LEVEL0_NAMES.get(prev_parent, prev_parent)
                labels.append((parent_name, group_start, i))
                group_start = i
            prev_parent = parent
        if prev_parent is not None:
            parent_name = LEVEL0_NAMES.get(prev_parent, prev_parent)
            labels.append((parent_name, group_start, len(matrix.index)))
        return boundaries, labels

    # ------------------------------------------------------------------
    # 4. Plot up direction
    # ------------------------------------------------------------------
    up_matrix = build_matrix(enrich_df, direction="up", sig_pathways=sig_pathway_ids)
    if len(up_matrix) > 0:
        col_bounds, col_labels = _col_groups(up_matrix)
        row_bounds, row_labels = _row_groups(up_matrix, enrich_df)
        make_heatmap(
            matrix=up_matrix,
            title="Pathway Enrichment — Up-regulated",
            out_path=RESULTS_DIR / "heatmap_enrichment_up.png",
            cmap=CMAP_UP,
            threshold=NEG_LOG10_THRESHOLD,
            log_fn=log,
            col_group_boundaries=col_bounds,
            row_group_boundaries=row_bounds,
            col_group_labels=col_labels,
            row_group_labels=row_labels,
        )
    else:
        log("No significant up-regulated pathways to plot")

    # ------------------------------------------------------------------
    # 5. Plot down direction
    # ------------------------------------------------------------------
    down_matrix = build_matrix(enrich_df, direction="down", sig_pathways=sig_pathway_ids)
    if len(down_matrix) > 0:
        col_bounds, col_labels = _col_groups(down_matrix)
        row_bounds, row_labels = _row_groups(down_matrix, enrich_df)
        make_heatmap(
            matrix=down_matrix,
            title="Pathway Enrichment — Down-regulated",
            out_path=RESULTS_DIR / "heatmap_enrichment_down.png",
            cmap=CMAP_DOWN,
            threshold=NEG_LOG10_THRESHOLD,
            log_fn=log,
            col_group_boundaries=col_bounds,
            row_group_boundaries=row_bounds,
            col_group_labels=col_labels,
            row_group_labels=row_labels,
        )
    else:
        log("No significant down-regulated pathways to plot")

    # ------------------------------------------------------------------
    # 6. Combined signed heatmap
    # ------------------------------------------------------------------
    log("Plotting signed enrichment heatmap...")
    plot_signed_heatmap(enrich_df, sig_pathway_ids, RESULTS_DIR, log,
                        _col_groups, _row_groups)

    # ------------------------------------------------------------------
    # 7. RNA vs protein discordance scatter
    # ------------------------------------------------------------------
    log("Plotting RNA vs protein discordance scatter...")
    plot_discordance_scatter(enrich_df, RESULTS_DIR, log)

    # ------------------------------------------------------------------
    # 7. Proteomics coculture trajectory
    # ------------------------------------------------------------------
    log("Plotting proteomics coculture trajectory...")
    plot_proteomics_trajectory(enrich_df, RESULTS_DIR, log)

    # ------------------------------------------------------------------
    # 8. Control separation bar chart
    # ------------------------------------------------------------------
    log("Plotting control separation bar chart...")
    plot_control_separation(enrich_df, RESULTS_DIR, log)

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
