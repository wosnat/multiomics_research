"""Step 4 exploration figures.

Three figures:
- Fig A: per-TP trajectory (score_up / score_down / final) per category × ontology,
  split by omics. Shows temporal dynamics: coculture early-vs-late, axenic_dying
  d14 magnitudes, axenic_dead decay pattern.
- Fig B: per-term contribution heatmap. Rows = 13 signature terms. Columns = T
  clusters (exp × tp × ontology), grouped by category + omics + ontology. Color =
  max_abs_contribution (diverging blue-red, display cap ±5 per D2 convention).
- Fig C: category-mean asymmetry bars. Per (ontology × omics × direction), bars
  for axenic_dying vs coculture vs axenic_dead. Shows the asymmetric-outcome
  pattern compactly.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = ANALYSIS_DIR / "exploration" / "qc"

DISPLAY_CAP = 5.0
CATEGORY_ORDER = ["axenic_dying", "axenic_dead", "coculture"]
ONTOLOGY_ORDER = ["cyanorak_role", "kegg"]
OMICS_ORDER = ["RNA", "Prot"]

# Category palette (biology-facing semantic: dying cells = orange, dead = grey,
# rescued = green).
CATEGORY_COLOR = {
    "axenic_dying": "#d9722f",
    "axenic_dead": "#8a8a8a",
    "coculture": "#2b8a3e",
}


def _tp_to_days(tp: str) -> float:
    """Convert timepoint string to numeric days for ordering."""
    if not isinstance(tp, str):
        return np.nan
    t = tp.strip().lower()
    if t.startswith("day "):
        try:
            return float(t.replace("day ", "").strip())
        except ValueError:
            pass
    if t.startswith("days "):
        # e.g. "days 60+89" — use mean
        body = t.replace("days ", "").strip()
        try:
            parts = [float(x) for x in body.replace("+", " ").split()]
            return float(np.mean(parts))
        except ValueError:
            pass
    return np.nan


def _tp_short(tp: str) -> str:
    if not isinstance(tp, str):
        return ""
    return (
        tp.replace("day ", "d")
          .replace("days ", "d")
          .replace("+", "+")
    )


# --------------------------------------------------------------------------
# Fig A — per-TP trajectory
# --------------------------------------------------------------------------

def fig_a_trajectory(scores_all: pd.DataFrame) -> None:
    t = scores_all[(scores_all["class"] == "T") & (scores_all["category"] != "")].copy()
    t["tp_days"] = t["timepoint"].apply(_tp_to_days)
    t = t.sort_values(["category", "omics", "ontology", "tp_days"])

    fig, axes = plt.subplots(
        nrows=len(ONTOLOGY_ORDER), ncols=len(OMICS_ORDER),
        figsize=(12, 7), sharey=True,
    )
    for i, ontology in enumerate(ONTOLOGY_ORDER):
        for j, omics in enumerate(OMICS_ORDER):
            ax = axes[i, j]
            sub = t[(t["ontology"] == ontology) & (t["omics"] == omics)]
            for cat in CATEGORY_ORDER:
                csub = sub[sub["category"] == cat].sort_values("tp_days")
                if csub.empty:
                    continue
                color = CATEGORY_COLOR[cat]
                ax.plot(csub["tp_days"], csub["score_up"], "-^", color=color,
                        alpha=0.85, markersize=7, label=f"{cat} ↑")
                ax.plot(csub["tp_days"], csub["score_down"], "-v", color=color,
                        alpha=0.55, markersize=7, linestyle="--",
                        label=f"{cat} ↓")
            ax.axhline(0, color="black", linewidth=0.5)
            ax.axhline(DISPLAY_CAP, color="grey", linewidth=0.3, linestyle=":")
            ax.set_title(f"{ontology} × {omics}", fontsize=10)
            ax.set_xlabel("days")
            if j == 0:
                ax.set_ylabel("dir score (↑ solid, ↓ dashed)")
            ax.grid(alpha=0.25)
            # Consolidated legend only in top-right.
            if i == 0 and j == len(OMICS_ORDER) - 1:
                ax.legend(loc="upper right", fontsize=7, ncol=2, frameon=True)
    fig.suptitle(
        "Fig A — per-TP signature engagement trajectories\n"
        "category colored (dying=orange, dead=grey, coculture=green); "
        "solid=up-direction, dashed=down-direction; dotted grey=|5| cap reference",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "step4_figA_trajectory.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Fig B — per-term contribution heatmap
# --------------------------------------------------------------------------

def fig_b_contribution_heatmap(
    decomp: pd.DataFrame, signature: pd.DataFrame
) -> None:
    # Column = T cluster identifier; row = signature term
    t = decomp.copy()
    t["tp_days"] = t["timepoint"].apply(_tp_to_days)
    t["col_key"] = (
        t["category"].str[:3] + "|"
        + t["omics"] + "|"
        + t["timepoint"].apply(_tp_short) + "|"
        + t["ontology"].str[:3]
    )
    # Order columns: category → omics → tp_days → ontology
    cat_rank = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    omics_rank = {o: i for i, o in enumerate(OMICS_ORDER)}
    ont_rank = {o: i for i, o in enumerate(ONTOLOGY_ORDER)}
    col_order_df = (
        t[["col_key", "category", "omics", "tp_days", "ontology"]]
        .drop_duplicates()
        .assign(
            cat_i=lambda d: d["category"].map(cat_rank),
            om_i=lambda d: d["omics"].map(omics_rank),
            on_i=lambda d: d["ontology"].map(ont_rank),
        )
        .sort_values(["cat_i", "om_i", "tp_days", "on_i"])
    )
    col_order = col_order_df["col_key"].tolist()

    # Row order: cyanorak up, cyanorak down, kegg up, kegg down
    sig = signature.copy()
    sig["ont_i"] = sig["ontology"].map(ont_rank)
    sig["dir_i"] = sig["direction"].map({"up": 0, "down": 1})
    sig = sig.sort_values(["ont_i", "dir_i", "term_id"])
    row_order = sig["term_id"].tolist()
    row_labels = (sig["term_id"].str.split(":", n=1).str[1] + " "
                  + sig["direction"].map({"up": "↑", "down": "↓"})).tolist()

    pivot = t.pivot_table(
        index="term_id", columns="col_key",
        values="max_abs_contribution", aggfunc="first",
    ).reindex(index=row_order, columns=col_order)

    fig, ax = plt.subplots(figsize=(max(12, 0.45 * pivot.shape[1]),
                                    max(5, 0.4 * pivot.shape[0] + 1.5)))
    display = pivot.clip(lower=-DISPLAY_CAP, upper=DISPLAY_CAP)
    im = ax.imshow(display.values, cmap="RdBu_r",
                   vmin=-DISPLAY_CAP, vmax=DISPLAY_CAP, aspect="auto")
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=90, fontsize=7)
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_labels, fontsize=8)

    # Annotate cells: saturation marker + sign-flag for anti-signature.
    for i, term_id in enumerate(row_order):
        for j, col_key in enumerate(col_order):
            v = pivot.iloc[i, j]
            if pd.isna(v) or v == 0:
                continue
            abs_v = abs(v)
            glyph = ""
            if abs_v >= 10:
                glyph = "**"
            elif abs_v >= 5:
                glyph = "*"
            # Anti-signature (negative) in red outline
            if v < 0:
                glyph = "-" + glyph
            if glyph:
                ax.text(j, i, glyph, ha="center", va="center",
                        fontsize=6, color="black")

    # Category dividers (vertical lines)
    prev_cat = None
    for j, k in enumerate(col_order):
        cat = col_order_df.iloc[j]["category"]
        if prev_cat is not None and cat != prev_cat:
            ax.axvline(j - 0.5, color="black", linewidth=1.2)
        prev_cat = cat

    # Ontology divider (horizontal)
    n_cyan = sum(1 for t in row_order if t.startswith("cyanorak"))
    ax.axhline(n_cyan - 0.5, color="black", linewidth=1.2)

    # Direction divider within ontology
    n_cyan_up = sum(1 for _, r in sig[sig["ontology"] == "cyanorak_role"].iterrows() if r["direction"] == "up")
    ax.axhline(n_cyan_up - 0.5, color="grey", linewidth=0.6, linestyle=":")
    n_cyan_plus_kegg_up = n_cyan + sum(1 for _, r in sig[sig["ontology"] == "kegg"].iterrows() if r["direction"] == "up")
    ax.axhline(n_cyan_plus_kegg_up - 0.5, color="grey", linewidth=0.6, linestyle=":")

    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02,
                 label=f"max-abs contribution (cap ±{DISPLAY_CAP:.0f}; saturation star * at |5|, ** at |10|)")
    ax.set_title(
        "Fig B — per-term contribution heatmap\n"
        "rows: signature terms (cyanorak ↑/↓ | kegg ↑/↓); columns: T clusters grouped by category | omics | tp | ont3"
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "step4_figB_contributions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Fig C — category-mean asymmetry bars
# --------------------------------------------------------------------------

def fig_c_category_means(cat_means: pd.DataFrame) -> None:
    # Filter to per-omics rows only (exclude omics='ALL'), direction in {up, down}.
    df = cat_means[
        (cat_means["omics"].isin(OMICS_ORDER))
        & (cat_means["direction"].isin(["up", "down"]))
    ].copy()

    fig, axes = plt.subplots(
        nrows=len(ONTOLOGY_ORDER), ncols=1, figsize=(11, 7),
        sharex=True,
    )
    x_base = np.arange(len(CATEGORY_ORDER))
    bar_width = 0.2
    for i, ontology in enumerate(ONTOLOGY_ORDER):
        ax = axes[i]
        sub = df[df["ontology"] == ontology]
        # Four bars per category: (RNA up, RNA down, Prot up, Prot down)
        for k, (omics, direction, offset, hatch) in enumerate([
            ("RNA", "up", -1.5, ""),
            ("RNA", "down", -0.5, "////"),
            ("Prot", "up", 0.5, ""),
            ("Prot", "down", 1.5, "////"),
        ]):
            vals = []
            for cat in CATEGORY_ORDER:
                row = sub[(sub["category"] == cat) & (sub["omics"] == omics)
                         & (sub["direction"] == direction)]
                vals.append(float(row.iloc[0]["mean"]) if not row.empty else 0.0)
            color = {"RNA": "#4f8ed6", "Prot": "#e07b3c"}[omics]
            ax.bar(
                x_base + offset * bar_width, vals,
                width=bar_width, color=color, hatch=hatch, alpha=0.85,
                edgecolor="black", linewidth=0.5,
                label=f"{omics} {'↑' if direction == 'up' else '↓'}",
            )
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xticks(x_base)
        ax.set_xticklabels(CATEGORY_ORDER)
        ax.set_ylabel("category-mean dir_score")
        ax.set_title(f"{ontology}", fontsize=10)
        ax.grid(axis="y", alpha=0.25)
        if i == 0:
            ax.legend(loc="upper left", fontsize=8, ncol=4, frameon=True)

    fig.suptitle(
        "Fig C — category-mean dir_scores per ontology × omics × direction\n"
        "solid = up-direction; hatched = down-direction; blue = RNA; orange = Prot",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "step4_figC_category_means.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    scores_all = pd.read_csv(ANALYSIS_DIR / "results" / "scores_all.csv")
    cat_means = pd.read_csv(ANALYSIS_DIR / "results" / "category_mean_scores.csv")
    decomp = pd.read_csv(
        ANALYSIS_DIR / "exploration" / "qc" / "step4_t_contribution_decomposition.csv"
    )
    signature = pd.read_csv(ANALYSIS_DIR / "data" / "reference_signature.csv")

    fig_a_trajectory(scores_all)
    fig_b_contribution_heatmap(decomp, signature)
    fig_c_category_means(cat_means)

    for f in ("step4_figA_trajectory.png", "step4_figB_contributions.png",
              "step4_figC_category_means.png"):
        print(f"wrote exploration/qc/{f}")


if __name__ == "__main__":
    main()
