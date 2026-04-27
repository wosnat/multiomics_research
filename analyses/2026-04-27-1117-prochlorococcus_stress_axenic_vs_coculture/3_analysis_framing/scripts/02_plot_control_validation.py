"""Heatmap of control DE values across experiments × TPs.

Inputs:
  - data/control_de_long.csv (from 01_pull_control_de.py)
  - data/control_panel.csv

Outputs:
  - figures/control_validation_heatmap.png
  - figures/control_validation_heatmap.pdf

Run from the multiomics_research repo root:
  uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/3_analysis_framing/scripts/02_plot_control_validation.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ANALYSIS_DIR / "data"
FIG_DIR = ANALYSIS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


AXIS_ORDER = ["n_stress", "photo", "proteotoxic", "oxidative", "cell_death", "negative"]
OMICS_ORDER = ["RNASEQ", "PROTEOMICS"]
CONDITION_ORDER = ["axenic", "coculture"]


def short_condition(bg: str) -> str:
    s = str(bg).lower()
    if "axenic" in s:
        return "axenic"
    if "coculture" in s:
        return "coculture"
    return s


def main() -> None:
    de = pd.read_csv(DATA_DIR / "control_de_long.csv")
    controls = pd.read_csv(DATA_DIR / "control_panel.csv")

    de["condition"] = de["background_factors"].apply(short_condition)
    de["tp_label"] = de["timepoint"].fillna("single")

    # Build column key: (omics, condition, tp_hours_or_label)
    # Sort columns by omics, condition, tp
    de["col_key"] = (
        de["omics_type"].astype(str)
        + " | "
        + de["condition"]
        + " | "
        + de["tp_label"].astype(str)
    )

    # Compute column ordering
    def col_sort_key(row: pd.Series) -> tuple:
        omics = OMICS_ORDER.index(row["omics_type"]) if row["omics_type"] in OMICS_ORDER else 99
        cond = CONDITION_ORDER.index(row["condition"]) if row["condition"] in CONDITION_ORDER else 99
        tp_h = pd.to_numeric(row.get("timepoint_hours"), errors="coerce")
        if pd.isna(tp_h):
            tp_sort = -1
        else:
            tp_sort = float(tp_h)
        return (omics, cond, tp_sort)

    col_meta = (
        de[["col_key", "omics_type", "condition", "timepoint", "tp_label", "timepoint_hours", "growth_phase"]]
        .drop_duplicates(subset=["col_key"])
        .copy()
    )
    col_meta["sort_key"] = col_meta.apply(col_sort_key, axis=1)
    col_meta = col_meta.sort_values("sort_key").reset_index(drop=True)
    col_order = col_meta["col_key"].tolist()

    # Row ordering: by axis (custom order), role (positive before negative), then locus_tag
    controls["axis_rank"] = controls["axis"].map({a: i for i, a in enumerate(AXIS_ORDER)})
    controls["role_rank"] = controls["role"].map({"positive": 0, "negative": 1})
    row_order = controls.sort_values(["axis_rank", "role_rank", "locus_tag"])
    row_keys = row_order["locus_tag"].tolist()
    row_labels = [
        f"{r['gene_name']} ({r['locus_tag']})"
        for _, r in row_order.iterrows()
    ]

    # Build wide log2fc matrix and significance mask
    wide_lfc = de.pivot_table(
        index="locus_tag", columns="col_key", values="log2fc", aggfunc="first"
    ).reindex(index=row_keys, columns=col_order)

    # Significance: any row whose expression_status contains "significant" but not "not_significant"
    de["is_sig"] = (
        de["expression_status"].astype(str).str.contains("significant", case=False, na=False)
        & ~de["expression_status"].astype(str).str.contains("not_significant", case=False, na=False)
    )
    wide_sig = de.pivot_table(
        index="locus_tag", columns="col_key", values="is_sig", aggfunc="first"
    ).reindex(index=row_keys, columns=col_order).astype(bool)

    # Color norm — symmetric around 0
    vmax = float(np.nanmax(np.abs(wide_lfc.values))) if wide_lfc.notna().any().any() else 4.0
    vmax = max(vmax, 1e-3)
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)

    fig_w = max(8, 0.5 * len(col_order) + 5)
    fig_h = max(7, 0.32 * len(row_keys) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(
        wide_lfc.values,
        cmap="RdBu_r",
        norm=norm,
        aspect="auto",
    )

    # Add * markers for significant cells
    for i, locus in enumerate(row_keys):
        for j, col in enumerate(col_order):
            val = wide_lfc.loc[locus, col]
            sig = wide_sig.loc[locus, col]
            if pd.isna(val):
                ax.text(j, i, "·", ha="center", va="center", color="lightgray", fontsize=10)
            elif sig:
                ax.text(j, i, "*", ha="center", va="center", color="black", fontsize=11, fontweight="bold")

    # Column labels: short (omics short + condition + day)
    col_short_labels = []
    for _, r in col_meta.iterrows():
        omics_short = "RNA" if r["omics_type"] == "RNASEQ" else "Prot"
        cond_short = "Ax" if r["condition"] == "axenic" else "Co"
        tp_short = str(r["tp_label"]).replace("day ", "d")
        if r.get("growth_phase") and pd.notna(r.get("growth_phase")) and "death" in str(r["growth_phase"]).lower():
            tp_short += "†"  # daggered = death phase
        col_short_labels.append(f"{omics_short}/{cond_short}\n{tp_short}")

    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_short_labels, rotation=0, fontsize=8)
    ax.set_yticks(range(len(row_keys)))
    ax.set_yticklabels(row_labels, fontsize=8)

    # Axis-group horizontal lines
    cumulative = 0
    for axis in AXIS_ORDER:
        n = (row_order["axis"] == axis).sum()
        if n == 0:
            continue
        if cumulative > 0:
            ax.axhline(cumulative - 0.5, color="black", lw=0.8)
        # Label on the right
        ax.text(
            len(col_order) - 0.4, cumulative + n / 2 - 0.5,
            axis, ha="left", va="center", fontsize=9, fontweight="bold",
            rotation=270,
        )
        cumulative += n

    # Vertical line between RNA and Proteomics
    n_rna = (col_meta["omics_type"] == "RNASEQ").sum()
    if 0 < n_rna < len(col_order):
        ax.axvline(n_rna - 0.5, color="black", lw=1.0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("log2 fold-change\n(starvation vs exponential, within-condition)", fontsize=9)

    title = (
        "Step 3 control panel — DE validation\n"
        "(* = significant; · = not measured; † = death phase)"
    )
    ax.set_title(title, fontsize=11)

    plt.tight_layout()
    out_png = FIG_DIR / "control_validation_heatmap.png"
    out_pdf = FIG_DIR / "control_validation_heatmap.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")
    print(f"Figure size: {len(row_keys)} genes x {len(col_order)} (omics, condition, TP) cells")


if __name__ == "__main__":
    main()
