"""Figure 5 — Secondary view: taxonomy-ordered rows × Cyanorak-category columns.

This is the A1+B1+C2 view deferred from step 3: instead of data-driven
clustering, order strains by taxonomy (HLI → HLII → LLI → LLII → LLIV
Pro → marine Syn → non-marine Syn → Thermosyn) and order columns by
pathway category. Reveals the urea-loss in SS120 + WH7803 as a vertical
band, and makes "HLI is the only Pro with cyanate" visible immediately.

Reuses the heatmap-cell + colorstrip rendering style from Figure 1
but skips the dendrograms (rows/cols are pre-ordered).
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle

HERE = Path(__file__).resolve()
ANALYSIS_ROOT = HERE.parents[2]
STEP_DIR = HERE.parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"
INVENTORY_DIR = ANALYSIS_ROOT / "2_kg_selection" / "data"
LOG_PATH = DATA_DIR / "02_secondary_taxonomy_view.log"

sys.path.insert(0, str(ANALYSIS_ROOT / "4_methods"))
from n_pathway_categories import CATEGORY_ORDER, CATEGORY_PALETTE, categorize  # noqa: E402

# Taxonomic strain order (HLI Pro → HLII Pro → LL Pro → marine Syn → non-marine cyano)
STRAIN_ORDER = [
    # Pro HLI
    "Prochlorococcus MED4", "Prochlorococcus RSP50",
    # Pro HLII
    "Prochlorococcus AS9601", "Prochlorococcus MIT9301", "Prochlorococcus MIT9312",
    # Pro LLI
    "Prochlorococcus MIT0801", "Prochlorococcus NATL1A", "Prochlorococcus NATL2A",
    # Pro LLII (SS120)
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)",
    # Pro LLIV
    "Prochlorococcus MIT9303", "Prochlorococcus MIT9313",
    # Marine Syn
    "Synechococcus WH8102", "Synechococcus WH7803",
    "Synechococcus CC9311", "Synechococcus sp. BL107",
    # Non-marine Syn + Thermosyn
    "Synechococcus PCC 7002",
    "Synechococcus elongatus PCC 7942", "Synechococcus elongatus UTEX 2973",
    "Thermosynechococcus vestitus BP-1",
]

SHORT = {
    "Prochlorococcus MED4": "MED4 (HLI)",
    "Prochlorococcus RSP50": "RSP50 (HLI)",
    "Prochlorococcus AS9601": "AS9601 (HLII)",
    "Prochlorococcus MIT9301": "MIT9301 (HLII)",
    "Prochlorococcus MIT9312": "MIT9312 (HLII)",
    "Prochlorococcus MIT0801": "MIT0801 (LLI)",
    "Prochlorococcus NATL1A": "NATL1A (LLI)",
    "Prochlorococcus NATL2A": "NATL2A (LLI)",
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)": "SS120 (LLII)",
    "Prochlorococcus MIT9303": "MIT9303 (LLIV)",
    "Prochlorococcus MIT9313": "MIT9313 (LLIV)",
    "Synechococcus WH8102": "WH8102 (III)",
    "Synechococcus WH7803": "WH7803 (V)",
    "Synechococcus CC9311": "CC9311 (I)",
    "Synechococcus sp. BL107": "BL107 (IV)",
    "Synechococcus PCC 7002": "PCC 7002",
    "Synechococcus elongatus PCC 7942": "PCC 7942 (elong.)",
    "Synechococcus elongatus UTEX 2973": "UTEX 2973 (elong.)",
    "Thermosynechococcus vestitus BP-1": "BP-1 (Thermosyn.)",
}


def render(matrix: pd.DataFrame, col_label: pd.Series, col_category: pd.Series,
           output_path: Path, title: str, copy_cap: int = 3) -> None:
    n_rows, n_cols = matrix.shape
    cmap = ListedColormap(["#ffffff", "#9ecae1", "#3182bd", "#08306b"])
    bounds = [-0.5, 0.5, 1.5, 2.5, copy_cap + 0.5]
    norm = BoundaryNorm(bounds, cmap.N)
    display_values = matrix.clip(upper=copy_cap)

    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(
        nrows=3, ncols=2,
        width_ratios=[2.0, 6.0],
        height_ratios=[0.4, 5.5, 0.15],
        hspace=0.45, wspace=0.02,
    )

    # Pathway color strip (top, above heatmap)
    ax_strip = fig.add_subplot(gs[0, 1])
    cat_values = [col_category.get(c, "other") for c in matrix.columns]
    for j, cat in enumerate(cat_values):
        ax_strip.add_patch(Rectangle((j, 0), 1, 1,
                                     facecolor=CATEGORY_PALETTE.get(cat, "#cccccc"),
                                     edgecolor="none"))
    # Pathway-category section labels above the strip
    prev = None
    section_start = 0
    for j, cat in enumerate(cat_values + [None]):
        if cat != prev:
            if prev is not None:
                section_end = j
                mid = (section_start + section_end) / 2
                ax_strip.text(mid, 1.5, prev.replace("_", " "), ha="center", va="bottom",
                              fontsize=8, rotation=0)
            section_start = j
            prev = cat
    ax_strip.set_xlim(0, n_cols); ax_strip.set_ylim(0, 1)
    ax_strip.set_xticks([]); ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)
    unique_cats = list(dict.fromkeys(cat_values))
    legend_handles = [Patch(facecolor=CATEGORY_PALETTE[c], label=c) for c in unique_cats]
    ax_strip.legend(handles=legend_handles, loc="center left",
                    bbox_to_anchor=(1.005, 0.5), fontsize=7, frameon=False, ncol=1)

    # Strain labels (left)
    ax_strain = fig.add_subplot(gs[1, 0])
    ax_strain.set_xlim(0, 1); ax_strain.set_ylim(n_rows, 0)
    for r, strain in enumerate(matrix.index):
        ax_strain.text(0.98, r + 0.5, SHORT.get(strain, strain), ha="right", va="center", fontsize=9)
    ax_strain.set_xticks([]); ax_strain.set_yticks([])
    for spine in ax_strain.spines.values():
        spine.set_visible(False)

    # Heatmap
    ax_main = fig.add_subplot(gs[1, 1])
    im = ax_main.imshow(display_values.values, aspect="auto", cmap=cmap, norm=norm,
                        interpolation="nearest")
    ax_main.set_xticks(np.arange(n_cols))
    ax_main.set_xticklabels(
        [col_label.get(c, c) for c in matrix.columns], rotation=90, fontsize=7,
    )
    ax_main.set_yticks([])
    raw = matrix.values
    for i in range(n_rows):
        for j in range(n_cols):
            v = raw[i, j]
            if v >= 2:
                ax_main.text(j, i, str(int(v)), ha="center", va="center",
                             fontsize=6, color="white")
    # Vertical dividers between pathway categories
    prev_cat = None
    for j, cat in enumerate(cat_values):
        if prev_cat is not None and cat != prev_cat:
            ax_main.axvline(j - 0.5, color="black", linewidth=1.2)
        prev_cat = cat
    # Horizontal divider between Pro and Syn
    # Pro count
    n_pro = sum(1 for s in matrix.index if "Prochlorococcus" in s)
    ax_main.axhline(n_pro - 0.5, color="black", linewidth=1.2)

    # Colorbar
    ax_cbar = fig.add_subplot(gs[2, 1])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation="horizontal",
                        ticks=[0, 1, 2, copy_cap])
    cbar.ax.set_xticklabels(["0", "1", "2", f"≥{copy_cap}"])
    cbar.set_label("copy count", fontsize=9)

    fig.suptitle(title, y=0.995, fontsize=11)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    matrix = pd.read_csv(INVENTORY_DIR / "05_inventory_matrix.csv", index_col=0)
    groups = pd.read_csv(INVENTORY_DIR / "04_ortholog_groups.csv")

    # Build col label series (consensus_gene_name with anchor_singleton "*" prefix)
    col_label = groups.set_index("group_id")["consensus_gene_name"].copy()
    for gid in col_label.index:
        v = col_label[gid]
        if pd.isna(v) or v == "":
            col_label[gid] = gid.split(":")[-1][:12]
    for gid in groups[groups["source"] == "anchor_singleton"]["group_id"]:
        if gid in col_label.index:
            col_label[gid] = "*" + str(col_label[gid])

    # Column categories from gene names
    col_category = groups.set_index("group_id")["consensus_gene_name"].apply(categorize)

    # Reorder rows by taxonomy
    matrix_rows = matrix.loc[STRAIN_ORDER]
    # Reorder cols by pathway category (CATEGORY_ORDER), then alphabetically within
    col_order = []
    for cat in CATEGORY_ORDER:
        cat_cols = [gid for gid in matrix.columns if col_category.get(gid, "other") == cat]
        # within category: sort by consensus_gene_name (or group_id)
        cat_cols_sorted = sorted(cat_cols, key=lambda g: (str(col_label.get(g, g)),))
        col_order.extend(cat_cols_sorted)
    matrix_ordered = matrix_rows[col_order]
    logging.info("Taxonomy-ordered matrix shape: %s", matrix_ordered.shape)
    logging.info("Pathway-block sizes (columns): %s",
                 {cat: sum(1 for g in col_order if col_category.get(g, "other") == cat)
                  for cat in CATEGORY_ORDER})

    out_png = FIG_DIR / "05_taxonomy_view.png"
    out_pdf = FIG_DIR / "05_taxonomy_view.pdf"
    title = ("Secondary view — taxonomy-ordered rows × pathway-ordered columns "
             "(A1+B1+C2; horizontal line separates Pro from Syn/Thermosyn; "
             "vertical lines separate pathway categories; * = anchor-singleton group)")
    render(matrix_ordered, col_label, col_category, out_png, title)
    render(matrix_ordered, col_label, col_category, out_pdf, title)
    logging.info("Wrote %s and %s", out_png, out_pdf)


if __name__ == "__main__":
    main()
