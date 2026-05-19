"""Build collapsed-by-pathway views of the inventory: Figures 2, 3, 4.

Reads the 19 × 61 inventory matrix and the pathway category mapping,
collapses 61 ortholog groups into 10 pathway categories, and produces:

  Figure 2 — Strain × Pathway heatmap (19 × 10), cell = % of pathway groups
             present in that strain. Strains ordered by Figure 1's clustering.
  Figure 3 — Clade-group × Pathway summary (8 clade groups × 10 pathways),
             cell = mean % across strains in the clade group.
  Figure 4 — Pathway composition stacked bars per strain (horizontal),
             strains ordered by total N-group count descending.

Outputs:
    data/02_strain_pathway_pct.csv         (19 × 10 percent-present matrix)
    data/02_strain_pathway_count.csv       (19 × 10 raw-count matrix)
    data/03_cladegroup_pathway_pct.csv     (8 × 10 mean-percent matrix)
    figures/02_strain_x_pathway.png/.pdf
    figures/03_cladegroup_x_pathway.png/.pdf
    figures/04_strain_pathway_composition.png/.pdf
    data/02_pathway_views.log
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

HERE = Path(__file__).resolve()
ANALYSIS_ROOT = HERE.parents[2]
sys.path.insert(0, str(ANALYSIS_ROOT / "4_methods"))
from n_pathway_categories import CATEGORY_ORDER, CATEGORY_PALETTE, categorize  # noqa: E402

STEP_DIR = HERE.parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"
INVENTORY_DIR = ANALYSIS_ROOT / "2_kg_selection" / "data"
LOG_PATH = DATA_DIR / "02_pathway_views.log"


# Short strain labels (kept consistent with Figure 1)
SHORT = {
    "Prochlorococcus AS9601": "AS9601",
    "Prochlorococcus MED4": "MED4",
    "Prochlorococcus MIT0801": "MIT0801",
    "Prochlorococcus MIT9301": "MIT9301",
    "Prochlorococcus MIT9303": "MIT9303",
    "Prochlorococcus MIT9312": "MIT9312",
    "Prochlorococcus MIT9313": "MIT9313",
    "Prochlorococcus NATL1A": "NATL1A",
    "Prochlorococcus NATL2A": "NATL2A",
    "Prochlorococcus RSP50": "RSP50",
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)": "SS120",
    "Synechococcus WH8102": "WH8102",
    "Synechococcus WH7803": "WH7803",
    "Synechococcus CC9311": "CC9311",
    "Synechococcus sp. BL107": "BL107",
    "Synechococcus PCC 7002": "PCC 7002",
    "Synechococcus elongatus PCC 7942": "PCC 7942 (elong.)",
    "Synechococcus elongatus UTEX 2973": "UTEX 2973 (elong.)",
    "Thermosynechococcus vestitus BP-1": "BP-1 (Thermosyn.)",
}

# Clade-group definitions for Figure 3 (rows ordered Pro HL → Pro LL → Syn → other)
CLADE_GROUPS = {
    "Pro HLI (n=2)":        ["Prochlorococcus MED4", "Prochlorococcus RSP50"],
    "Pro HLII (n=3)":       ["Prochlorococcus AS9601", "Prochlorococcus MIT9301",
                             "Prochlorococcus MIT9312"],
    "Pro LLI (n=3)":        ["Prochlorococcus MIT0801", "Prochlorococcus NATL1A",
                             "Prochlorococcus NATL2A"],
    "Pro LLII (n=1, SS120)": ["Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)"],
    "Pro LLIV (n=2)":       ["Prochlorococcus MIT9303", "Prochlorococcus MIT9313"],
    "Marine Syn (n=4)":     ["Synechococcus WH8102", "Synechococcus WH7803",
                             "Synechococcus CC9311", "Synechococcus sp. BL107"],
    "Non-marine Syn (n=3)": ["Synechococcus PCC 7002",
                             "Synechococcus elongatus PCC 7942",
                             "Synechococcus elongatus UTEX 2973"],
    "Thermosyn (n=1, BP-1)": ["Thermosynechococcus vestitus BP-1"],
}


def collapse_to_pathway(matrix: pd.DataFrame, groups: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """Return (count_per_pathway, pct_present_per_pathway, n_groups_per_pathway).

    count_per_pathway: 19 × 10, cell = number of present ortholog groups in that
                      pathway for that strain (paralogs collapsed; cell ≥1 means
                      strain has at least one member of at least one group in
                      that pathway).
    pct_present_per_pathway: 19 × 10, cell = count_per_pathway / n_groups_per_pathway.
    n_groups_per_pathway: 10 values, total ortholog groups in each pathway.
    """
    # group_id → category
    cat_map = groups.set_index("group_id")["consensus_gene_name"].apply(categorize)
    # presence (>0) per cell
    binary = (matrix > 0).astype(int)
    # For each pathway, sum the binary values across columns of that pathway
    per_pathway = pd.DataFrame(0, index=binary.index, columns=CATEGORY_ORDER, dtype=int)
    n_groups_per = pd.Series(0, index=CATEGORY_ORDER, dtype=int)
    for category in CATEGORY_ORDER:
        cols = [gid for gid in binary.columns if cat_map.get(gid, "other") == category]
        n_groups_per[category] = len(cols)
        if cols:
            per_pathway[category] = binary[cols].sum(axis=1)
    pct = per_pathway.div(n_groups_per.replace(0, np.nan), axis=1).fillna(0)
    return per_pathway, pct, n_groups_per


# ---------- Figure 2 — Strain × Pathway heatmap ----------
def render_strain_pathway(pct: pd.DataFrame, count: pd.DataFrame,
                          n_groups_per: pd.Series, strain_order: list[str],
                          output_path: Path, title: str) -> None:
    """Heatmap of % present per (strain, pathway). Strain order = Figure 1 clustering order."""
    pct = pct.loc[strain_order]
    count = count.loc[strain_order]

    fig, ax = plt.subplots(figsize=(10, 7))
    im = ax.imshow(pct.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)

    n_rows, n_cols = pct.shape
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(
        [f"{cat}\n(n={int(n_groups_per[cat])})" for cat in pct.columns],
        rotation=45, ha="right", fontsize=9,
    )
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels([SHORT.get(s, s) for s in pct.index], fontsize=9)

    # Cell text: "k/N" where k = present groups, N = total in pathway
    for i in range(n_rows):
        for j, cat in enumerate(pct.columns):
            n_tot = int(n_groups_per[cat])
            n_pres = int(count.iat[i, j])
            if n_tot == 0:
                continue
            color = "white" if pct.iat[i, j] > 0.55 else "black"
            ax.text(j, i, f"{n_pres}/{n_tot}", ha="center", va="center",
                    fontsize=8, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("fraction of pathway groups present", fontsize=9)
    ax.set_title(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------- Figure 3 — Clade-group × Pathway summary ----------
def render_cladegroup_pathway(pct: pd.DataFrame, n_groups_per: pd.Series,
                              output_path: Path, title: str) -> pd.DataFrame:
    """Mean % per (clade_group, pathway). Returns the underlying clade-group matrix."""
    # Build clade-group mean-percent matrix
    rows = []
    row_labels = []
    for group_label, strains in CLADE_GROUPS.items():
        strains_in_scope = [s for s in strains if s in pct.index]
        if not strains_in_scope:
            continue
        mean_pct = pct.loc[strains_in_scope].mean(axis=0)
        rows.append(mean_pct.values)
        row_labels.append(group_label)
    cg_pct = pd.DataFrame(rows, index=row_labels, columns=pct.columns)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    im = ax.imshow(cg_pct.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)

    n_rows, n_cols = cg_pct.shape
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(
        [f"{cat}\n(n={int(n_groups_per[cat])})" for cat in cg_pct.columns],
        rotation=45, ha="right", fontsize=9,
    )
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(cg_pct.index, fontsize=9)

    for i in range(n_rows):
        for j in range(n_cols):
            v = cg_pct.iat[i, j]
            color = "white" if v > 0.55 else "black"
            ax.text(j, i, f"{int(round(v * 100))}%", ha="center", va="center",
                    fontsize=8, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("mean fraction of pathway groups present", fontsize=9)
    ax.set_title(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return cg_pct


# ---------- Figure 4 — Pathway composition stacked bars per strain ----------
def render_strain_composition(count: pd.DataFrame, output_path: Path, title: str) -> None:
    """Horizontal stacked bars: each strain's count per pathway. Strains ordered by total descending."""
    totals = count.sum(axis=1)
    order = totals.sort_values(ascending=True).index.tolist()  # ascending so largest at top
    count = count.loc[order]
    totals = totals.loc[order]

    fig, ax = plt.subplots(figsize=(11, 8))
    n_strains = len(count)
    bottoms = np.zeros(n_strains)
    ys = np.arange(n_strains)
    for cat in CATEGORY_ORDER:
        values = count[cat].values
        color = CATEGORY_PALETTE.get(cat, "#cccccc")
        ax.barh(ys, values, left=bottoms, color=color, edgecolor="white",
                linewidth=0.5, label=cat)
        # Inline count for slices ≥ 3
        for y, v, b in zip(ys, values, bottoms):
            if v >= 3:
                ax.text(b + v / 2, y, str(int(v)), ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
        bottoms = bottoms + values

    ax.set_yticks(ys)
    ax.set_yticklabels([SHORT.get(s, s) for s in count.index], fontsize=9)
    # Totals at the right end of each bar
    for y, t in zip(ys, totals):
        ax.text(t + 0.5, y, f" {int(t)}", ha="left", va="center", fontsize=8, color="black")
    ax.set_xlabel("number of present ortholog groups", fontsize=10)
    ax.set_xlim(0, totals.max() + 5)
    ax.set_title(title, fontsize=11)
    ax.legend(loc="lower right", fontsize=8, frameon=False, ncol=2, title="pathway")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
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
    # 01_row_order.txt holds SHORT names; map back to full strain names for indexing.
    short_to_full = {v: k for k, v in SHORT.items()}
    strain_order = [
        short_to_full[s] for s in (DATA_DIR / "01_row_order.txt").read_text().strip().splitlines()
    ]

    count, pct, n_groups_per = collapse_to_pathway(matrix, groups)
    logging.info("Pathway-collapse complete. Pathway → n_groups: %s", dict(n_groups_per))
    logging.info("Per-strain total present groups (matches Figure 1):")
    for strain, n in count.sum(axis=1).items():
        logging.info("  %-55s %3d", strain, n)

    count.to_csv(DATA_DIR / "02_strain_pathway_count.csv")
    pct.to_csv(DATA_DIR / "02_strain_pathway_pct.csv")

    # Figure 2
    out2 = FIG_DIR / "02_strain_x_pathway.png"
    out2_pdf = FIG_DIR / "02_strain_x_pathway.pdf"
    title2 = "N-pathway coverage per cyano strain (19 × 10)\nCells: # present groups / total groups in pathway. Strain order from Figure 1 clustering."
    render_strain_pathway(pct, count, n_groups_per, strain_order, out2, title2)
    render_strain_pathway(pct, count, n_groups_per, strain_order, out2_pdf, title2)

    # Figure 3
    out3 = FIG_DIR / "03_cladegroup_x_pathway.png"
    out3_pdf = FIG_DIR / "03_cladegroup_x_pathway.pdf"
    title3 = "N-pathway coverage by clade group (mean across strains in each group)"
    cg_pct = render_cladegroup_pathway(pct, n_groups_per, out3, title3)
    render_cladegroup_pathway(pct, n_groups_per, out3_pdf, title3)
    cg_pct.to_csv(DATA_DIR / "03_cladegroup_pathway_pct.csv")
    logging.info("Clade-group × pathway summary (mean %% present):")
    for r in cg_pct.index:
        logging.info("  %-30s %s", r, dict(zip(cg_pct.columns, [round(v, 2) for v in cg_pct.loc[r]])))

    # Figure 4
    out4 = FIG_DIR / "04_strain_pathway_composition.png"
    out4_pdf = FIG_DIR / "04_strain_pathway_composition.pdf"
    title4 = "N-machinery composition per cyano strain (ortholog groups, stacked by pathway)"
    render_strain_composition(count, out4, title4)
    render_strain_composition(count, out4_pdf, title4)

    logging.info("Wrote 3 figures (PNG + PDF) and 3 data CSVs.")


if __name__ == "__main__":
    main()
