"""Step 2 diagnostic heatmap v2 — two panels per ontology, class-grouped.

Per ontology, one figure with two stacked panels:
    top panel    = non-T clusters  (R | PC | CTX | NC with dividers)
    bottom panel = T clusters      (Weissberg; Prot | RNA with divider)

Rows per ontology:
    key pathways (from exploration/key_pathways.csv)
    ∪ discovered pathways: max(n_up, n_down) >= 3 in R clusters at hours > 3
                           (matches Step 3 signature rule)

Color scale: display-capped at ±5 (scoring cap remains ±10 per spec).
Cells at |signed_score| >= 5 get stars:
    *    : 5 <= |s| < 10
    **   : 10 <= |s| < 20
    ***  : |s| >= 20

Cell size is uniform across both panels; panels' widths scale with their
cluster count (short panel doesn't stretch).

Figures:
    exploration/qc/step2_heatmap_cyanorak_role.png
    exploration/qc/step2_heatmap_kegg.png
"""
from __future__ import annotations

import os
import re
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
SIGNATURE_THRESHOLD = 3       # spec §5 Step 3 core rule
TIMEPOINT_CUTOFF = 3.0        # decision #1: exclude hours <= 3 from signature

DISPLAY_CAP = 5.0             # visualization color cap
STAR_TIERS = [
    (5.0, "*"),
    (10.0, "**"),
    (20.0, "***"),
]

AUTHOR_TRUNC = 6              # truncate author last name to this many chars

NON_T_ORDER = ["R", "PC", "CTX", "NC"]
DIRECTION_MARK = {"up": "↑", "down": "↓", "ambiguous": "~"}

OMICS_SHORT = {
    "RNASEQ": "RNA",
    "PROTEOMICS": "Prot",
    "MICROARRAY": "Array",
}

# Cell sizing (inches). Used for uniform cell size across panels.
CELL_W_IN = 0.22
CELL_H_IN = 0.36


def short_author(full: str, max_chars: int = AUTHOR_TRUNC) -> str:
    if not isinstance(full, str) or not full.strip():
        return "?"
    last = full.strip().split()[-1]
    return last[:max_chars]


def short_tp(tp: str) -> str:
    if tp is None or (isinstance(tp, float) and np.isnan(tp)):
        return "—"
    s = str(tp).strip()
    low = s.lower()
    if low == "steady state":
        return "ss"
    if low == "darkness":
        return "dark"
    if low == "stationary":
        return "stat"
    m = re.match(r"^days?\s+(.+)$", low)
    if m:
        return f"{m.group(1).replace(' ', '')}d"
    m = re.match(r"^(\d+)\s*min$", low)
    if m:
        return f"{m.group(1)}m"
    return s


def star_for_score(s: float) -> str:
    abs_s = abs(s) if np.isfinite(s) else 0.0
    mark = ""
    for thresh, star in STAR_TIERS:
        if abs_s >= thresh:
            mark = star
    return mark


def discover_pathways(df: pd.DataFrame, r_exp_ids: set[str]) -> dict[str, list[str]]:
    r = df[df["experiment_id"].isin(r_exp_ids)].copy()
    r = r[r["timepoint_hours"] > TIMEPOINT_CUTOFF]
    sig = r[r["p_adjust"] < 0.05]
    agg = (
        sig.groupby(["ontology", "term_id", "direction"])
           .agg(n_clusters=("cluster", "nunique"))
           .reset_index()
    )
    agg_any = (
        agg.groupby(["ontology", "term_id"])
           .agg(n_max=("n_clusters", "max"))
           .reset_index()
    )
    picks = agg_any[agg_any["n_max"] >= SIGNATURE_THRESHOLD]
    return {
        ont: sorted(sub["term_id"].tolist())
        for ont, sub in picks.groupby("ontology")
    }


def build_label(row: pd.Series, include_omics: bool) -> str:
    author = short_author(row["first_author"])
    tp = short_tp(row["timepoint"])
    if include_omics:
        omics = OMICS_SHORT.get(str(row["omics_type"]).upper(), "?")
        return f"{author}|{omics}|{tp}"
    return f"{author}|{tp}"


def condition_from_exp_id(exp_id: str) -> str:
    """Parse axenic / coculture from the T experiment_id suffix."""
    if not isinstance(exp_id, str):
        return ""
    low = exp_id.lower()
    if low.endswith("_coculture"):
        return "coculture"
    if low.endswith("_axenic"):
        return "axenic"
    return ""


def order_within_class(cm_sub: pd.DataFrame, include_omics: bool) -> pd.DataFrame:
    cm_sub = cm_sub.copy()
    cm_sub["_author"] = cm_sub["first_author"].apply(short_author)
    cm_sub["_tp_sort"] = cm_sub["timepoint_hours"].fillna(1e9)
    cm_sub["_dir_sort"] = cm_sub["direction"].map({"down": 0, "up": 1})
    sort_cols = []
    if include_omics:
        # T panel: group by axenic/coculture (the biological contrast).
        # Within each condition, keep omics+tp ordering.
        cm_sub["_cond"] = cm_sub["experiment_id"].apply(condition_from_exp_id)
        sort_cols += ["_cond", "omics_type"]
    sort_cols += ["_author", "_tp_sort", "_dir_sort"]
    return cm_sub.sort_values(sort_cols)


def panel_layout(
    cm: pd.DataFrame, include_omics: bool
) -> tuple[list[str], list[str], list[tuple[str, int, int]]]:
    """Return (cluster order, xlabels, class/omics-group dividers)."""
    if len(cm) == 0:
        return [], [], []
    cm = cm.copy()
    if not include_omics:
        cm["_cls_order"] = cm["class"].map(
            {c: i for i, c in enumerate(NON_T_ORDER)}
        )
        cm = cm.sort_values("_cls_order")
        frames = []
        groups = []
        start = 0
        for cls in NON_T_ORDER:
            sub = cm[cm["class"] == cls]
            if len(sub) == 0:
                continue
            sub_ord = order_within_class(sub, include_omics=False)
            frames.append(sub_ord)
            groups.append((cls, start, start + len(sub_ord)))
            start += len(sub_ord)
        ordered = pd.concat(frames)
    else:
        # T panel: split by axenic/coculture (the biological contrast).
        ordered = order_within_class(cm, include_omics=True)
        groups = []
        start = 0
        for cond in ordered["_cond"].drop_duplicates():
            sub = ordered[ordered["_cond"] == cond]
            label = cond if cond else "?"
            groups.append((label, start, start + len(sub)))
            start += len(sub)
    col_order = ordered["cluster"].tolist()
    xlabels = [build_label(r, include_omics) for _, r in ordered.iterrows()]
    return col_order, xlabels, groups


def draw_panel(
    ax,
    M: np.ndarray,
    xlabels: list[str],
    groups: list[tuple[str, int, int]],
    n_rows: int,
    row_labels: list[str],
    row_is_key: list[bool],
    row_expected: list[str],
    cmap,
    norm,
    title: str,
) -> plt.cm.ScalarMappable | None:
    if M.shape[1] == 0:
        ax.set_visible(False)
        return None
    im = ax.imshow(np.clip(M, -DISPLAY_CAP, DISPLAY_CAP),
                   cmap=cmap, norm=norm, aspect="auto")
    # Stars on saturated cells (use original uncapped signed_score).
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if not np.isfinite(v):
                continue
            mark = star_for_score(v)
            if mark:
                ax.text(j, i, mark, ha="center", va="center",
                        fontsize=6, color="white" if abs(v) >= 3 else "black")
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(xlabels, rotation=90, fontsize=7)
    ax.set_yticks(range(n_rows))
    ytlabels = []
    for lab, is_key, exp_dir in zip(row_labels, row_is_key, row_expected):
        mark = DIRECTION_MARK.get(exp_dir, "")
        prefix = f"{mark} " if mark else "  "
        ytlabels.append(f"{prefix}{lab}")
    ax.set_yticklabels(ytlabels, fontsize=8)
    for tick, is_key in zip(ax.get_yticklabels(), row_is_key):
        if is_key:
            tick.set_fontweight("bold")
    # Dividers + class/omics labels.
    for gname, s, e in groups:
        if s > 0:
            ax.axvline(s - 0.5, color="black", linewidth=1.0)
        mid = (s + e - 1) / 2
        ax.text(mid, -0.8, gname, ha="center", va="bottom",
                fontsize=11, fontweight="bold")
    ax.set_ylim(n_rows - 0.5, -1.7)
    ax.set_title(title, fontsize=10, pad=4)
    return im


def render_ontology_figure(
    matrix: pd.DataFrame,
    non_t_cm: pd.DataFrame,
    t_cm: pd.DataFrame,
    row_labels: list[str],
    row_is_key: list[bool],
    row_expected: list[str],
    title: str,
    outpath: Path,
) -> None:
    top_cols, top_xlabels, top_groups = panel_layout(non_t_cm, include_omics=False)
    bot_cols, bot_xlabels, bot_groups = panel_layout(t_cm, include_omics=True)

    M_top = matrix.reindex(index=row_labels, columns=top_cols).to_numpy(dtype=float)
    M_bot = matrix.reindex(index=row_labels, columns=bot_cols).to_numpy(dtype=float)
    M_top = np.where(np.isfinite(M_top), M_top, np.nan)
    M_bot = np.where(np.isfinite(M_bot), M_bot, np.nan)

    n_rows = len(row_labels)
    n_top = len(top_cols)
    n_bot = len(bot_cols)

    # Uniform cell sizing: each axis width = n_cols * CELL_W_IN.
    label_pad_in = 5.0        # left margin for y-tick labels
    cbar_pad_in = 1.0         # right margin for colorbar
    panel_h_in = max(3.0, n_rows * CELL_H_IN + 1.2)
    panel_gap_in = 1.6
    title_pad_in = 2.4  # headroom for title + 3-line legend
    xaxis_pad_in = 1.4

    fig_w = label_pad_in + max(n_top, n_bot) * CELL_W_IN + cbar_pad_in
    fig_h = title_pad_in + panel_h_in + panel_gap_in + panel_h_in + xaxis_pad_in

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Top panel (non-T).
    top_left = label_pad_in / fig_w
    top_w = (n_top * CELL_W_IN) / fig_w if n_top else 0
    top_bottom_y = (xaxis_pad_in + panel_h_in + panel_gap_in) / fig_h
    top_h = panel_h_in / fig_h
    ax_top = fig.add_axes([top_left, top_bottom_y, top_w, top_h]) \
        if n_top else None

    # Bottom panel (T).
    bot_left = label_pad_in / fig_w
    bot_w = (n_bot * CELL_W_IN) / fig_w if n_bot else 0
    bot_bottom_y = xaxis_pad_in / fig_h
    bot_h = panel_h_in / fig_h
    ax_bot = fig.add_axes([bot_left, bot_bottom_y, bot_w, bot_h]) \
        if n_bot else None

    # Colorbar on the right.
    cbar_left = (fig_w - cbar_pad_in + 0.1) / fig_w
    cbar_w = 0.02
    cax = fig.add_axes([cbar_left, bot_bottom_y, cbar_w, top_bottom_y + top_h - bot_bottom_y])

    norm = TwoSlopeNorm(vmin=-DISPLAY_CAP, vcenter=0.0, vmax=DISPLAY_CAP)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#f0f0f0")

    im = None
    if ax_top is not None:
        im = draw_panel(
            ax_top, M_top, top_xlabels, top_groups, n_rows,
            row_labels, row_is_key, row_expected, cmap, norm,
            title="Non-T clusters — signature + calibration",
        )
    if ax_bot is not None:
        im2 = draw_panel(
            ax_bot, M_bot, bot_xlabels, bot_groups, n_rows,
            row_labels, row_is_key, row_expected, cmap, norm,
            title="T clusters — Weissberg 2025 (scored)",
        )
        im = im or im2

    if im is not None:
        cbar = fig.colorbar(im, cax=cax, extend="both")
        cbar.set_label(f"signed_score (display ±{int(DISPLAY_CAP)}; *≥5 **≥10 ***≥20)",
                       fontsize=9)

    fig.suptitle(title, fontsize=11, y=0.995)

    # Legend block below the title — multi-line so text stays readable.
    legend_y = 0.972
    legend_text = (
        "Cells: signed_score = sign × −log10(padj).   "
        "Color-cap ±5 (preserves biology-range resolution).   "
        "Blank cell = term not tested.\n"
        "Saturation stars (uncapped |signed_score|):   "
        "* ≥ 5 (padj ≤ 1e-5)      "
        "** ≥ 10 (padj ≤ 1e-10)      "
        "*** ≥ 20 (padj ≤ 1e-20, statistical ceiling).\n"
        "Bold row label = a priori key-pathway panel;  regular = discovered-strong "
        "(n_sig ≥ 3 in signature-eligible R).   "
        "↑↓~ beside row label = expected direction."
    )
    fig.text(
        0.5, legend_y, legend_text,
        ha="center", va="top", fontsize=9, color="#333333",
    )
    fig.savefig(outpath, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")


def main() -> None:
    df = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    key_panel = pd.read_csv(ANALYSIS_DIR / "exploration" / "key_pathways.csv")

    # enrichment_all.csv has omics_type NaN for Weissberg T experiments. Pull
    # from classified instead (which is populated).
    df = df.drop(columns=["omics_type"], errors="ignore")
    meta_cols = ["experiment_id", "class", "first_author", "omics_type"]
    df = df.merge(
        classified[meta_cols].drop_duplicates(subset=["experiment_id"]),
        on="experiment_id", how="left",
    )

    cluster_meta = (
        df[["cluster", "experiment_id", "class", "first_author", "omics_type",
            "timepoint", "timepoint_hours", "direction", "ontology"]]
          .drop_duplicates(subset=["cluster", "ontology"])
          .reset_index(drop=True)
    )

    r_exp_ids = set(classified.loc[classified["class"] == "R", "experiment_id"])
    discovered = discover_pathways(df, r_exp_ids)

    for ontology in ["cyanorak_role", "kegg"]:
        key_rows = key_panel[key_panel["ontology"] == ontology].copy()
        key_ids = list(key_rows["term_id"])
        discovered_ids = discovered.get(ontology, [])
        all_ids = list(dict.fromkeys(key_ids + discovered_ids))

        ont_df = df[df["ontology"] == ontology]
        matrix = ont_df.pivot_table(
            index="term_id", columns="cluster",
            values="signed_score", aggfunc="first",
        )

        names = dict(zip(ont_df["term_id"], ont_df["term_name"]))
        key_expected = dict(
            zip(key_rows["term_id"], key_rows["expected_direction"])
        )
        row_labels = []
        row_is_key = []
        row_expected = []
        for tid in all_ids:
            name = names.get(tid, "")
            short = name.split(" > ")[-1][:48]
            row_labels.append(f"{tid.split(':', 1)[-1]}  {short}")
            row_is_key.append(tid in key_ids)
            row_expected.append(key_expected.get(tid, ""))
        matrix = matrix.reindex(all_ids).copy()
        matrix.index = row_labels

        cm_ont = cluster_meta[cluster_meta["ontology"] == ontology].copy()
        non_t_cm = cm_ont[cm_ont["class"].isin(NON_T_ORDER)]
        t_cm = cm_ont[cm_ont["class"] == "T"]

        outpath = (
            ANALYSIS_DIR / "exploration" / "qc"
            / f"step2_heatmap_{ontology}.png"
        )
        render_ontology_figure(
            matrix=matrix,
            non_t_cm=non_t_cm,
            t_cm=t_cm,
            row_labels=row_labels,
            row_is_key=row_is_key,
            row_expected=row_expected,
            title=(f"Step 2 — {ontology}   bold row = key panel; "
                   f"↑↓~ expected; display-capped ±{int(DISPLAY_CAP)}"),
            outpath=outpath,
        )


if __name__ == "__main__":
    main()
