"""Clustered-heatmap utilities for the cyano N-gene inventory analysis.

Method:
    Binarize the copy-count matrix → cell > 0 becomes 1, else 0.
    Compute pairwise Jaccard distance on rows (strains) and columns (groups).
    Build UPGMA (average-linkage) hierarchical clustering on each axis.
    Reorder the matrix from each axis's leaf ordering.
    Render a heatmap of the original copy-count matrix in the reordered layout,
    with row dendrogram + column dendrogram + row annotation strips.

This module is methodology-first: minimal interface, no plotting framework
beyond matplotlib/seaborn defaults.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist


@dataclass
class ClusteredMatrix:
    matrix: pd.DataFrame              # original copy-count matrix, reordered
    row_linkage: np.ndarray            # scipy linkage matrix for rows
    col_linkage: np.ndarray            # scipy linkage matrix for columns
    row_order: list[str]
    col_order: list[str]


def binarize(matrix: pd.DataFrame) -> pd.DataFrame:
    """Convert copy-count matrix to binary presence/absence (>0 → 1)."""
    return (matrix > 0).astype(int)


def jaccard_linkage(binary_matrix: pd.DataFrame, axis: str = "rows") -> np.ndarray:
    """Return UPGMA (average-linkage) hierarchical clustering on Jaccard distance.

    axis='rows': cluster rows (rows × cols).
    axis='cols': cluster columns (transpose first).
    """
    data = binary_matrix.values if axis == "rows" else binary_matrix.T.values
    # scipy's pdist with metric='jaccard' returns 1 - (intersection / union) for binary input
    distances = pdist(data, metric="jaccard")
    return linkage(distances, method="average")


def cluster_inventory(matrix: pd.DataFrame) -> ClusteredMatrix:
    """Full pipeline: binarize → cluster rows + cols → reorder."""
    binary = binarize(matrix)
    row_lk = jaccard_linkage(binary, axis="rows")
    col_lk = jaccard_linkage(binary, axis="cols")
    row_order = [matrix.index[i] for i in leaves_list(row_lk)]
    col_order = [matrix.columns[i] for i in leaves_list(col_lk)]
    reordered = matrix.loc[row_order, col_order]
    return ClusteredMatrix(
        matrix=reordered, row_linkage=row_lk, col_linkage=col_lk,
        row_order=row_order, col_order=col_order,
    )


def render_heatmap(
    clustered: ClusteredMatrix,
    row_annotations: pd.DataFrame,
    col_labels: pd.Series,
    output_path,
    *,
    col_category: pd.Series | None = None,
    col_category_palette: dict | None = None,
    title: str | None = None,
    copy_cap: int = 3,
    figsize: tuple[float, float] = (20, 11),
) -> None:
    """Render the clustered heatmap.

    Layout (left → right):  row dendro | strain labels | row annot | heatmap
    Layout (top → bottom): col dendro | col-category strip | (row 2 above) | colorbar

    `row_annotations` indexed by strain name, columns = annotation tracks.
    `col_labels` indexed by group_id, value = readable label (e.g., gene_name).
    `col_category` (optional) indexed by group_id, value = pathway category string;
        each category drawn as a colored strip above the heatmap.
    `col_category_palette` (optional) dict[category]→hex color.
    """
    matrix = clustered.matrix
    n_rows, n_cols = matrix.shape

    cmap = ListedColormap(["#ffffff", "#9ecae1", "#3182bd", "#08306b"])
    bounds = [-0.5, 0.5, 1.5, 2.5, copy_cap + 0.5]
    norm = BoundaryNorm(bounds, cmap.N)
    display_values = matrix.clip(upper=copy_cap)

    n_row_tracks = row_annotations.shape[1]
    show_col_strip = col_category is not None

    fig = plt.figure(figsize=figsize)
    # 4 columns: row dendro | strain labels | row annot tracks | heatmap
    # 4 rows: col dendro | col-category strip | main | colorbar
    gs = fig.add_gridspec(
        nrows=4, ncols=4,
        width_ratios=[1.0, 1.6, 0.6 * max(n_row_tracks, 1), 6.0],
        height_ratios=[0.7, 0.18, 5.5, 0.12],
        hspace=0.55, wspace=0.05,
    )

    # Column dendrogram (top, above heatmap area)
    ax_coldendro = fig.add_subplot(gs[0, 3])
    dendrogram(clustered.col_linkage, ax=ax_coldendro, no_labels=True, color_threshold=0,
               above_threshold_color="black")
    ax_coldendro.set_xticks([]); ax_coldendro.set_yticks([])
    for spine in ax_coldendro.spines.values():
        spine.set_visible(False)

    # Column category strip (above heatmap, below col dendro)
    if show_col_strip:
        ax_colstrip = fig.add_subplot(gs[1, 3])
        cat_values = [col_category.get(c, "other") for c in clustered.col_order]
        unique_cats = list(dict.fromkeys(cat_values))
        if col_category_palette is None:
            base_cmap = plt.get_cmap("tab10", len(unique_cats))
            col_category_palette = {c: base_cmap(i) for i, c in enumerate(unique_cats)}
        # Draw as rectangles
        for j, cat in enumerate(cat_values):
            color = col_category_palette.get(cat, "#cccccc")
            ax_colstrip.add_patch(plt.Rectangle((j, 0), 1, 1, facecolor=color, edgecolor="none"))
        ax_colstrip.set_xlim(0, n_cols); ax_colstrip.set_ylim(0, 1)
        ax_colstrip.set_xticks([]); ax_colstrip.set_yticks([])
        for spine in ax_colstrip.spines.values():
            spine.set_visible(False)
        # Legend on the right of the strip
        from matplotlib.patches import Patch
        legend_handles = [Patch(facecolor=col_category_palette[c], label=c) for c in unique_cats]
        ax_colstrip.legend(
            handles=legend_handles, loc="center left", bbox_to_anchor=(1.005, 0.5),
            fontsize=7, frameon=False, handlelength=1.0, ncol=1,
        )

    # Row dendrogram (left)
    ax_rowdendro = fig.add_subplot(gs[2, 0])
    dendrogram(clustered.row_linkage, ax=ax_rowdendro, orientation="left",
               no_labels=True, color_threshold=0, above_threshold_color="black")
    ax_rowdendro.set_xticks([]); ax_rowdendro.set_yticks([])
    for spine in ax_rowdendro.spines.values():
        spine.set_visible(False)
    ax_rowdendro.invert_yaxis()

    # Strain labels (column 1 — left of annotation tracks and heatmap)
    ax_strain = fig.add_subplot(gs[2, 1])
    ax_strain.set_xlim(0, 1); ax_strain.set_ylim(n_rows, 0)
    for r, strain in enumerate(clustered.row_order):
        ax_strain.text(0.98, r + 0.5, strain, ha="right", va="center", fontsize=9)
    ax_strain.set_xticks([]); ax_strain.set_yticks([])
    for spine in ax_strain.spines.values():
        spine.set_visible(False)

    # Row annotation tracks
    ax_annot = fig.add_subplot(gs[2, 2])
    ann = row_annotations.loc[clustered.row_order].copy().astype(str).fillna("—")
    n_tracks = ann.shape[1]
    for t, track in enumerate(ann.columns):
        values = ann[track].tolist()
        unique = sorted(set(values))
        codes = np.array([unique.index(v) for v in values])
        palette_cmap = plt.get_cmap("tab20", max(len(unique), 1))
        ax_annot.imshow(
            codes.reshape(-1, 1), aspect="auto", cmap=palette_cmap,
            extent=(t, t + 1, n_rows, 0),
        )
        for r, v in enumerate(values):
            ax_annot.text(t + 0.5, r + 0.5, v, ha="center", va="center", fontsize=7, color="black")
    ax_annot.set_xlim(0, n_tracks); ax_annot.set_ylim(n_rows, 0)
    ax_annot.set_xticks(np.arange(n_tracks) + 0.5)
    ax_annot.set_xticklabels(ann.columns, rotation=90, fontsize=8)
    ax_annot.set_yticks([])
    for spine in ax_annot.spines.values():
        spine.set_visible(False)

    # Heatmap
    ax_main = fig.add_subplot(gs[2, 3])
    im = ax_main.imshow(display_values.values, aspect="auto", cmap=cmap, norm=norm,
                        interpolation="nearest")
    ax_main.set_xticks(np.arange(n_cols))
    ax_main.set_xticklabels(
        [col_labels.get(c, c) for c in clustered.col_order],
        rotation=90, fontsize=7,
    )
    ax_main.set_yticks([])  # strain labels are on the dedicated left axis now
    raw = matrix.values
    for i in range(n_rows):
        for j in range(n_cols):
            v = raw[i, j]
            if v >= 2:
                ax_main.text(j, i, str(int(v)), ha="center", va="center",
                             fontsize=6, color="white")

    # Horizontal colorbar below heatmap
    ax_cbar = fig.add_subplot(gs[3, 3])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation="horizontal",
                        ticks=[0, 1, 2, copy_cap])
    cbar.ax.set_xticklabels(["0", "1", "2", f"≥{copy_cap}"])
    cbar.set_label("copy count", fontsize=9)

    if title:
        fig.suptitle(title, y=0.995, fontsize=11)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------- toy-data verification (run when invoked as a script) ----------
def _toy_verify() -> None:
    """Hand-computable toy: 3 strains × 4 genes.

    s1 = {g1, g2, g3}        binary 1,1,1,0
    s2 = {g1, g2}            binary 1,1,0,0
    s3 = {g3, g4}            binary 0,0,1,1

    Hand-computed Jaccard distances (1 - intersection/union):
      d(s1, s2) = 1 - 2/3 = 1/3 ≈ 0.333
      d(s1, s3) = 1 - 1/4 = 3/4 = 0.75
      d(s2, s3) = 1 - 0/4 = 1.0

    Expected: s1+s2 cluster first (closest), then s3 joins.
    """
    matrix = pd.DataFrame(
        [[1, 1, 1, 0],
         [1, 1, 0, 0],
         [0, 0, 1, 1]],
        index=["s1", "s2", "s3"],
        columns=["g1", "g2", "g3", "g4"],
    )
    binary = binarize(matrix)
    distances = pdist(binary.values, metric="jaccard")
    expected = np.array([1/3, 3/4, 1.0])
    print("Hand-computed expected pdist:", expected)
    print("scipy pdist (Jaccard):       ", distances)
    assert np.allclose(distances, expected), "Jaccard distance mismatch on toy data"

    clustered = cluster_inventory(matrix)
    print("Row order:", clustered.row_order)
    print("Col order:", clustered.col_order)
    # Expected row order: s1, s2 adjacent (they merge first); s3 farther.
    assert clustered.row_order.index("s1") - clustered.row_order.index("s2") in (-1, 1), \
        "s1 and s2 should be adjacent in the clustered ordering"
    print("Toy verification: PASS")


if __name__ == "__main__":
    _toy_verify()
