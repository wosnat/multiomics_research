# Step 4 — Methods

## Context

Step 3 locked the framing: A3 + B3 + C2 (hierarchical clustering on both axes via Jaccard distance over the binarized matrix; copy count in cells). This step produces the ad-hoc Python module that implements the clustering + heatmap rendering. The module is imported by step 5's analysis script.

## What I did

Wrote `4_methods/heatmap_clustering.py` — single module providing:

| Function | Purpose |
|---|---|
| `binarize(matrix)` | Copy-count → presence/absence (>0 → 1). |
| `jaccard_linkage(binary_matrix, axis)` | Pairwise Jaccard distance (via `scipy.spatial.distance.pdist`) + UPGMA average linkage (`scipy.cluster.hierarchy.linkage`). |
| `cluster_inventory(matrix)` | Full pipeline: binarize → cluster rows + cols → return `ClusteredMatrix` dataclass with reordered matrix + both linkage matrices + leaf orders. |
| `render_heatmap(clustered, row_annotations, col_labels, output_path, ...)` | Render the clustered heatmap with row dendrogram, row-annotation strip, column dendrogram, discrete copy-count colormap (capped at 3), and cell-text overlays for cells with copy_count ≥ 2. |

Module is methodology-first per Code Lifecycle Phase 1: minimal API, no scaffolding beyond what's needed, no premature productization.

## Toy-data verification

Ran `uv run python 4_methods/heatmap_clustering.py` which invokes `_toy_verify()`:

3-strain × 4-gene synthetic input:
- `s1 = {g1, g2, g3}` → binary `1, 1, 1, 0`
- `s2 = {g1, g2}` → binary `1, 1, 0, 0`
- `s3 = {g3, g4}` → binary `0, 0, 1, 1`

Hand-computed Jaccard distances (1 − intersection / union):
- d(s1, s2) = 1 − 2/3 = 0.333
- d(s1, s3) = 1 − 1/4 = 0.75
- d(s2, s3) = 1 − 0/4 = 1.0

scipy returned exactly those values. The clustered row order was `['s3', 's1', 's2']` — s1 and s2 adjacent (they merge first at d=0.333), s3 at the outer leaf. Matches expectation. **Verification: PASS.**

## Results

Module is ready for use by step 5's analysis script. No data files produced — this step builds infrastructure, not analysis output.

## Decisions

**2026-05-19 — Discrete colormap with cap at 3.** The longest tail in copy_count is PCC 7002 amt1 × 3. Higher cells would extend the colormap range and compress contrast; capping at 3 (with cell-text overlay for cells ≥ 2) preserves both the binary signal and the paralog signal. If later analyses surface cells > 3, revisit.

**2026-05-19 — Average linkage (UPGMA) confirmed.** Same as the step-3 framing decision; implemented in `jaccard_linkage`.

## Advance rationale

Methods module built, toy-tested, and PASS. Step 5 can call `cluster_inventory()` on the real inventory matrix and `render_heatmap()` for the figure.

---

## Decide-gate checklist

- **Outputs produced.** `4_methods/heatmap_clustering.py` (module). No data, no figures.
- **Results presented.** Toy-data verification output shown inline above; expected matches scipy.
- **QC gate.**
  - Toy Jaccard distances match hand computation → PASS.
  - Toy clustering row order matches expectation (closest pair adjacent) → PASS.
  - Module imports cleanly when invoked as a script → PASS.
- **Decisions made this step.** Two: colormap cap at 3; UPGMA linkage confirmed.
- **Advance rationale.** Module verified on a toy; step 5 can run it on real data.
