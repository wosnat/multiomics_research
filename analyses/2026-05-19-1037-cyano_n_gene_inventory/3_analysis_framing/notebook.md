# Step 3 — Analysis framing

## Context

Step 2 produced the 19 strain × 61 ortholog-group inventory matrix (`2_kg_selection/data/05_inventory_matrix.csv`). Step 3 locks how to present and interpret it. Selection is already done (19 cyano strains, 61 groups = 54 eggnog Cyanobacteria + 7 anchor singletons); this step's job is the framing of the cross-strain comparison.

## Locked framing

Researcher chose **A3 + B3 + C2**, with a secondary A1+B1 view deferred to a later figure pass.

### Primary figure: data-driven clustered heatmap

| Axis | Choice | Details |
|---|---|---|
| Row ordering (A3) | Hierarchical clustering of strains | Jaccard distance on the binarized matrix (cell > 0 → 1, else 0); average linkage (UPGMA). Synthetic singleton groups participate in the clustering — they'll naturally push their host strain toward terminal leaves but won't dominate (only 7 of 61 columns). |
| Column ordering (B3) | Hierarchical clustering of ortholog groups | Same distance + linkage: Jaccard on transposed binarized matrix. Groups that co-occur across the same set of strains cluster together. |
| Cell encoding (C2) | Copy count | Integer values from `05_inventory_matrix.csv`. Long-tail handled by capping the colormap at 3 (cells with copy_count ≥ 3 use the same darkest color); annotations on those cells optional. |
| Row annotations | `clade_canonical`, `genus`, Cyanorak `Pigment type`, `gene_count_kg` | All from `02_strain_table.csv`. Pigment is `NaN` for non-Cyanorak strains; display as "—". |
| Column annotations | `consensus_gene_name`, `gene_category` (when consistent), `source` flag (`eggnog` vs `anchor_singleton`) | From `04_ortholog_groups.csv`. |

### Hypothesis (expected outcome)

`[interpretation]` Pre-registering one named expectation, kept minimal per just-in-time formalization rule:

- **(taxa-clustering)** Strains should cluster largely by taxonomy: HL Prochlorococcus together, LL Prochlorococcus together, marine Synechococcus together, non-marine cyano (PCC 7002 + 2× S. elongatus + Thermosynechococcus) together. If clustering deviates strongly from taxonomy (e.g., SS120 clusters with the S. elongatus pair because both lack urea), that's an interpretable finding for step 6 Discussion.

No other predictions registered — we're letting the heatmap speak.

### Negative-result control

If the heatmap shows random / non-taxonomic clustering, that would suggest the N-gene set is too small (61 groups) to recover meaningful structure, or that the orthology grouping is too coarse. Either would be interpretable; neither rules out the analysis.

## Secondary figure (deferred, planned)

A1 + B1 + C2 view: taxonomy-ordered rows (Pro by clade → marine Syn by clade → non-Cyanorak), Cyanorak-category-ordered columns (Regulation → Assimilation → Ammonium → Urea → Cyanate → Nitrate/Nitrite → Mo cofactor → Met biosynthesis → Other). Same C2 cell encoding. Reads "function category by lineage" and makes the urea-loss in SS120 and WH7803 visually obvious as vertical bands of zeros. Plan to add after the primary figure is reviewed.

## Methods needed (step 4)

- A `heatmap_clustering.py` module with: (i) binarize, (ii) compute pairwise Jaccard, (iii) UPGMA linkage, (iv) reorder matrix from a leaf ordering, (v) heatmap render with row+column annotations + dendrograms.
- Toy-data verification: a 3-strain × 4-gene synthetic input with known Jaccard distances, hand-computed, to validate the module before running on real data.

## Decisions

**2026-05-19 — Cluster on binarized data, encode copy-count in cells.** Splits the role: clustering uses presence/absence (the simplest comparable signal); cells display copy count (preserves the paralog story we worked to recover). Avoids the question "what distance metric for paralog-aware clustering" entirely.

**2026-05-19 — Synthetic singletons participate in clustering.** They are real anchor genes in the inventory; excluding them from clustering would reintroduce the "throw out the orphans" mistake at a different layer. They will pull their host strain slightly toward terminal leaves (each singleton group has 1 member, so contributes 1/19 to that strain's distance from any other). 7 out of 61 columns is a small enough fraction not to dominate.

**2026-05-19 — Average linkage (UPGMA) for both axes.** Standard for biological similarity dendrograms; balanced between single (chaining) and complete (compactness) linkage.

**2026-05-19 — Secondary A1+B1 figure deferred.** Built after primary heatmap reviewed.

## Surprises

None — step 3 is short for this analysis. Framing was straightforward given the inventory is already in hand.

## Advance rationale

Framing is locked. Clustering choices (Jaccard binary + UPGMA), encoding (copy count), and annotations (clade, genus, pigment, gene_count) are all explicit. Step 4 can implement the methods module against a toy example, then step 5 runs it on the real inventory.

---

## Decide-gate checklist

- **Outputs produced.** This `notebook.md` only — no scripts, data, or figures (step 3 is framing prose).
- **Results presented.** Framing tables shown inline above (axis choices + annotations + hypothesis).
- **QC gate.** Framing internal consistency:
  - C2 copy-count encoding pairs with binary clustering (avoids paralog-distance ambiguity) → consistent.
  - Synthetic singletons included in clustering matches the step-2 decision to recover them → consistent.
  - A3 + B3 are data-driven (clustering); A1 + B1 deferred to secondary figure → no conflict.
- **Decisions made this step.** Four: cluster-binary + display-count split; singletons included; UPGMA average linkage; secondary figure deferred.
- **Advance rationale.** Framing complete and minimal per just-in-time formalization. Step 4 (methods module) can begin.
