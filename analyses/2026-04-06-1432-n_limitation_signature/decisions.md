# Decision Log

Design decisions made during brainstorming and implementation, with rationale.
These explain WHY the analysis was done this way, not what it does (that's in methods.md).

## Signature building

### Directional rank as primary metric (not log2fc)

**Decision:** Use `rank_up`/`rank_down` (directional rank within experiment x timepoint) as the primary cross-platform metric. log2fc is secondary, valid only within a single experiment.

**Rationale:** log2fc magnitudes are not comparable across platforms (microarray vs RNA-seq vs proteomics) or statistical pipelines (Goldenspike vs Rockhopper vs DESeq2). Rank is platform-independent — a gene at rank_up 3 is the 3rd most upregulated gene regardless of the fold-change scale.

**Status:** Implemented. Both rank_score and mean_signed_log2fc computed; which produces better axenic/coculture separation is an open question flagged for the walkthrough.

### Peak timepoint selected by rank, not log2fc

**Decision:** When summarizing a multi-timepoint study to one row per gene, pick the timepoint with the best (lowest) directional rank, not the highest |log2fc|.

**Rationale:** Consistent with rank-as-primary. A gene can have a large fold-change at an early timepoint (when few genes are DE, so rank is low) and a smaller fold-change at a later timepoint (when many genes are DE, but rank is better).

### Concordant direction required for core signature

**Decision:** A gene must be significantly DE in the SAME direction (both up or both down) in Tolonen AND Read to be in the core signature. Discordant genes (up in one, down in the other) are excluded entirely — not placed in core or extended.

**Rationale:** Discordant genes are ambiguous — the signal isn't reproducible in direction. Including them would add noise. They could be interesting biologically (e.g., a gene whose response depends on the severity or timing of N-stress) but that's a question for Approach C, not A.

### Majority direction for mixed-timepoint genes

**Decision:** When a gene is significant_up at some timepoints and significant_down at others within a single study, use the majority direction (more timepoints in that direction).

**Rationale:** Some genes show transient responses that reverse. The majority direction captures the dominant pattern. Ties are possible in theory (equal up/down timepoints) but pandas `idxmax` breaks them arbitrarily — acceptable for the signature building step.

## Signature scoring

### Mean signed log2FC formula

**Decision:** `score = mean(sign_i * log2FC_i)` where sign_i is +1 for reference-up genes, -1 for reference-down genes.

**Rationale:** Converts all genes to a common "N-limitation direction" so up-regulated N-transport genes and down-regulated photosynthesis genes both contribute positively when the condition is N-limited. A score > 0 means the condition matches the N-limitation pattern; < 0 means reversal.

**Open question:** Whether to use log2fc or normalized rank for this metric. Both are computed. The spec flags this as "resolve empirically during exploration."

### Rank score formula

**Decision:** `score = mean(concordance_i * normalized_rank_i)` where concordance_i = +1 (same direction as reference), -1 (opposite), or 0 (not significant). normalized_rank_i = 1 - (directional_rank / total_genes).

**Rationale:** Cross-platform safe. Non-significant genes contribute 0 (not noise). Top-ranked concordant genes contribute near +1, top-ranked reversed genes near -1.

### Tolonen 0h and 3h excluded

**Decision:** Filter out Tolonen timepoints 0h and 3h before summarizing.

**Rationale:** 0h has 5 up / 8 down (essentially baseline noise); 3h has 6 up / 0 down (barely responding). Including them would dilute the per-gene summary — a gene might get its "peak" assigned to a noisy early timepoint.

### Extended signature scored in parallel

**Decision:** Run all three metrics for both core AND extended signatures, then compare.

**Rationale:** The core signature (intersection) is conservative. Some genes are in "extended" not because they failed to replicate, but because Read 2017's filtered_subset excluded low-expression genes. Scoring both and comparing tests whether the core is representative or whether the extended adds signal. If they tell the same story, the core is sufficient.

### Axenic RNA-seq plotted at day 18

**Decision:** The axenic RNA-seq experiment (labeled `timepoint=single`) is plotted at day 18 on trajectory figures.

**Rationale:** The experiment is a starvation-vs-exponential comparison at ~day 11 in PRO99-lowN. Day 18 is the first coculture timepoint, so plotting there enables visual comparison. The actual experimental day is earlier, but for the trajectory narrative ("what does axenic look like when coculture measurement begins?") day 18 is the right position. This is a visualization choice, not an analytical one — the scores are the same regardless of x-position.

### Reference baselines computed for the signature's own studies

**Decision:** Score the Tolonen and Read experiments (which defined the signature) per-timepoint using the same metrics, as reference baselines.

**Rationale:** Provides an "acute N-deprivation trajectory" for comparison. Not circular because: (a) hit rate is NOT 100% per-timepoint — many signature genes aren't yet significant at early timepoints; (b) log2fc and rank vary by timepoint; (c) the signature was defined across ALL timepoints, but per-timepoint scoring shows the progression.

## Alteromonas check

### Broad search terms for N-recycling genes

**Decision:** Used 7 broad search terms ("ammonium AND transport", "peptidase OR protease", etc.) to find 171 candidate genes, rather than a curated list.

**Rationale:** We don't know which specific N-recycling pathways Alteromonas uses. A broad net captures candidates for exploratory analysis. The 65 significant genes can be curated later for the publication.
