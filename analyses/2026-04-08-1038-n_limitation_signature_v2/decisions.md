# Decision Log

Design decisions with rationale — WHY the analysis was done this way, not what it does (that's in methods.md).

## Study selection

### Data-driven study selection from KG
**Decision:** Query all MED4 experiments from the KG rather than hardcoding Tolonen + Read.
**Rationale:** Part of understanding the method from the ground up. Discovered that only two N-deprivation studies exist for MED4 in the KG, confirming Tolonen + Read are the right (and only) pair.

### Negative control selection
**Decision:** Four negative controls: Tolonen cyanate/urea (N-sufficient), Aharonovich coculture (coculture without N-limitation), Steglich high white light (non-N stress).
**Rationale:** Cyanate/urea test whether the signature is N-deprivation-specific vs general N-metabolism. Aharonovich separates coculture effect from N effect. Light tests non-N stress specificity. Phosphate (34 genes) and iron (112 genes) experiments skipped — too few genes for meaningful scoring.
**Status:** Steglich proved unreliable (only 43/189 genes matched). Cyanate/urea/Aharonovich are clean controls.

## Timepoint exclusions

### Tolonen 0h and 3h excluded
**Decision:** Filter out Tolonen timepoints 0h and 3h before per-gene summarization.
**Rationale:** 0h has 5 up / 8 down (baseline noise); 3h has 6 up / 0 down (barely responding). Including them dilutes per-gene summaries.

### Read 3h excluded
**Decision:** Filter out Read timepoint 3h.
**Rationale:** Read 3h has 101 DE genes with 2:1 up/down bias. Photosynthesis genes (psbX, rbcL) show transient upregulation at 3h that reverses to downregulation by 12h. This caused rbcL (canonical N-limitation down gene) to be classified as discordant due to the tie-breaker. Excluding Read 3h parallels the Tolonen 0h/3h exclusion — both are transient early responses, not sustained N-deprivation. Result: rbcL correctly enters core signature as "down."

## Signature building

### Tie-breaker by best directional rank
**Decision:** When a gene has equal up and down significant timepoints, break the tie by choosing the direction with the better (lower) directional rank.
**Rationale:** v1 used pandas `idxmax` (arbitrary). Rank-based tie-breaking picks the direction with the stronger signal. Affected rbcL classification until Read 3h was excluded.

### Symmetric extended classification
**Decision:** Both study-A-only and study-B-only genes get `_ns` (present but not significant) vs `_absent` (not in dataset) classification.
**Rationale:** v1 only classified A-only genes relative to B. Symmetric classification is necessary when either study has limited scope (Read's filtered_subset).

### Discordant genes saved separately
**Decision:** Save discordant genes (significant in both studies, opposite direction) to a separate file rather than silently dropping them.
**Rationale:** Discordant genes are biologically interesting — their response may depend on stress severity, timing, or platform. Saving them enables future investigation.

## Scoring

### Single primary metric: rank score
**Decision:** Rank score as the sole primary metric, not three metrics (v1 used hit rate, mean signed log2FC, and rank score).
**Rationale:** One understood metric > three black boxes. Rank score is platform-independent. Hit rate reported alongside as interpretive context.

### Normalization by n_significant, not total genes
**Decision:** Normalized rank = 1 - (dir_rank / n_significant_in_direction) instead of 1 - (dir_rank / total_genes).
**Rationale:** With total_genes normalization (e.g., 1849 for RNA-seq), all ranks compressed to ~0.95+, making the score nearly identical to hit rate. Normalizing by n_significant gives meaningful spread: rank 10 of 300 → 0.967, rank 200 of 300 → 0.333. Discovered during interactive review of v2 results.

### Three scoring tiers
**Decision:** Score with top (38 genes, rank ≤ 10), core (189 genes), and extended (531 genes) tiers.
**Rationale:** If all three agree, the signature is robust. If they diverge, the strongest markers behave differently from weaker ones. Top > Core > Extended was observed consistently — the signal concentrates in the best markers.

### Permutation test on all experiments, not just targets
**Decision:** Run permutation test for every experiment (references, controls, targets), not just Weissberg.
**Rationale:** Reference self-scoring p-values (low by construction) serve as a sanity check. Control p-values test whether the negative controls truly score low. All reported together.

## Marker genes

### Six marker genes traced through every step
**Decision:** glnA (PMM0920), cynA (PMM0370), rbcL (PMM0550), atpD (PMM1452), PMM0030, PMM0346.
**Rationale:** Span up/down directions, named/unnamed, RNA-seq-only/proteomics-detectable. PMM0030 (up, rank 1, RNA-seq only) and PMM0346 (down, rank 1, in proteomics) added as edge-case tracers for platform coverage.
