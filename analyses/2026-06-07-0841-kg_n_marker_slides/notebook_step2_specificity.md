# Step 2 — Is the up-regulation specific to N stress? (statistical test)

## Question
The flowers (step 1) showed these genes go UP under nitrogen stress. But they
also move under other conditions. **Is up-regulation actually enriched under N
relative to other stresses, or do these genes just respond to everything?**

## Method (locked framing)
- **Test:** per gene, Fisher's exact on a 2×2 of
  `(experiment is N vs non-N) × (gene significantly UP vs not-up)`.
- **Unit:** experiment (timepoints collapsed — "up" if significant up in ≥1 timepoint).
- **"Tested":** gene is *detected* (has a DE row) in that experiment.
- **No fold-change magnitude is used** — only the DE/not-DE call — so the test is
  robust to platform dynamic-range differences (RNA-seq vs microarray vs proteomics).
  [feedback: compare rank or DE/not-DE, not FC]
- One-sided (`greater`) = up-regulation enriched under N. BH-FDR within each scope.
- **Scopes (both):** MED4 per-gene panel (cynS, cynD, urtA–E, amt1); and pooled
  across all genomes at the gene-family level (urtA, cynS ortholog groups).
- Script: `scripts/n_specificity_fisher.py` (rerunnable against frozen CSVs in `data/`).

## Results [KG]  (`data/n_specificity_results.csv`, figure `figures/fig8_n_specificity_forest.*`)

| scope | gene/family | up/tested (N) | up/tested (other) | OR | q (BH) |
|---|---|---|---|---|---|
| MED4 | cynS | 4/8 | 1/4 | 3.0 | 0.82 |
| MED4 | cynD | 5/7 | 1/2 | 2.5 | 0.82 |
| MED4 | amt1 | 5/8 | 2/6 | 3.3 | 0.82 |
| MED4 | urtA | 5/8 | 3/11 | 4.4 | 0.82 |
| MED4 | urtB | 4/8 | 1/5 | 4.0 | 0.82 |
| MED4 | urtC | 2/7 | 1/5 | 1.6 | 0.82 |
| MED4 | urtD | 2/7 | 1/2 | 0.4 | 0.92 |
| MED4 | urtE | 3/8 | 1/3 | 1.2 | 0.82 |
| **POOLED** | **urtA family** | **8/11** | **6/27** | **9.3** | **0.011 \*** |
| POOLED | cynS family | 4/8 | 1/9 | 8.0 | 0.11 |

## Interpretation
**Plain English:** Pooled across all genomes, urea-transporter (urtA) up-regulation
is significantly more frequent under nitrogen stress than under other conditions —
up in **8 of 11** N experiments vs only **6 of 27** non-N (odds ratio ≈ 9, q = 0.011).
The cyanate-hydratase (cynS) family points the same way (4/8 vs 1/9, OR 8) but does
**not** reach significance (q = 0.11). Within MED4 alone, **7 of 8** genes have OR > 1
(direction = N-specific) but **none is individually significant** (all q ≈ 0.82).

`[interpretation]` Two things are true at once: (1) the up-regulation *is*
preferentially nitrogen-associated — strongly so when the whole dataset is pooled;
(2) it is **not exclusive** — these genes (urtA especially) also go up under some
other stresses, which is why per-gene MED4 odds ratios are modest. The marker signal
is real but **partial**: a panel read together (and pooled across strains) detects N
stress; no single gene in a single genome is a clean, statistically-proven N-only
switch given the data now in the KG.

`[interpretation]` The limiting factor is **power, not effect**: most genes are tested
in only a handful of non-N experiments, so per-gene MED4 tests can't resolve a 3–4×
odds ratio. Pooling across genomes (urtA) is what crosses significance.

## Caveats `[KG]` / `[interpretation]`
- **Underpowered per gene/per strain** — non-N tested counts are small (2–11 in MED4);
  consistent OR>1 across 7/8 genes is itself suggestive but each alone is n.s.
- **Pseudoreplication** — timepoints collapsed to experiment; experiments within one
  study are not fully independent. A study-level random effect would be stricter.
- **Per-publication thresholds** — "significant" uses each study's own DE call, not a
  uniform padj. This is the *point* (platform-robust) but means the unit of "up" varies.
- **tested-absent / table_scope** — for `significant_only` / `top_n` tables the gene is
  only recorded when significant, so "tested-but-not-up" is undercounted there (pooled set
  includes 12 significant_only + 7 top_n + 18 filtered_subset rows; MED4 Read 2017 is
  top-50% filtered_subset). This can inflate up-rates for those experiments.
- **"non-N" is a heterogeneous union** (carbon, P, light, viral, coculture, salt, iron,
  darkness) — the test asks "N vs everything else", not "N vs any one specific stress".
- **Direction = up only** (down and n.s. both count as "not-up").

## Decision
Headline for the slide: **"Up-regulation is significantly N-enriched at the family level
(urtA: 8/11 vs 6/27, OR≈9, q=0.011); cynS trends the same; the signal is a partial,
panel-level marker, not a single-gene on/off switch."** Use `fig8` next to the flowers.
