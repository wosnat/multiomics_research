# 2026-04-06: Signature building

## Question
How many MED4 genes form a reproducible N-limitation signature across
Tolonen 2006 (microarray) and Read 2017 (RNA-seq)?
Type: exploratory

## Approach
Extract all DE genes from both studies, summarize per gene (peak timepoint
by best directional rank), intersect with concordant direction requirement.

## Findings

- [KG] Tolonen N-deprivation (6-48h, after filtering 0h/3h): 409 significant genes (196 up, 213 down)
- [KG] Read N-depleted (3-24h): 388 significant genes (177 up, 211 down)
- [KG] Core signature (concordant intersection): 198 genes (83 up, 115 down)
- [KG] Extended signature: 367 genes (173 read_only, 122 tolonen_only_read_ns, 72 tolonen_only_read_absent)
- [KG] Top core genes by rank (locus_tag, gene_name, direction):
  - PMM0030, -, up (rank 1)
  - PMM0370, cynA, up (rank 1)
  - PMM0346, -, down (rank 1)
  - PMM0337, -, up (rank 1)
  - PMM0550, rbcL, down (rank 1)
  - PMM1452, atpD, down (rank 1)
  - PMM1453, atpF, down (rank 1)
  - PMM0690, hli, up (rank 2)
  - PMM1404, hli, up (rank 2)
  - PMM0920, glnA, up (rank 2)

## Assessment

- **Established:** The core signature of 198 genes is well-powered (well above the 30-gene minimum for permutation tests). The signature captures known N-limitation biology: glnA (glutamine synthetase) and cynA (cyanate transporter) are canonical N-limitation markers; hli (high-light inducible) genes are up, consistent with N stress triggering photoacclimation; rbcL and ATP synthase subunits (atpD, atpF, atpG, atpH) are down, reflecting reduced biosynthetic capacity. The down-bias (115 vs 83) matches expectations: N starvation broadly suppresses translation and photosynthesis.
- **Preliminary:** The 122 tolonen_only_read_ns genes (present but not significant in Read) may include real signal lost to dataset coverage differences; worth revisiting for scoring with relaxed thresholds if core scoring is insufficient.

## Gaps and friction

- 72 Tolonen-only genes are absent from Read entirely — these are likely genes below the 50th-percentile expression cutoff applied during Read extraction. They cannot be scored against Read data regardless of threshold.
- Many top-ranked genes lack a gene_name (NaN) — these require locus_tag-based KG lookup to interpret functionally.

## Next
Score core and extended signatures against Weissberg 2025 experiments.
