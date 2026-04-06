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
- Many top-ranked genes lack a gene_name (NaN) — annotated via gene_overview + gene_ontology_terms (see below).

## Follow-up: Unnamed gene annotation

47 of 198 core signature genes have no gene_name. Looked up via `gene_overview` and `gene_ontology_terms`.

**Key findings:**

- [KG] Most (16/20 top unnamed) are COG category S ("Function unknown"), genuinely uncharacterized.
- [KG] **DUF1651 family (Pfam PF07864):** 3 genes all upregulated — PMM0684 (rank 11), PMM0819 (rank 9), PMM1134 (rank 14). [interpretation] This uncharacterized domain family may have a specific role in N-stress response; worth investigating conservation across Prochlorococcus strains.
- [KG] PMM0365 (up, rank 6): DsrE/DsrF-like protein (Pfam PF02635, sulfur relay). [interpretation] May indicate N/S metabolism cross-talk under N-stress.
- [KG] PMM1391 (up, rank 17): ribbon helix-helix domain, GO: regulation of DNA-templated transcription, CyanoRak: DNA interactions. [interpretation] A transcriptional regulator upregulated under N-stress — candidate NtcA regulon member.
- [KG] PMM0996 (up, rank 19): DUF3303, CyanoRak: phosphorus adaptation. [interpretation] P/N cross-talk under nutrient limitation.
- [KG] PMM0605 (down, rank 16): DUF760, CyanoRak: trace metal detoxification. [interpretation] Metal homeostasis affected by N-stress.
- [KG] PMM0699 (down, rank 4): membrane protein, CyanoRak: adaptation/detoxification. [KG] Conserved in Prochlorococcus only (ortholog group size 8).

Full annotation table: `results/unnamed_genes_annotated.csv`

## Next
Score core and extended signatures against Weissberg 2025 experiments.
