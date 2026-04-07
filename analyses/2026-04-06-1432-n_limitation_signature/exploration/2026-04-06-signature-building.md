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

## Scoring results (core signature, all conditions)

| Study | Platform | Condition | Timepoint | Hit rate | Reversed | Mean signed log2FC | Rank score | Concordant | Reversed | Not sig |
|-------|----------|-----------|-----------|----------|----------|--------------------|------------|------------|----------|---------|
| Weissberg | rnaseq | axenic | single | 0.818 | 0.015 | 2.680 | 0.743 | 162 | 3 | 33 |
| Weissberg | rnaseq | coculture | day 18 | 0.015 | 0.030 | -0.077 | -0.015 | 3 | 6 | 189 |
| Weissberg | rnaseq | coculture | day 31 | 0.182 | 0.121 | 0.059 | 0.058 | 36 | 24 | 138 |
| Weissberg | rnaseq | coculture | day 60 | 0.177 | 0.106 | 0.139 | 0.068 | 35 | 21 | 142 |
| Weissberg | rnaseq | coculture | day 89 | 0.136 | 0.071 | 0.028 | 0.061 | 27 | 14 | 157 |
| Weissberg | proteomics | axenic | day 14 | 0.121 | 0.005 | 0.415 | 0.113 | 24 | 1 | 173 |
| Weissberg | proteomics | axenic | day 31 | 0.187 | 0.182 | 0.313 | 0.012 | 37 | 36 | 125 |
| Weissberg | proteomics | axenic | day 89 | 0.202 | 0.167 | 0.228 | 0.037 | 40 | 33 | 125 |
| Weissberg | proteomics | coculture | day 18 | 0.066 | 0.000 | 0.470 | 0.065 | 13 | 0 | 185 |
| Weissberg | proteomics | coculture | day 31 | 0.283 | 0.005 | 0.759 | 0.268 | 56 | 1 | 141 |
| Weissberg | proteomics | coculture | day 60 | 0.152 | 0.040 | 0.414 | 0.107 | 30 | 8 | 160 |
| Weissberg | proteomics | coculture | day 89 | 0.167 | 0.040 | 0.434 | 0.121 | 33 | 8 | 157 |
| Tolonen | microarray | reference | 6h | 0.586 | 0.000 | 1.205 | 0.567 | 116 | 0 | 82 |
| Tolonen | microarray | reference | 12h | 0.778 | 0.000 | 1.660 | 0.745 | 154 | 0 | 44 |
| Tolonen | microarray | reference | 24h | 0.712 | 0.000 | 1.477 | 0.688 | 141 | 0 | 57 |
| Tolonen | microarray | reference | 48h | 0.348 | 0.000 | 1.625 | 0.344 | 69 | 0 | 129 |
| Read | rnaseq | reference | 3h | 0.081 | 0.040 | 0.093 | 0.038 | 16 | 8 | 174 |
| Read | rnaseq | reference | 12h | 0.747 | 0.000 | 2.435 | 0.700 | 148 | 0 | 50 |
| Read | rnaseq | reference | 24h | 0.859 | 0.000 | 2.525 | 0.792 | 170 | 0 | 28 |

## Next
Score core and extended signatures against Weissberg 2025 experiments. ← DONE. See `results/signature_scores_core.csv`.

Next priority: walkthrough and verification of these results, then RNA/protein discordance investigation.
