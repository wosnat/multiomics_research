# N-Limitation Signature Analysis: Methods

## Research question

Can we quantify the degree of nitrogen limitation in Prochlorococcus MED4
molecularly, and show that coculture with Alteromonas alleviates it over
time, using a gene signature derived from independent reference studies?

## Data scope

### Reference studies
- Tolonen et al. 2006 (DOI: 10.1038/msb4100087) — MED4 N-deprivation time course (microarray, 0-48h)
- Read et al. 2017 (DOI: 10.1038/ismej.2017.88) — MED4 N-depleted vs N-replete (RNA-seq, 3-24h)

### Target study
- Weissberg et al. 2025 (DOI: 10.1101/2025.11.24.690089) — MED4 in PRO99-lowN, axenic vs coculture with Alteromonas HOT1A3 (RNA-seq + proteomics, d14-89)

### Inclusions/exclusions
- Tolonen timepoints 0h and 3h excluded (nearly zero DE genes)
- Read dataset filtered to top 50% of genes by expression level (table_scope: filtered_subset)
- Weissberg d60+89 combined timepoint included in data but excluded from trajectory plots

## Gene selection

Genes were selected based on differential expression (DE) in both Tolonen 2006
and Read 2017, with concordant direction of change (both up-regulated or both
down-regulated under N limitation).

**Tolonen 2006 (microarray):** 414 significant DE genes across timepoints
6h–48h (timepoints 0h and 3h excluded: near-zero DE genes at those timepoints).
Genes were summarized to a single direction and best directional rank across
timepoints.

**Read 2017 (RNA-seq):** 388 significant DE genes across timepoints 3h–24h.
Dataset is filtered to top 50% of genes by expression level (table_scope:
filtered_subset); genes absent from this filtered set are not testable and
are classified as "tolonen_only_read_absent" in the extended signature.

**Core signature (n=198):** Genes significant in both studies with concordant
direction — 83 up-regulated and 115 down-regulated under N limitation. These
represent the high-confidence reproducible N-limitation response.

**Extended signature (n=367):** Adds genes significant in Read only (173
read_only), Tolonen-only genes that are present but non-significant in Read
(122 tolonen_only_read_ns), and Tolonen-only genes absent from the Read
filtered subset (72 tolonen_only_read_absent). The extended signature
captures the broader response but with lower cross-study confidence.

## Statistical methods

Three metrics were computed to score Weissberg 2025 timepoints against the
N-limitation signature:

1. **Hit rate:** Fraction of signature genes that are DE in the same direction
   as defined in the signature (concordant fraction). Values near 1 indicate
   strong N-limitation signal.

2. **Mean signed log2FC:** Mean log2 fold change within the target dataset,
   sign-corrected so that positive values indicate concordance with the
   signature direction. Computed only over genes present in the target dataset.

3. **Rank score:** Cross-platform metric based on the rank of signature genes
   within the target dataset's DE ranking. Accounts for the ordering of
   differential expression rather than just significance.

**Permutation test:** 1000 permutations of random gene sets of the same size
as the core signature (n=198). Empirical p-value computed as the fraction of
permutations with a rank score equal to or greater than the observed value.
The permutation test is applied to the core signature only (n=198 > 30
minimum threshold for stable permutation statistics).

**Reference baselines:** The same three metrics were computed for the Tolonen
and Read reference datasets per timepoint, to provide an upper bound on
expected scores under true N-limitation conditions.

## Results summary

**RNA-seq axenic MED4 (day 14):** Strong N-limitation signal. Core signature
hit_rate=0.818, mean_signed_log2fc=2.68, rank_score=0.743, permutation
p=0.000. [interpretation] MED4 is strongly N-limited at day 14 in axenic
PRO99-lowN culture.

**RNA-seq coculture MED4 (days 18–89):** Near-zero scores at all timepoints,
not significant (p=0.062 to 0.672). [interpretation] No transcriptomic
N-limitation signal in coculture — the transcriptome does not reflect N
limitation when Alteromonas is present.

**Proteomics axenic MED4:** Significant at day 31 and day 89 (rank scores
0.23–0.42). [interpretation] Protein-level N-limitation signal is detectable
in axenic conditions at later timepoints.

**Proteomics coculture MED4:** Significant at all timepoints (rank scores
0.41–0.76), consistently higher than proteomics axenic. [interpretation]
Paradoxically, the proteome of coculture MED4 shows a stronger N-limitation
signature than axenic MED4.

**RNA/protein discordance in coculture:** [interpretation] The coculture
transcriptome shows no N-limitation signal while the proteome shows
persistent, strong N-limitation. This discordance is the main biological
finding: Alteromonas suppresses the transcriptional N-stress response in
MED4 (possibly by providing N compounds that relieve transcription-level
regulation) while protein turnover still reflects the underlying low-N
environment.

**Alteromonas N-recycling:** 171 candidate N-recycling genes searched in
HOT1A3; 65 were significantly DE across Weissberg 2025 timepoints. Both
up- and down-regulated, with activity increasing from day 18 through day 60,
then plateauing. [KG] Alteromonas engages N-recycling pathways progressively
over the coculture time course.

## Limitations

- **Read 2017 filtered_subset caveat:** Read 2017 is filtered to top 50% of
  genes by expression level (table_scope: filtered_subset). Genes below this
  threshold are absent from the Read DE data and cannot be confirmed in the
  cross-study intersection. These genes are classified as
  tolonen_only_read_absent in the extended signature but cannot be confirmed
  as true negatives in Read.

- **Microarray vs RNA-seq comparison:** The core signature merges Tolonen
  (microarray) and Read (RNA-seq) data. Platform differences in sensitivity
  and dynamic range may bias the intersection toward genes with large
  fold changes.

- **Proteomics coverage:** Proteomics data covers a subset of the genome
  (proteins detectable by mass spectrometry). Signature scores for
  proteomics are computed only over detected proteins, which may be a
  non-representative subset.

- **RNA/protein discordance as unexpected finding:** [gap] The strong
  N-limitation proteome signal in coculture, accompanied by no
  transcriptomic signal, was not anticipated. This may indicate
  post-transcriptional regulation, protein stability differences, or a
  genuine decoupling of transcription from protein abundance in coculture.
  Further investigation is needed to distinguish these mechanisms.

- **Single time point for axenic RNA-seq:** Only day 14 is available for
  axenic RNA-seq, limiting the ability to track temporal dynamics
  transcriptomically in the axenic condition.
