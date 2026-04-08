# N-Limitation Signature Analysis v2: Methods

## Research question

Can we quantify the degree of nitrogen limitation in Prochlorococcus MED4 molecularly, using a gene signature derived from independent N-stress reference studies in the KG, and show that axenic vs coculture conditions differ?

## Data scope

### Reference studies (signature building)
- Tolonen et al. 2006 (DOI: 10.1038/msb4100087) — MED4 N-deprivation time course (microarray, 6h-48h). Continuous light, 50 µmol m⁻² s⁻¹, 22°C. Timepoints 0h and 3h excluded (near-zero DE). Statistical test: Goldenspike. Table scope: all_detected_genes (1,697 genes).
- Read et al. 2017 (DOI: 10.1038/ismej.2017.88) — MED4 N-depleted vs N-replete Pro99 (RNA-seq, 12h-24h). Continuous light, ~55 µmol m⁻² s⁻¹, 25°C. Timepoint 3h excluded (transient up-bias reverses by 12h). Statistical test: Rockhopper. Table scope: filtered_subset (top 50% by expression, 857 genes).

### Target study
- Weissberg et al. 2025 (DOI: 10.1101/2025.11.24.690089) — MED4 in PRO99-lowN, axenic vs coculture with Alteromonas HOT1A3. Continuous light, 20 µmol m⁻² s⁻¹, 24°C. RNA-seq (1,849 genes) + proteomics (1,424 genes). Axenic: day 14 (RNA-seq single tp), days 14/31/89 (proteomics). Coculture: days 18/31/60/89 (both platforms).

### Negative controls
- Tolonen 2006 cyanate as sole N source (1,697 genes, all_detected_genes) — N-sufficient
- Tolonen 2006 urea as sole N source (1,697 genes, all_detected_genes) — N-sufficient
- Aharonovich & Sher 2016 (DOI: 10.1038/ismej.2016.70) — MED4 coculture with HOT1A3 under normal growth (1,714 genes, all_detected_genes) — tests coculture effect without N-limitation
- Steglich et al. 2006 (DOI: 10.1128/JB.01097-06) — MED4 high white light 55 µmol (198 genes, filtered_subset) — non-N stress specificity test

### Inclusions/exclusions
- Tolonen 0h/3h and Read 3h excluded: transient early responses, not sustained N-deprivation
- Weissberg "days 60+89" combined timepoint included in data but excluded from trajectory plots
- All studies used continuous light — diel confound investigated and ruled out

## Gene selection

### Core signature (189 genes)
Genes significant in BOTH Tolonen (6h-48h) and Read (12h-24h) with concordant direction (both up or both down under N limitation). 74 up-regulated, 115 down-regulated.

**Per-gene summarization:** For each gene in a multi-timepoint study, majority direction across significant timepoints was assigned. Ties broken by best directional rank (lower rank = stronger signal wins). Best directional rank = minimum rank_up (for up genes) or rank_down (for down genes) across timepoints.

**Cross-study intersection:** Genes present and significant in both studies with same direction → core. Genes with opposite directions → discordant (4 genes, excluded). Genes significant in only one study → extended.

**Extended signature classification (symmetric):** One-study-only genes classified as `_ns` (present but not significant in the other study) or `_absent` (not in the other study's dataset). Applied symmetrically to both study-A-only and study-B-only genes.

### Extended signature (342 genes)
- tolonen_only_read_ns: 143 genes (present in Read but not significant)
- read_only_tolonen_ns: 114 genes (present in Tolonen but not significant)
- tolonen_only_read_absent: 73 genes (absent from Read's filtered subset)
- read_only_tolonen_absent: 12 genes (absent from Tolonen's microarray)

### Discordant (4 genes)
Significant in both studies but opposite direction. Excluded from signature.

## Three scoring tiers
- **Top (38 genes):** Core genes with cross_study_best_dir_rank ≤ 10. Strongest markers.
- **Core (189 genes):** Full concordant intersection. Primary metric.
- **Extended (531 genes):** Core + all extended genes.

## Statistical methods

### Rank score
For each signature gene in the target dataset:
1. **Concordance:** +1 if DE direction matches signature, -1 if opposite, 0 if not significant
2. **Normalized rank:** 1 - (dir_rank / n_significant_in_direction). Normalizing by n_significant (not total genes) ensures rank 10 of 300 significant_up genes (0.967) is meaningfully different from rank 200 of 300 (0.333).
3. **Contribution:** concordance × normalized_rank (0 for not-significant genes)
4. **Score:** mean(contribution) over all matched genes (absent genes excluded)

Score range: [-1, +1]. Positive = N-limitation signal. Zero = no signal.

### Hit rate
n_concordant / n_matched. Fraction of detectable signature genes responding in the expected direction. Captures breadth; rank score captures strength.

### Permutation test
1000 random gene sets of equal size, same direction ratio (up/down). Empirical p-value = fraction of null |scores| ≥ observed |score|. Two-tailed. Minimum 30 matched genes required. For reference studies (self-scoring), low p-values are expected circularity.

### Breakdown reported per condition
n_matched, n_absent, n_concordant, n_reversed, n_not_significant, hit_rate — essential context for interpreting scores across platforms with different gene coverage.

## Results summary

### RNA-seq axenic (day 14)
Strong N-limitation signal. Core rank_score=0.583, hit_rate=0.85, p=0.000. 158 concordant, 3 reversed, 24 not significant. MED4 transcriptome strongly mirrors the reference N-deprivation response.

### RNA-seq coculture (days 18-89)
No transcriptomic N-limitation signal. Core scores 0.01-0.07, p=0.01-0.42. Most signature genes not significant. Alteromonas eliminates the MED4 transcriptional N-stress response.

### Proteomics axenic
Only day 14 significant (score=0.066, p=0.000, hit_rate=0.15). Days 31/89: equal concordant and reversed (~35 each) — no directional signal.

### Proteomics coculture
Significant at all timepoints (p=0.000). Day 31 peak: score=0.217, hit_rate=0.36, 53 concordant, 1 reversed. Stronger and more asymmetric signal than proteomics axenic.

### RNA/protein discordance
The coculture transcriptome shows no N-limitation signal while the proteome shows persistent activation. Tested with the same 147 genes detectable on both platforms: RNA-seq axenic scores 0.548 vs proteomics 0.066 for the same gene set. This is genuine DE difference, not a coverage artifact.

### Negative controls
- Tolonen cyanate: score=0.054 (weak — shares N-assimilation biology)
- Tolonen urea: score=0.019 (clean negative)
- Aharonovich coculture: score=-0.013 (clean negative)
- Steglich high white light: score=0.400 (but only 43/189 genes matched — unreliable)

### Tier comparison
Top > Core > Extended consistently across conditions. The strongest 38 markers amplify the signal; the 342 extended genes dilute it.

## Limitations

- **Read 2017 filtered_subset:** 857 of ~1900 MED4 genes. 73 Tolonen-only genes absent from Read may be below the 50th percentile expression cutoff, not truly absent.
- **Proteomics coverage bias:** 42 of 189 core signature genes undetected by proteomics. The missing set is biased toward small/membrane proteins (hli, photosystem subunits, ribosomal) with better-than-average signature ranks (median 16 vs 29). Proteomics scores are computed over a non-representative subset.
- **Signature from acute stress, applied to chronic:** Tolonen (6-48h) and Read (12-24h) measure acute N-deprivation. Weissberg measures chronic N-limitation over weeks-months. Genes specific to chronic stress would be missed.
- **Single axenic RNA-seq timepoint:** Cannot assess transcriptomic temporal dynamics in axenic condition.
- **Steglich light control unreliable:** Only 43 of 189 signature genes in the 198-gene filtered subset. Score is based on too few genes.
