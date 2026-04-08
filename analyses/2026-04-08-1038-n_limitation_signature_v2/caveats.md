# Caveats for Interpretation

Things a reader of these results needs to know before drawing conclusions.

## Platform coverage differs and is biased

- RNA-seq covers 1,849 MED4 genes; proteomics covers 1,424 genes.
- Of 189 core signature genes, 185 are in RNA-seq datasets and 147 in proteomics (42 missing).
- The 42 missing proteomics genes are NOT random — they are biased toward small/membrane proteins (hli, photosystem subunits, ribosomal proteins) with better-than-average signature ranks (median 16 vs 29 for detected genes).
- This means proteomics scores are computed over a non-representative subset that underrepresents small, strongly-responsive stress proteins.
- Scores across platforms are NOT directly comparable in magnitude. Use direction, significance, and hit rate for cross-platform comparison.

## RNA/protein discordance is genuine, not a coverage artifact

- Tested by scoring RNA-seq with only the 147 proteomics-detectable genes: RNA-seq axenic still scores 0.548 vs proteomics 0.066 for the same gene set.
- The discordance reflects real differences in DE at the transcript vs protein level, not differential gene coverage.

## Signature derived from acute stress, applied to chronic

- Tolonen (6-48h) and Read (12-24h) measure acute N-deprivation responses.
- Weissberg measures chronic N-limitation over weeks to months.
- Genes that respond only to chronic N-stress (not acute) would be missed by this signature.
- Low scores in coculture RNA-seq could reflect either (a) genuine absence of N-stress transcription or (b) return to baseline of acute-response genes under chronic conditions.

## Axenic RNA-seq is a single timepoint

- Axenic MED4 RNA-seq has one contrast (~day 14 starvation vs exponential), not a time course.
- Cannot assess whether the axenic transcriptomic N-limitation signal changes over time.
- Proteomics axenic has 3 timepoints (day 14, 31, 89) providing some temporal resolution.

## Reference study self-scoring is circular

- Tolonen and Read per-timepoint scores define the upper bound (positive control range) but are circular — the signature was derived from these experiments.
- P-values for reference experiments will be very low by construction — this is expected, not evidence of independent validation.
- Use Weissberg scores (independent data) for biological conclusions.

## Read 2017 filtered_subset

- Read's DE data filtered to top 50% of genes by expression level (table_scope: filtered_subset).
- 73 Tolonen-only genes absent from Read may be expressed but below the expression threshold, not truly absent.
- Core signature (189 genes) includes only genes present and significant in both studies — this caveat affects the extended classification only.

## Steglich high white light control is unreliable

- Only 43 of 189 core signature genes detected (198-gene filtered_subset).
- Score of 0.400 is based on too few genes and 146 absent genes.
- Cannot meaningfully test N-specificity with this dataset. Included for completeness but should not be used for conclusions.

## Early timepoints excluded from both references

- Tolonen 0h/3h and Read 3h excluded for showing transient responses that reverse by later timepoints.
- This means the signature captures the sustained (6h+) N-deprivation response, not the immediate (0-3h) response.
- Genes with exclusively early response patterns are absent from the signature.

## Permutation test limitations

- 1000 permutations — smallest possible non-zero p is 0.001.
- Two-tailed test — detects both activation and reversal as significant.
- Requires ≥ 30 matched genes — top tier (38 genes) may have fewer after absent genes, returning NaN p-values for some experiments.
- Random gene sets use the same up/down direction ratio as the real signature — tests gene identity, not direction assignment.

## "Days 60+89" combined timepoint

- Weissberg coculture includes a combined "days 60+89" timepoint alongside individual day 60 and day 89.
- Included in score tables, excluded from trajectory plots to avoid visual confusion.
