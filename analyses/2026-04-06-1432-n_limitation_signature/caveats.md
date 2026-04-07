# Caveats for Interpretation

Things a reader of these results needs to know before drawing conclusions.

## Platform coverage differs

- RNA-seq covers 1,849 MED4 genes; proteomics covers 1,424 genes.
- Of the 198 core signature genes, not all are detected in both platforms.
- Scores across platforms are NOT directly comparable in magnitude — only direction and rank are comparable.
- The RNA/protein discordance finding must be evaluated in the context of gene coverage: are the genes driving the proteomics signal even detectable in the RNA-seq data?

## Axenic RNA-seq is a single timepoint

- Axenic MED4 RNA-seq has one contrast (starvation vs exponential, ~day 11), not a time course.
- Cannot assess whether the axenic transcriptomic N-limitation signal changes over time.
- Proteomics axenic has 3 timepoints (day 14, 31, 89), providing some temporal resolution.
- The RNA-seq coculture trajectory (d18-89) has no axenic trajectory to compare against — only a single reference point.

## Signature is defined from acute N-stress, applied to chronic

- Tolonen (0-48h) and Read (3-24h) measure acute N-deprivation responses.
- Weissberg measures chronic N-limitation over weeks to months.
- Some acute-response genes may have returned to baseline by day 18 (the first Weissberg coculture timepoint). This would show up as low hit rate — which is informative, not a failure.
- Genes that respond only to chronic N-stress (not acute) would be missed entirely by this signature. These are candidates for Approach C's Tier 2.

## Permutation test uses the gene's own direction labels

- The permutation test shuffles gene labels but keeps the direction assignments from the signature.
- This means the null distribution is generated under the assumption that roughly equal numbers of genes are labeled "up" and "down" — which is true (83 up, 115 down).
- For the reference baselines (Tolonen, Read), the test is somewhat circular because the signature directions were DERIVED from these same experiments. The test is valid for Weissberg experiments (independent data).

## The "days 60+89" combined timepoint

- Weissberg coculture experiments include a combined "days 60+89" timepoint alongside individual day 60 and day 89.
- This timepoint is included in data tables but excluded from trajectory plots to avoid visual confusion.
- It is included in score tables — scores should be similar to the average of day 60 and day 89.

## Read 2017 filtered subset

- Read's DE data is filtered to the top 50% of genes by expression level (table_scope: filtered_subset).
- Genes absent from this dataset may be genuinely unexpressed or may be expressed but below the 50th percentile.
- This affects extended signature classification: "tolonen_only_read_absent" (72 genes) may include genes that WOULD have been significant in Read if the full dataset were available.
- The core signature is not affected (core genes are present and significant in both studies).

## Alteromonas N-recycling is exploratory

- The 171 candidate N-recycling genes were found by broad text search, not curated pathway annotation.
- "Significant" means DE relative to exponential growth — not necessarily related to N-recycling specifically.
- The bidirectional pattern (both up and down genes) may reflect general metabolic remodeling, not specifically N-recycling.
- The temporal pattern (increasing activity day 18 → 60) is suggestive but needs pathway-level analysis (Approach B) to confirm which functions are driving it.
