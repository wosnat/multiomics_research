# Statistical rigor

## What the KG provides

- log2FC and padj from the original DESeq2 analysis
- Per-study statistics — do not combine p-values across studies
- Precomputed summary counts (gene counts, significant counts)

## What you must compute in scripts

- Enrichment analysis (Fisher's exact, hypergeometric) with proper
  background set from the KG
- Multiple testing correction when comparing across conditions
- Effect size reporting (not just significance)
- Volcano plots, heatmaps with proper clustering

## Worked examples in specs

Every formula in a spec must include at least one worked example
with concrete numbers showing input → computation → output.

Shorthand like "reference direction as expected sign" is not
sufficient — a researcher must be able to verify the formula
produces the expected result by reading the example. If there are
edge cases (e.g., what happens when a gene is not significant?),
include an example for those too.

## What to flag

- Studies without adjusted p-values (report as fold-change only)
- Borderline significance (padj near 0.05)
- Duplicate KG entries for the same gene/condition (different
  contrasts — investigate, don't cherry-pick)
- Small sample sizes or missing replicates
- Conflation of paralogs or orthologs in aggregated results

## Strength of language vs strength of evidence

- padj < 0.001: "strongly significant" / "highly significant"
- padj < 0.01: "significant"
- padj < 0.05: "nominally significant" or "borderline significant"
- padj = 0.05: "at the conventional threshold" — flag explicitly
- No padj: "fold-change of X (no statistical test available)"
- Never write "consistent with [mechanism]" for data without
  statistical support — use "suggestive of" or "potentially
  related to"

## Cross-experiment comparison caveats

- Compare direction (up/down) and rank across platforms, not
  magnitude
- Compare log2FC magnitudes only within the same platform and study
- When presenting cross-experiment tables, include a platform
  column and note the caveat
- Use `gene_response_profile` rank fields (comparable across
  platforms) instead of raw fold changes for cross-study comparisons
- Expression data from different platforms (microarray, RNA-seq)
  or different statistical tests (Goldenspike, Rockhopper, DESeq2)
  are NOT directly comparable

## Cross-study p-value comparison

P-values from different studies have different designs, sample
sizes, and statistical power. They are not comparable:
- Never compare p-values across studies
- Compare effect sizes (log2FC) across studies if needed, but
  note the caveat
- For cross-study claims, use direction and rough magnitude,
  not precise p-value ranking

## Summary statistics

Summary statistics must be computed in scripts over complete data,
not eyeballed from MCP output:
- If you report a summary statistic, cite the script that computed it
- For quick estimates in chat, explicitly say "rough estimate from
  the top N results" and never present as a finding
