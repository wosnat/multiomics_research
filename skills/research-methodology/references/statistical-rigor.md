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

### Enrichment: background set construction

The background set is the universe of genes that *could* have been
detected. Getting this wrong inflates or deflates enrichment
p-values — it's the most common enrichment error.

**Which background to use depends on the question:**

| Question | Background set | How to get it |
|----------|---------------|---------------|
| "Are nitrogen-responsive genes enriched for GO term X?" | All genes tested in the experiment | `differential_expression_by_gene(organism=..., experiment_ids=[...], limit=None)` — use all returned locus tags, not just significant ones |
| "Are genes in cluster Y enriched for function X?" | All genes in the organism's genome | `genes_by_function("*", organism=..., limit=None)` — returns all annotated genes |
| "Is ontology term X over-represented in my gene set?" | All genes with any annotation in that ontology | `genes_by_ontology(term_id=<root_term>, organism=..., limit=None)` |

**Common mistakes:**
- Using only significant genes as background (guarantees
  inflated enrichment)
- Using the whole genome when only a subset was measured
  (e.g., proteomics covers fewer genes than RNA-seq)
- Not filtering background to the same organism as the test set
- Forgetting that `differential_expression_by_gene` returns only
  genes with DE edges — genes with no expression data at all are
  not in the background

**Worked example (Fisher's exact):**

Suppose 800 MED4 genes were tested in experiment E1. Of these,
120 are significantly upregulated. Of the 800 total, 40 are
annotated to GO term "response to oxidative stress". Of the 120
significant genes, 12 have this annotation.

```
                  In GO term    Not in GO term
Significant          12              108         = 120
Not significant      28              652         = 680
                     40              760         = 800

Fisher's exact test on this 2×2 table.
Expected: 120 × 40/800 = 6.0
Observed: 12 → enrichment ratio = 2.0×
```

Report: "response to oxidative stress" enriched 2.0-fold among
upregulated genes (Fisher's exact p = X, BH-corrected q = Y,
12/40 annotated genes significant vs 6.0 expected).

Always report: fold enrichment, raw p-value, corrected q-value,
counts (k/K from n/N), and the background set definition.

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
