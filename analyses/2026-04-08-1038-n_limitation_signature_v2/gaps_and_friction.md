# Gaps and friction points

## KG gaps

1. **No well-powered non-N stress experiments for MED4.** Phosphate (34 genes), iron (112 genes), salt (28 genes), viral (25-45 genes) — all too sparse for signature-based scoring. Carbon experiments have no non-significant rows. The only usable non-N control (Steglich high white light) has 198 genes. Limits our ability to test N-specificity of the signature.

2. **Read 2017 filtered_subset scope.** Top 50% of genes by expression — 857 vs 1,697 genes in Tolonen. 73 signature-relevant genes absent from Read may be expressed but below the threshold. No way to recover them without the full dataset.

## MCP/API friction

1. **`list_experiments` doesn't return matching_genes.** Had to call `differential_expression_by_gene(summary=True)` per experiment to get gene counts. A `matching_genes` field in `list_experiments` results would save N+1 queries during scoping.

## Methodology insights

1. **Rank normalization matters.** Normalizing by total_genes (v1 approach) compresses all ranks to ~0.95+, making rank score ≈ hit rate. Normalizing by n_significant_in_direction gives meaningful spread. This was discovered during interactive review — would have been missed without the walkthrough.

2. **Early timepoint exclusion is symmetric.** Both Tolonen 0h/3h and Read 3h show transient responses. v1 excluded only Tolonen's early timepoints. Excluding Read 3h fixed the rbcL misclassification.

3. **Proteomics gene coverage bias is systematic.** Missing 42 genes aren't random — they're small/membrane proteins (hli, photosystem subunits, ribosomal) with better-than-average ranks. This biases proteomics scoring toward larger, more detectable proteins.

## Process retrospective

### What worked
1. **Methodology skill loaded before brainstorming.** Artifacts structure, notebook discipline, and source tagging shaped the design from the start.
2. **Interactive exploration at every step.** Caught the rbcL discordance, Read 3h transient, rank normalization issue, proteomics coverage bias — all through walking through specific genes.
3. **Toy-data verification before real data.** 31 tests caught the normalization change cleanly.
4. **Marker gene tracing.** Six genes across every step made the method transparent and caught edge cases.
5. **Rank normalization fix.** The walkthrough surfaced a fundamental metric issue that would have been invisible in v1's three-metric approach.

### What didn't work
1. **Steglich light control too sparse.** Should have checked gene counts before including. The 198-gene filtered_subset is unusable for 189-gene signature scoring.
