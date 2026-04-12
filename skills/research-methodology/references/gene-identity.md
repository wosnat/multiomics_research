# Gene identity rules

## Locus tags, not gene names

Gene names are ambiguous. Always:
- Resolve gene names to locus tags early (`resolve_gene`,
  `genes_by_function`)
- Report locus tags in all tables and analysis outputs
- Use gene names as labels alongside locus tags, never as the sole
  identifier

### Gene discovery tools

- By name/keyword: `resolve_gene`, `genes_by_function`
- By pathway/ontology: `search_ontology` → `genes_by_ontology`
- By homology: `search_homolog_groups` → `genes_by_homolog_group`
- By prior knowledge: start with `resolve_gene` for known names

## Paralog handling

When paralogs exist (multiple locus tags share a gene name):
- Treat each locus tag as a separate entity throughout
- Never aggregate expression data across locus tags that share
  a name without explicitly noting they are paralogs
- In tables: always include locus_tag column; add a suffix like
  "katA-1", "katA-2" when presenting results
- Use `gene_overview` on all locus tags to confirm they are
  distinct genes

See [Anti-hallucination — Paralog conflation](anti-hallucination.md#11-paralog-conflation)
for a concrete failure example.

## Ortholog cluster rules

Genes in the same ortholog cluster are NOT interchangeable:
- Ortholog clusters are for cross-organism comparison of
  *conserved function*, not for within-organism gene equivalence
- When reporting per-gene expression, use locus tags, not cluster
  membership
- If a gene has no expression data, say so — don't substitute
  a cluster-mate's data

See [Anti-hallucination — Ortholog cluster conflation](anti-hallucination.md#12-ortholog-cluster-conflation)
for a concrete failure example.

## Gene function claims

Gene function claims must come from `gene_overview`,
`gene_ontology_terms`, or `genes_by_function` output — not from
training data. If intrinsic knowledge disagrees with the KG
annotation, report the KG annotation and note the discrepancy.

Prefix intrinsic knowledge with: "Based on general knowledge
(not from this KG)..."
