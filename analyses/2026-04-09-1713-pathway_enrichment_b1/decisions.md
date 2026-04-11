# Decision Log

Design decisions with rationale — WHY the analysis was done this way, not what it does (that's in methods.md).

---

## D1 — Ontology selection: CyanoRak level 1

**Decision:** Use CyanoRak level 1 (110 terms, 80% genome coverage, median 9 genes/term) as the enrichment ontology for B1.

**Alternatives considered:**
- CyanoRak level 2: only 18% genome coverage — most genes are annotated at level 1, not level 2. Level 2 terms are niche subdivisions (vitamin subclasses, specific pigments) that don't exist for most functional categories.
- KEGG subcategory (level 1): 41% genome coverage (constant across all levels — a structural gap where KEGG-annotated genes lack pathway linkage). Better coverage than GO but much lower than CyanoRak. Retained as a candidate for B2.
- GO BP: DAG structure creates ancestor terms that absorb most genes (637g for "primary metabolic process" at level 3). No level combines good coverage + manageable maximum term size.
- EC numbers: enzyme-only, 43% genome coverage. Misses transporters, structural proteins, regulators.
- COG category, GO CC, Pfam, TIGR role: not competitive on genome coverage or term granularity.

**Critical method note:** The genome_coverage metric (genes reachable at this level / total genome) is essential for hierarchy level selection. Before adding this metric, the initial analysis selected CyanoRak level 2 based on appealing term sizes (median 12). The genome_coverage calculation revealed level 2 covers only 18% of the genome — the correct level is 1. This bug was caught during the characterization step (step 2) before enrichment was run.

---

## D2 — Background set: per experiment × timepoint

**Decision:** Use the set of genes present in each experiment × timepoint DE table as the background for Fisher's exact tests, rather than the full 1,976-gene genome.

**Rationale:** Enrichment tests should assess whether DE genes are enriched relative to the genes that *could* have been DE — i.e., genes that were quantified and passed QC in that experiment. Using the full genome as background would include genes that were never measured, inflating the denominator and underestimating enrichment. The `table_scope` (genes in the DE table) varies by experiment and is the correct biological background.

---

## D3 — Minimum pathway size: 5 genes

**Decision:** Include only pathways with ≥5 genes in the experiment background.

**Rationale:** Fisher's exact test is unreliable for very small gene sets — a single DE gene in a 2-gene pathway produces OR = ∞. 5 genes is a common minimum for enrichment analyses. This reduced the 110-term CyanoRak level 1 to 69 testable pathways.

---

## D4 — Signed enrichment score

**Decision:** Summarize each pathway × experiment × timepoint with a signed score: sign × -log10(padj), where sign is determined by the dominant enrichment direction (up or down). Use this as the primary visualization metric.

**Rationale:** Separate up and down heatmaps require side-by-side comparison to see discordance and direction. The signed score encodes magnitude (size of the cell) and direction (sign) in a single value, enabling a single heatmap to show the full enrichment landscape. When both up and down are significant, the direction with the smaller padj wins — this occurs rarely in these data (see caveats.md for the information loss caveat).

---

## D5 — Cyanate reclassification as partial negative control

**Decision:** Treat Tolonen cyanate as a "partial negative control" rather than a clean negative control.

**Rationale:** The original design treated cyanate (a poor N source that MED4 can use but does not prefer) as a negative control for N-limitation. The enrichment results show that cyanate activates photosynthesis downregulation (J.1, J.2, J.7 enriched down, padj < 0.05) without activating N-metabolism upregulation (E.4 not enriched). This is biologically coherent — cyanate causes mild N-stress at the photosynthesis level but does not produce the full N-limitation response. Cyanate is therefore useful as a partial control: it confirms that photosynthesis downregulation is not unique to true N-deprivation, while confirming that N-metabolism enrichment distinguishes true N-deprivation from alternative N sources.
