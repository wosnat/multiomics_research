# Caveats for Interpretation

Things a reader of these results needs to know before drawing conclusions.

---

## C1 — Single ontology (CyanoRak only)

All results are relative to the CyanoRak functional role ontology. CyanoRak was selected for its high genome coverage and cyanobacterial specificity (see decisions.md D1), but using a single ontology creates blind spots. Pathways not defined in CyanoRak — or defined differently — will not appear in these results. KEGG subcategory-level enrichment is planned for B2 and may surface different or complementary findings. Results should not be described as "pathway-level" analysis in general; they are CyanoRak-specific.

---

## C2 — 20% genome unannotated

Approximately 470 MED4 genes (24%) have no CyanoRak annotation and are excluded from all enrichment tests. These unannotated genes are predominantly uncharacterized hypothetical proteins. If the unannotated fraction is enriched among DE genes in any condition, that signal will be invisible in this analysis. The annotation gap is not random — hypothetical proteins are likely enriched in condition-specific or poorly conserved functions.

---

## C3 — Catch-all categories (R.2 and D.1)

R.2 "Conserved hypothetical proteins" (340 genes) and D.1 "Adaptation/acclimation" (213 genes) are the two largest CyanoRak level 1 categories. Their large size gives them high statistical power — even modest over-representation will produce small padj values. Enrichment in R.2 or D.1 should be interpreted with caution: these categories are heterogeneous and do not represent a coherent biological function. They are most useful as indicators that a condition produces a broad, diffuse DE response.

---

## C4 — No FDR correction across experiments

Benjamini-Hochberg FDR correction is applied within each experiment × timepoint × direction group. The 133 significant results are from 10 experiments tested independently, with no cross-experiment correction. With 5,589 total tests and a 5% threshold, ~280 false positives would be expected by chance if all tests were independent. The fact that significant results cluster in biologically coherent pathways (E.4, J.7, J.8, K.2) and replicate across independent experiments (Tolonen and Read both show the same N-metabolism and photosynthesis patterns) provides confidence, but this is biological replication, not statistical correction.

---

## C5 — Steglich unreliable (low power)

The Steglich high light experiment has only 198 genes in its DE universe (10% pathway coverage for most terms). Absence of enrichment in Steglich reflects low statistical power, not biological absence of signal. Steglich is retained as a visual control to confirm that the overall analysis does not produce spurious enrichments, but absence of signal in Steglich cannot be interpreted as "high light does not affect these pathways."

---

## C6 — Read N-dep background caveat

The Read N-deprivation experiment uses a `filtered_subset` DE table. The background set (genes in the DE table) is smaller than the full genome, which may reduce sensitivity for some pathways. The Read results should be interpreted as valid but potentially incomplete relative to Tolonen, which uses the full genome background.

---

## C7 — Signed score loses information when both directions are significant

The signed enrichment score assigns a single direction (up or down) to each pathway × experiment × timepoint based on which direction has the smaller padj. When both up and down are independently significant, one direction's signal is discarded. This edge case was uncommon in these data, but could occur in large heterogeneous pathways (R.2, D.1) where up-regulated and down-regulated subsets both exceed the threshold. The full test results (both directions) are available in `results/enrichment_all.csv` for any case requiring both directions.
