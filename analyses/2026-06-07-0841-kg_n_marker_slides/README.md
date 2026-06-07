# KG N-stress marker slides — cynS featured, urtA contrast

Conference figures built from the multiomics KG. Question driving the slides:
**"How are the N-acquisition genes (amt / urt / cynS) differentially expressed
across conditions, and can we use them to detect N stress?"**

Organisms: Prochlorococcus (MED4 as the data-rich anchor) + Synechococcus
(via cross-genus ortholog groups). All values are `[KG]`-sourced; interpretation
is tagged `[interpretation]`.

## The figures

| File | What it shows |
|---|---|
| `figures/fig6_urtA_family_flower.*` | **★ Hero / impact slide (radial "flower").** urtA family at the centre, genomes as petals, experiments at the rim; **nitrogen experiments highlighted in gold**. Graphviz `twopi` (mermaid can't do radial). Best for a title/impact slide. |
| `figures/fig7_cynS_family_flower.*` | **cynS flower (contrast).** Same radial style; hub reads "DE data in only 3 of 10 genomes." Makes the "cleanest marker, least measured" point next to fig6. |
| `figures/fig8_n_specificity_forest.*` | **Step 2 — statistical specificity test.** Forest plot of Fisher odds ratios (DE/not-DE, no fold-change) asking "is up-regulation enriched under N?" Pooled urtA family is significant (OR≈9, q=0.011); see `notebook_step2_specificity.md`. |
| `figures/fig5_urtA_family_integration.*` | **★ "Power of the KG" detail slide.** One gene family → 7 genomes → each genome's fan of experiments, edges colored up/down/mixed/not-significant, real author-year citations. Data integrated across **38 experiments and 18 publications** (rim shows a representative subset). No locus tags / IDs. Best for the walk-through. |
| `figures/fig1_cynS_med4.*` | **Featured single gene.** cynS (PMM0373) DE edges to every MED4 condition tested. One green arrow to N stress; red/gray everywhere else → the *specific* detector. |
| `figures/fig2_urtA_med4_contrast.*` | **Contrast.** urtA (PMM0970) responds in all directions to many stresses → *sensitive but promiscuous*. |
| `figures/fig3_cynS_orthologs.*` | cynS ortholog group (CK_00001552) — DE edges in the 2 strains that have data; 8 orthologs have no DE experiments in the KG. |
| `figures/fig4_urtA_orthologs.*` | urtA ortholog group (CK_00000076) — N-induction conserved in MED4 + MIT9313 (HL + LL ecotypes); 13 orthologs lack DE data. |

Each figure is a `.mmd` (mermaid source), `.svg` (vector, for slides), and `.png`
(3× raster). Edit the `.mmd` and re-render with:
```bash
node_modules/.bin/mmdc -i figures/fig1_cynS_med4.mmd -o figures/fig1_cynS_med4.svg -b transparent
```

## Edge colour legend (unified across all figures)

- **Red** = significant **up** · **Blue** = significant **down**
- **Light purple, dashed** = **mixed** (up and down across experiments/timepoints)
- **Gray, thin** = tested but **not significant** (a real "no response" from
  `all_detected_genes` tables — an honest absence, not missing data)
- **Gold node outline** (flowers, fig6/fig7) = **nitrogen-stress experiment**
- **Dotted faded** (fig3/fig4) = orthologs / strains with no DE experiment in the KG

## Slide narrative (suggested)

1. **The KG as a graph** — a gene is a node; every differential-expression result
   is an edge to a condition. Colour the edges by direction and the biology reads
   off the picture.
2. **fig1 — cynS is a clean N-stress reporter.** Up to **+4.7 log2FC**, significant
   in 4 of 8 N experiments, and essentially flat (gray) under iron, salt, viral,
   and dark stress. The only non-N signal is a modest down-shift under low-CO₂.
3. **fig2 — not every N gene is a good marker.** urtA has the *biggest* N response
   (+5.4) but also moves under carbon, phosphate, light/dark, viral, and coculture,
   and even flips direction under N depending on the experiment. High sensitivity,
   low specificity.
4. **fig3 / fig4 — orthologs and the data gap.** The cynS N-response is only
   measured in MED4; the urtA N-response is conserved across two ecotypes
   (MED4 + MIT9313). [interpretation] The cleanest marker is also the least
   measured across strains — the KG makes that gap explicit, and it's a concrete
   "what to measure next" ask.
5. **Take-home.** `[interpretation]` A small N-acquisition panel can report N stress,
   but marker choice is a specificity/sensitivity trade-off: **cynS/cynD** (specific),
   **amt1** (intermediate), **urtA** (sensitive, promiscuous). Reading the panel as
   a set — not any single gene — is the robust detector.

## Key numbers (all `[KG]`, Prochlorococcus MED4)

cynS (PMM0373), tested in 8/8 N experiments:
- N stress: **UP in 4** (best rank #1; max **+4.67** log2FC, proteomics coculture),
  not significant in the other 4 (the Tolonen microarrays + one RNA-seq), **0 down**.
- Low-CO₂ carbon: down (−1.0, −1.1). Coculture: +1.7 (Weissberg) and −13.4* (ismej.2016.70).
- Iron / salt / viral / dark: tested, not significant.

urtA (PMM0970):
- N stress: 14 timepoints up (max **+5.37**), but **2 timepoints down** → inconsistent.
- Also: carbon (+2.1 / −3.1), phosphate (−1.9), dark (−1.9), viral (+1.4), coculture (+2.5/mixed).

cynD (PMM0372, cyanate transporter, same operon as cynS): N up in 5 experiments
(max +3.27), including the Tolonen N-deprivation microarray where cynS itself was
n.s. → the *cyn operon as a unit* is a more robust N reporter than cynS alone.

\* The cynS −13.4 log2FC in the standalone coculture RNA-seq (ismej.2016.70) is a
single extreme value, almost certainly a detection-floor / near-absence artifact,
not a 10,000-fold biological repression. Flagged, not used for the marker claim.

## Caveats `[KG]` / `[interpretation]`

- "Significant" uses each publication's own threshold, not a uniform padj cutoff
  (cross-study magnitude comparisons are approximate).
- Platforms are mixed (RNA-seq, microarray, proteomics); microarrays are less
  sensitive, which partly explains cynS n.s. in the Tolonen N experiments.
- Ortholog DE coverage is uneven: cynS has DE data in 2/10 orthologs, urtA in 6/19.
  "No DE in KG" means no experiment, not no response.
- amt1 (PMM0263) was characterised in discovery but not drawn here — available as a
  third panel if wanted (intermediate specificity: N up, but also light & viral up).

## Data completeness / QA (checked 2026-06-07)

- **Family totals** (from `differential_expression_by_ortholog`, `significant_only=false`):
  - urtA group `CK_00000076`: **7 genomes · 38 experiments · 18 publications** (89 rows: 27 up / 21 down / 41 n.s.).
  - cynS group `CK_00001552`: **3 genomes · 17 experiments · 8 publications** (44 rows: 11 up / 6 down / 27 n.s.).
- **Pitfall fixed:** an earlier pass used `significant_only=true`, which dropped the n.s. rows and hid **Synechococcus WH7803** (urtA + cynS both tested, both not-significant, Christie-Oleza 2017 co-culture proteomics). WH7803 is now shown as a gray "tested, no change" genome in fig5/fig6/fig7.
- **Truncation:** the full pull returns 133 rows, `returned=133, truncated=false` — complete. (The large-result file was read via `jq`, not partial reads.)
- **Paralog / sibling groups not used as the family anchor (deliberate):** cyanate hydratase also has a Synechococcus-only group `CK_00003051` (cynH); ammonium has `CK_00008701` (amt2) beside amt1 `CK_00000244`. urtB–E are their own groups (`CK_00001365/1366/1367/8074`). The flowers anchor on the single curated cross-genus group per gene (cynS, urtA, amt1), which is the right "gene family" unit; eggnog/COG-level groups are broader (COG0004 pulls in Alteromonas etc.) and were intentionally excluded.
- The flower rims show a **representative** subset of experiments; the hub counts are the full-family totals.

## Provenance

- Underlying edges + sources: `data/de_edges.csv`
- KG release: dev build `0.0.0-dev`, built 2026-06-03 (120,416 genes / 197 experiments / 43 papers)
- Tools used: `genes_by_function`, `list_experiments`, `gene_response_profile`,
  `gene_homologs`, `genes_by_homolog_group`, `differential_expression_by_ortholog`
- Primary N publications in play: Weissberg 2025 (10.1101/2025.11.24.690089),
  Tolonen 2006 (10.1038/msb4100087), Read 2017 (10.1038/ismej.2017.88).
