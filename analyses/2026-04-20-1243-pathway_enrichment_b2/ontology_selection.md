# Ontology + level selection for pathway enrichment B2

**Date:** 2026-04-20
**Selected by:** interactive review with researcher during Step 1b MCP orientation
**Locked at:** Task 3 Commit 1 (do-phase)

All numbers come from `ontology_landscape(organism=..., experiment_ids=[classified set])` — see `data/landscape_<org>.csv`. N-term counts come from `search_ontology("nitrogen", ontology=X)` — see `data/nitrogen_ontology_search.csv`. MED4-gene counts per key-pathway term come from `genes_by_ontology(term_ids=[...], organism="MED4")` — logged in `exploration/key_pathways.csv`.

## Ranked landscape (MED4, 13 selected experiments: T+R+PC+NC)

Top 10 `(ontology, level)` pairs by `relevance_rank`:

| rank | ontology | level | tree | n_terms | genome_cov | median_genes/term |
|---|---|---|---|---|---|---|
| 1 | tigr_role | 0 | — | 77 | 87% | 13 |
| 2 | **cyanorak_role** | **1** | — | **69** | **73%** | **14** |
| 3 | go_mf | 2 | — | 35 | 58% | 25 |
| 4 | go_mf | 3 | — | 72 | 56% | 19.5 |
| 5 | go_bp | 3 | — | 68 | 54% | 13 |
| 6 | go_bp | 4 | — | 132 | 50% | 13 |
| 7 | go_mf | 4 | — | 85 | 48% | 16 |
| 8 | cyanorak_role | 0 | — | 17 | 76% | 77 (too broad) |
| 9 | brite | 0 | transporters | 3 | 4% | 30 |
| … | | | | | | |
| 32 | **kegg** | **2** | — | **97** | **38%** | **12** |

## Ranked landscape (other organisms with ≥1 selected experiment)

**Prochlorococcus MIT9313** (CTX2, 1 experiment):

| rank | ontology | level | n_terms | genome_cov | median_genes/term |
|---|---|---|---|---|---|
| 1 | tigr_role | 0 | 82 | 71% | 16 |
| 2 | cyanorak_role | 1 | 79 | 55% | 16 |
| 5 | go_bp | 3 | 76 | 43% | 14.5 |
| 31 | kegg | 1 | 30 | 28% | ~ |

**Prochlorococcus SS120** (CTX3, 1 experiment):

| rank | ontology | level | n_terms | genome_cov | median_genes/term |
|---|---|---|---|---|---|
| 1 | tigr_role | 0 | 73 | 83% | 13 |
| 9 | cyanorak_role | 1 | 72 | 68% | 14 |
| 7 | go_bp | 3 | 71 | 55% | 12 |
| 21 | kegg | 0 | 5 | 15% | ~ |

## Nitrogen term density per ontology (from `search_ontology("nitrogen")`)

| ontology | N-term matches | best-scoring term (id, name, level) |
|---|---|---|
| cyanorak_role | 2 | `cyanorak.role:E.4` Nitrogen metabolism (L1) |
| tigr_role | 1 | `tigr.role:160` Nitrogen metabolism (L0) |
| go_bp | 14 | `go:0009399` nitrogen fixation (L4), `go:0019740` nitrogen utilization (L3) |
| go_mf | 12 | mostly enzyme-activity terms |
| kegg | 12 | `kegg.pathway:ko00910` Nitrogen metabolism (L2), `ko01310` Nitrogen cycle (L2) |
| ec | 13 | N-related enzyme classes (6.7.-.- etc.) |
| pfam | 2 | `pfam:PF00795` Carbon-nitrogen hydrolase, `pfam:PF00543` P-II regulatory protein |

## Organism-ontology compatibility matrix

| Ontology | MED4 | MIT9313 | SS120 |
|---|---|---|---|
| **cyanorak_role** (selected) | ✓ 73% | ✓ 55% | ✓ 68% |
| **kegg** (selected) | ✓ 38% | ✓ 28% | ✓ 15% |
| tigr_role | ✓ 87% | ✓ 71% | ✓ 83% |
| go_bp | ✓ 54% | ✓ 43% | ✓ 55% |

All 3 organisms are cyanobacteria, so both selected ontologies cover the full set. (Step 1a originally included HOT1A3 as CTX1, which has no cyanorak/tigr annotation — dropped during this review; see `exploration/notebook.md` Step 1a redo entry and `data/experiments_classified.csv`.)

---

## Selected ontologies (2)

### Pick 1 — `cyanorak_role` level 1 (relevance-based, MED4-optimal, cyanobacteria-specific)

**Why this pick.** cyanorak_role is the cyanobacteria-specific functional-role ontology. Level 1 is the recommended depth by `relevance_rank` for MED4 (rank 2 overall, rank 1 for any non-tigr ontology). All 4 canonical N-limitation key-pathway anchors present at this level with MED4 gene counts above min-5:

- `cyanorak.role:E.4` **Nitrogen metabolism** (28 genes): `glnA`, `ntcA`, `amt1`, `cynA-S`, `ureA-D`, `glsF`
- `cyanorak.role:J.7` Photosystem I (17 genes): `psaA-M`, `ycf3/4/37`
- `cyanorak.role:J.8` Photosystem II (31 genes): `psbA-Z`, `ycf48`
- `cyanorak.role:J.2` CO2 fixation (22 genes): `rbcL/S/R`, carboxysome shell, Calvin cycle
- `cyanorak.role:K.2` Ribosomal proteins (61 genes): full `rpl*/rpm*/rps*` set
- `cyanorak.role:A.3` AA biosynthesis Glu family (16 genes): `argB-J`, `glnA`, `proA-S`

**Coverage:** MED4 73%, MIT9313 55%, SS120 68%. All remaining CTX organisms covered.

### Pick 2 — `kegg` (pathway) level 2 (complementary pathway-map framing)

**Why this pick.** KEGG pathway is an orthology-based pathway-map ontology (not a DAG). Level 2 is the pathway level (between category at L1 and KO at L3). `relevance_rank` is lower (32) because genome coverage is only 38%, but the median gene-set size (12) sits in the sweet spot and the N anchor is clean:

- `kegg.pathway:ko00910` **Nitrogen metabolism** (6 genes): `glnA`, `cynA-S`, `glsF` (small but direct)
- `kegg.pathway:ko00195` Photosynthesis (51 genes): combined PS I + PS II + cyt b6f + ATP synthase
- `kegg.pathway:ko03010` Ribosome (54 genes)
- `kegg.pathway:ko00250` Ala/Asp/Glu metabolism (15 genes): GS-GOGAT cycle
- `kegg.pathway:ko00260` Gly/Ser/Thr metabolism (25 genes)

**Coverage:** MED4 38%, MIT9313 28%, SS120 15% (lowest — expected due to SS120's reduced KEGG annotation). All organisms covered at min-5 for the N/photosynthesis/ribosome anchors.

**Complementarity with Pick 1.** Different framing: cyanorak is functional role ("what does this gene do?"); KEGG is metabolic pathway ("what network is this gene in?"). Independent categorizations → cross-ontology agreement is an informative stability check (spec §5 Step 4 M3 check).

---

## Considered and rejected alternatives

### `tigr_role` level 0 (rank 1 by coverage)

Relevance rank 1 for MED4 (87% coverage, 77 terms, median 13). Rejected because:
- Only 1 N-specific term (`tigr.role:160` Nitrogen metabolism) vs cyanorak's more granular categorization (E.4 core + D.1.3 adaptation + A.3 Glu family AA biosynthesis).
- Photosynthesis collapses to a single term rather than PSI/PSII/CO2 fixation as separate anchors — less diagnostic for the expected-direction QC.
- TIGR is maintained, but cyanorak's cyanobacteria-specific hierarchy is more aligned with MED4 biology.

### `go_bp` level 3 (rank 5)

Considered as an "organism-universal" Pick 2 when HOT1A3 was still in the CTX set. After dropping CTX1, the organism-universality requirement disappeared.

Evaluated with corrected MED4 annotation data: GO BP is annotation-sparse for MED4 N-limitation biology. Only **2 N-specific terms pass min-5 gene-set size** across all go_bp levels:

| GO BP term | Level | MED4 genes |
|---|---|---|
| `go:0071941` nitrogen cycle metabolic process | L3 | 13 |
| `go:0071705` nitrogen compound transport | L4 | 18 |
| `go:0006995` **cellular response to nitrogen starvation** | L5 | <5 (ideal term, but MED4-sparse) |
| `go:0009399` nitrogen fixation | L4 | 0 |
| `go:0043562` cellular response to nitrogen levels | L4 | <5 |
| `go:0019740` nitrogen utilization | L3 | 2 |
| ... (all other N-labeled terms) | L3–L6 | <5 |

Going deeper doesn't surface more N-specific terms — MED4's GO BP is thin on N-limitation specifics regardless of level.

### `go_bp` level 3–5 union via `term_ids` (considered late)

`pathway_enrichment` accepts a hand-curated `term_ids` list, so we could enrich against a cross-level GO BP panel. Rejected because:

1. **No new N-specific signal.** The two N-terms passing min-5 (`go:0071941` L3 + `go:0071705` L4) are already accessible at a single level; the rest of the N-labeled panel is annotation-sparse in MED4.
2. **DAG parent-child redundancy.** GO is a DAG; genes annotated to `go:0006995` (L5) also roll up to its parents. Enriching against L3∪L4∪L5 counts the same genes in multiple overlapping terms, inflating BH-corrected tests (343 non-independent tests per cluster instead of 68 L3-peer tests).
3. **Fisher null ambiguity.** At `level=3` the null is "is this pathway enriched among the 68 L3 processes that partition the genome?" — clean. At `term_ids=[L3∪L4∪L5]` the null becomes "is this pathway enriched among these 343 overlapping terms?" — parent-child overlap makes the background distribution of "typical" enrichment harder to interpret.
4. **Curated `term_ids` = our key-pathway panel by another name.** If we manually pick ~20 GO terms covering the biology, we've effectively hand-constructed a panel — losing the unbiased-enrichment property that motivated a second ontology.

Captured in `gaps_and_friction.md` as feature-scope thoughts for a future "DAG-aware pathway_enrichment" mode (parent-child deduplication, subtree-aware level selection, annotation-sparsity handling).

### BRITE subtrees (transporters, enzymes, ribosome, …)

Considered as potential Pick 3. Rejected because:
- BRITE subtrees are small (3–34 terms at L1, 4–32% genome coverage for enzymes, <5% for others).
- `transporters` is the most N-relevant (has `Phosphate and amino acid transporters`), but mixes P+N and covers only 4% of MED4 genome.
- `enzymes` L1 is usable (32%, 34 terms) but largely redundant with KEGG pathway membership.
- BRITE trees each require separate `pathway_enrichment` calls with `tree=<name>` — adds analysis complexity without a strong unique-signal story for B2.

BRITE drill-down, if needed, is better done as a targeted `explore_*.py` script on specific trees after Step 2.

### 3 ontologies total

Considered adding BRITE transporters or tigr_role as Pick 3. Rejected because:
- Cross-ontology agreement (spec §5 Step 4 M3) is more informative with 2 well-framed orthogonal ontologies than 3 overlapping ones.
- Fig 1 legibility (spec §8 risk concerns) is better with fewer ontology panels.
- Additional picks would mostly reinforce signal already visible in cyanorak_role + kegg rather than add independent evidence.

---

## Key-pathway panel

See `exploration/key_pathways.csv` for the full 11-row panel. Summary:

| ontology | term_id | term_name | expected | marker gene | n_MED4_genes |
|---|---|---|---|---|---|
| cyanorak_role | `cyanorak.role:E.4` | Nitrogen metabolism | up | glnA | 28 |
| cyanorak_role | `cyanorak.role:J.7` | Photosystem I | down | psaA | 17 |
| cyanorak_role | `cyanorak.role:J.8` | Photosystem II | down | psbA | 31 |
| cyanorak_role | `cyanorak.role:J.2` | CO2 fixation | down | rbcL | 22 |
| cyanorak_role | `cyanorak.role:K.2` | Ribosomal proteins | down | rplA | 61 |
| cyanorak_role | `cyanorak.role:A.3` | AA biosynthesis Glu family | up | glnA | 16 |
| kegg | `kegg.pathway:ko00910` | Nitrogen metabolism | up | glnA | 6 |
| kegg | `kegg.pathway:ko00195` | Photosynthesis | down | psbA | 51 |
| kegg | `kegg.pathway:ko03010` | Ribosome | down | rplA | 54 |
| kegg | `kegg.pathway:ko00250` | Ala/Asp/Glu metabolism | up | glnA | 15 |
| kegg | `kegg.pathway:ko00260` | Gly/Ser/Thr metabolism | ambiguous | glyA | 25 |

All term IDs validated via `genes_by_ontology(term_ids=..., organism="MED4")` — canonical marker genes (`glnA`, `cynA-S`, `psbA`, `psaA`, `rbcL`, `rplA`) confirmed present in each term's MED4 gene set. `ko00260` direction left ambiguous — Gly/Ser/Thr pathway includes both N-scavenging precursors and protein-biosynthesis sinks, so sign is data-driven rather than prior-constrained.

The key-pathway panel is the **biological anchor for Step 2 QC** (spec §5 Step 2 explore): R clusters must show UP enrichment at N-metabolism and DOWN enrichment at photosynthesis/ribosome. If any canonical marker fails to appear in the expected direction, the signal/ontology combination is wrong and Step 1b gets redone.
