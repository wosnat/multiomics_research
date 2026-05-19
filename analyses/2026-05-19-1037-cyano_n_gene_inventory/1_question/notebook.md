# Step 1 — Research question

## Context

User prompt: *"I want an overview of N-related genes in the KG. cyano only."*

Two upstream forks were resolved before scope-locking:

1. **Quick browse vs new analysis.** Researcher chose to start a formal 6-step analysis. The active MED4 axenic up-hypotheticals dossier (`analyses/2026-04-29-1025-axenic_up_hypotheticals_med4/`, end of step 5) is paused.
2. **Angle of "overview".** Among (annotation landscape only / behavior under N perturbation / behavior across all conditions / cross-strain comparison / other), researcher chose **cross-strain comparison** — how does the N-related gene set differ across HL/LL *Prochlorococcus* clades and *Synechococcus*.

A third fork emerged during step-1 grounding: how to handle cyano strains lacking Cyanorak annotation. Resolved by switching to an **orthology-bridged** gene-set definition (Cyanorak anchor → eggnog Cyanobacteria-level ortholog group expansion), which extends coverage from the 13 Cyanorak-annotated strains to all 19 cyano genome_strains in the KG.

## Locked question

> **What N-related genes — anchored by Cyanorak roles `E.4 ∪ D.1.3` across the Cyanorak-annotated strains and expanded via their eggnog Cyanobacteria-level ortholog groups — exist across the 19 cyanobacterial genome_strains in the KG, and how does the inventory differ across *Prochlorococcus* clades (HL vs LL) and across *Synechococcus* lineages?**

### In scope (this analysis)

**19 cyano genome_strains.** Clade assignments below come from the **Cyanorak source CSVs**, not the KG `clade` field (two KG bugs documented in `gaps_and_friction.md`: Synechococcus clades are null, and NATL1A/NATL2A clades are mislabeled).

- ***Prochlorococcus* (11 strains, all in Cyanorak Pro table except RSP50):**
  - HLI: MED4, RSP50 (RSP50 not in Cyanorak Pro table → clade from KG)
  - HLII: AS9601, MIT9301, MIT9312
  - LLI: MIT0801, NATL1A, NATL2A
  - LLII: SS120
  - LLIV: MIT9303, MIT9313
- ***Synechococcus* marine (4 strains in Cyanorak Syn table, SubCluster 5.1):**
  - WH8102 — Clade III, SubClade IIIa, Pigment 3c
  - WH7803 — Clade V, SubClade V, Pigment 3a
  - CC9311 — Clade I, SubClade Ia, Pigment 3dA
  - BL107 — Clade IV, SubClade IVa, Pigment 3dA
- **Other cyano (4 strains, outside Cyanorak picocyano scope):**
  - Synechococcus PCC 7002
  - Synechococcus elongatus PCC 7942
  - Synechococcus elongatus UTEX 2973
  - Thermosynechococcus vestitus BP-1

**Gene-set definition (two-step):**

1. **Anchor set.** Cyanorak roles `cyanorak.role:E.4 ∪ cyanorak.role:D.1.3` across all 13 Cyanorak-annotated strains. Per-strain anchor sizes range from 12 (SS120) to 43 (WH8102).
2. **Ortholog expansion.** For each anchor gene, take its eggnog ortholog group at `taxonomic_level='Cyanobacteria'` (specificity_rank=2). Union of group IDs across the anchor set = the analysis's gene-set, expressed as ortholog groups rather than per-strain locus tags.
3. **Cross-strain inventory.** For each group, list members across all 19 strains → 19 × N(groups) presence/copy-number matrix.

### Deferred to step 2 (KG entries)
- Build the anchor set from the 13 Cyanorak-annotated strains and resolve each anchor gene to its eggnog Cyanobacteria-level group.
- Enrich the strain table with Cyanorak source-CSV fields (SubCluster, Clade, SubClade, Pigment for marine Syn; species naming consistency).
- Identify the 4 N-perturbation publications in the KG as Background context for `paper.md` (Tolonen 2006, Read 2017, Domínguez-Martín 2017, Weissberg 2025) — surfaced even though their experiments are not analyzed in this inventory.

### Deferred to step 3 (framing)
- *How* to organize the comparison: presence/absence matrix only, copy-number matrix, per-Cyanorak-category breakdown, or some combination.
- Whether to layer on additional ortholog tiers (curated/Prochloraceae/Bacteria) for genes where the Cyanobacteria-level group splits an ecologically meaningful lineage (e.g., Prochlorococcus vs Synechococcus amt1 are in separate Cyanobacteria-level groups).

### Out of scope (not this analysis)
- Expression behavior across experiments.
- Functional response of N-machinery under N perturbation.
- Predictions about N preference of strains in the field (interpretive — Discussion-only if warranted at step 6).
- Habitat / ecology labels for cyano strains beyond what the Cyanorak source CSVs provide.

## KG context

Grounding queries run before the question was locked. Counts and structural notes that informed scope decisions.

### Cyano strains available in KG (19 total, genome_strain type)

`list_organisms` returned 37 organisms across all types; manually filtered to cyanobacterial genome_strains. **Pagination matters** — default `limit=5` and even my `limit=25` first call missed 7 cyano strains on page 2; I had to paginate to `offset=25` to find them.

| Strain | Genus | Genes | Pubs | Experiments | Clade (Cyanorak source) |
|---|---|---|---|---|---|
| AS9601 | Prochlorococcus | 1951 | 1 | 1 | HLII |
| MED4 | Prochlorococcus | 1976 | 17 | 114 | HLI |
| MIT0801 | Prochlorococcus | 2358 | 1 | 6 | LLI |
| MIT9301 | Prochlorococcus | 1935 | 3 | 8 | HLII |
| MIT9303 | Prochlorococcus | 3114 | 1 | 4 | LLIV |
| MIT9312 | Prochlorococcus | 1978 | 5 | 27 | HLII |
| MIT9313 | Prochlorococcus | 2948 | 9 | 38 | LLIV |
| NATL1A | Prochlorococcus | 2226 | 1 | 2 | LLI *(KG bug: KG says LLII)* |
| NATL2A | Prochlorococcus | 2214 | 6 | 17 | LLI *(KG bug: KG says LLII)* |
| RSP50 | Prochlorococcus | 1870 | 0 | 0 | HLI (KG only; not in Cyanorak Pro CSV) |
| SS120 | Prochlorococcus | 1964 | 3 | 64 | LLII |
| WH8102 | Synechococcus | 2881 | 4 | 74 | 5.1 / III / IIIa / 3c |
| WH7803 | Synechococcus | 2621 | 4 | 69 | 5.1 / V / V / 3a |
| CC9311 | Synechococcus | 3052 | 1 | 11 | 5.1 / I / Ia / 3dA |
| BL107 | Synechococcus | 2595 | 1 | 60 | 5.1 / IV / IVa / 3dA |
| PCC 7002 | Synechococcus | 3207 | 1 | 6 | not in Cyanorak picocyano CSVs |
| PCC 7942 | Synechococcus elongatus | 2717 | 1 | 4 | not in Cyanorak picocyano CSVs |
| UTEX 2973 | Synechococcus elongatus | 2714 | 1 | 2 | not in Cyanorak picocyano CSVs |
| BP-1 | Thermosynechococcus vestitus | 2449 | 1 | 1 | not in Cyanorak picocyano CSVs |

### N-perturbation publications in KG

`list_publications(search_text="nitrogen")` returned 8 hits; 4 have `treatment_type='nitrogen'`:

| DOI | Year | Lead | Organism(s) | Omics | KG experiments |
|---|---|---|---|---|---|
| 10.1038/msb4100087 | 2006 | Tolonen | MED4, MIT9313 | MICROARRAY | 6 (time-course, N depletion) |
| 10.1038/ismej.2017.88 | 2017 | Read | MED4 | RNASEQ | 1 (N-deprived TSS mapping) |
| 10.1128/mSystems.00008-17 | 2017 | Domínguez-Martín | SS120 | PROTEOMICS | 1 (N-limited) |
| 10.1101/2025.11.24.690089 | 2025 | Weissberg | MED4, HOT1A3 | RNASEQ, PROTEOMICS | 10 (coculture, N recycling) |

These are not the target of this inventory analysis but form **Background** for `paper.md`.

### Candidate N-gene definitions sized on MED4

`genes_by_ontology` and `genes_by_function` were run for each candidate. Distinct-gene counts in MED4:

| Definition | MED4 gene count | Notes |
|---|---|---|
| KEGG `ko00910` "Nitrogen metabolism" | 6 | cynABDS, glnA, glsF only. Misses NtcA, GlnB, urease, urea transport, AMT1. Too narrow — see `gaps_and_friction.md`. |
| KEGG `ko01310` "Nitrogen cycle" | 0 | Filtered out for MED4. |
| KEGG N-regulator KOs (8 queried) | 2 unique | NtcA + GlnB only — the rest aren't in MED4. |
| `genes_by_function` "nitrogen" | 35 | Real machinery + false-positives ("apolipoprotein **N**-acyltransferase"). Lucene noise. |
| GO BP (13 N-terms union) | 35 | Real set + noisy edges (secA/B/D, ccsA, bioY tagged as "nitrogen compound transport"). |
| **Cyanorak `E.4 ∪ D.1.3`** | **28** | Clean cyano-specific curation. Selected as anchor for ortholog expansion. |
| TIGR `role:160` | 7 | Strict subset; too narrow. |

### Cross-strain sanity check on Cyanorak union

| Strain | Cyanorak `E.4 ∪ D.1.3` count |
|---|---|
| MED4 (HLI) | 28 |
| AS9601 (HLII) | 24 |
| MIT9301 (HLII) | 25 |
| MIT9312 (HLII) | 24 |
| MIT9303 (LLIV) | 26 |
| MIT9313 (LLIV) | 26 |
| NATL1A (LLI) | 27 |
| NATL2A (LLI) | 27 |
| SS120 (LLII) | 12 *(flag: unusually low, investigate at step 2)* |
| WH8102 (Syn) | 43 |
| WH7803 (Syn) | 27 *(flag: lower than other marine Syn, investigate at step 2)* |
| CC9311 (Syn) | 42 |
| BL107 (Syn) | 37 |
| MIT0801, RSP50, PCC 7002, PCC 7942, UTEX 2973, T. vestitus BP-1 | **0 / no Cyanorak annotation** |

Six strains lack Cyanorak annotation entirely. This motivates the ortholog bridge.

### Ortholog-bridge feasibility check (4 representative N-genes via eggnog Cyanobacteria-level groups)

Probed `gene_homologs` on PMM0246 (ntcA), PMM0920 (glnA), PMM0963 (ureC), PMM0263 (amt1). Each gene maps to 4 tiers of ortholog group (cyanorak curated → eggnog Prochloraceae → eggnog Cyanobacteria → eggnog Bacteria). The **eggnog Cyanobacteria tier** (`taxonomic_level='Cyanobacteria'`, `specificity_rank=2`) is the bridge — it includes the 6 strains that lack Cyanorak.

| Gene | eggnog Cyanobacteria group | Members | Strains covered (of 19) |
|---|---|---|---|
| NtcA | `eggnog:1G07U@1117` | 20 | **all 19** including MIT0801, RSP50, PCC 7002, PCC 7942, UTEX 2973, T. vestitus BP-1 |
| GlnA | `eggnog:1G255@1117` | 19 | **all 19** |
| UreC | `eggnog:1G12D@1117` | 15 | 15/19 — missing in WH7803, PCC 7942, UTEX 2973, plus one not in the returned set (biologically plausible: freshwater *S. elongatus* often lacks urease) |
| Amt1 | `eggnog:1GD8R@1117` | 9 | 9 Prochlorococcus only — **Synechococcus amt1 is in a separate Cyanobacteria-level group**; must seed from both Pro and Syn anchors to capture lineage-specific paralogs |

The ortholog bridge works. Two nuances:
- **Lineage-specific paralog groups exist even at Cyanobacteria level** (Pro amt1 vs Syn amt1). Seeding from both Cyanorak Pro and Cyanorak Syn anchors handles this.
- **Some absences are biologically real** (no urease in PCC 7942 / UTEX 2973). The ortholog bridge separates real biology from annotation gaps.

## Decisions

**2026-05-19 — Gene-set definition.** Cyanorak `E.4 ∪ D.1.3` as the **anchor set**, expanded via eggnog ortholog groups at `taxonomic_level='Cyanobacteria'` (specificity_rank=2). KEGG `ko00910` rejected for narrowness; GO BP / Lucene rejected for noise; pure Cyanorak rejected because it covers only 13 of 19 cyano strains in the KG.

**2026-05-19 — Strain set: all 19 cyano genome_strains in the KG.** The 6 strains lacking Cyanorak annotation (MIT0801, RSP50, PCC 7002, PCC 7942, UTEX 2973, T. vestitus BP-1) are brought into scope via the ortholog bridge. SS120 and WH7803 are kept but flagged for coverage investigation at step 2.

**2026-05-19 — Clade assignments use Cyanorak source CSVs, not the KG `clade` field.** Two KG bugs (Syn clade null; NATL1A/NATL2A mislabeled as LLII) make the KG `clade` field unreliable here. Both bugs are logged in `gaps_and_friction.md` with the loader-fix recommendation.

**2026-05-19 — Comparison axes deferred to step 3.** The presence/absence vs copy-number framing, per-Cyanorak-category breakdown, and use of additional ortholog tiers are step-3 framing decisions made after step 2 surfaces the actual group inventory.

**2026-05-19 — MED4 axenic up-hypotheticals analysis paused.** That analysis is at end of step 5; resuming requires reading `analyses/2026-04-29-1025-axenic_up_hypotheticals_med4/` and `memory/project_med4_axenic_up_hypotheticals.md`.

## Surprises

- **KEGG `ko00910` is much narrower than its name implies** for Prochlorococcus — 6 genes in MED4, missing NtcA, urease, urea transport. Logged in `gaps_and_friction.md`.
- **D.1.3 ⊆ E.4 in MED4 but not in larger genomes.** In MED4 all 16 D.1.3 genes are in E.4. In WH8102 the two sets contribute different genes (D.1.3 covers nitrate/nitrite + Mo cofactor; E.4 covers cyanate + Met biosynthesis). Union is the right choice.
- **Met-biosynthesis genes are in Cyanorak E.4** (cysK, metB, MetB-like, metY, metC). They *use* N rather than handle the N cycle proper. Curators' choice; flagged so it's not surprising in step 3.
- **Ortholog groups split lineage-specific paralogs even at Cyanobacteria level** (Pro amt1 vs Syn amt1). Method handles by seeding from both Cyanorak Pro and Syn anchors.
- **Two KG data-quality bugs surfaced during scoping** — Syn clade is null in KG (source CSV has it); NATL1A/NATL2A are mislabeled as LLII in KG (source CSV says LLI). Both in `gaps_and_friction.md`.
- **One anti-hallucination slip caught by the researcher** — I labeled Synechococcus habitats from intrinsic knowledge before checking the KG. Corrected by deferring to the Cyanorak source CSVs. Logged.

## Advance rationale

Question is locked. The anchor + ortholog-bridge gene-set definition is grounded against alternatives on MED4 and against 4 representative N-genes across all 19 strains. Strain set is enumerated with clade assignments from the authoritative source. KG bugs that affect this analysis are documented separately and don't block progress.

---

## Decide-gate checklist

- **Outputs produced.**
  - `analyses/2026-05-19-1037-cyano_n_gene_inventory/.gitignore` (template)
  - `analyses/2026-05-19-1037-cyano_n_gene_inventory/paper.md` (skeleton; Question populated)
  - `analyses/2026-05-19-1037-cyano_n_gene_inventory/gaps_and_friction.md` (4 entries: KEGG narrowness, anti-hallucination slip, Syn-clade KG bug, NATL clade KG bug)
  - `analyses/2026-05-19-1037-cyano_n_gene_inventory/1_question/notebook.md` (this file)
  - No scripts, data, or figures (step 1 is a conversation).
- **Results presented.** All grounding tables (cyano organisms, N publications, candidate gene-set comparison, Cyanorak per-strain counts, ortholog-bridge feasibility on 4 genes) shown inline above as the same tables presented in chat.
- **QC gate.**
  - Candidate N-gene definitions sized against each other on MED4 → Cyanorak anchor chosen.
  - Cyanorak coverage checked across all 19 cyano → 13 covered, 6 not.
  - Ortholog-bridge feasibility checked on 4 representative N-genes (ntcA, glnA, ureC, amt1) → bridge works for the 6 uncovered strains; lineage-specific paralog quirk identified and handled by seeding from both Pro and Syn anchors.
  - Clade assignments cross-checked against the Cyanorak source CSVs → 2 KG bugs identified and logged.
- **Decisions made this step.** Five (logged above): gene-set anchor + ortholog expansion; 19-strain scope; clade assignments from source CSVs; comparison axes deferred to step 3; MED4 axenic up-hypotheticals analysis paused.
- **Advance rationale.** Locked question, grounded scope, explicit deferrals, two KG bugs documented but not blocking. Step 2 (KG entries) can begin — its first task is to build the anchor set across the 13 Cyanorak-annotated strains and resolve each anchor gene to its eggnog Cyanobacteria-level group(s), producing the ortholog-group inventory matrix.
