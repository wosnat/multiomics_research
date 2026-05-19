# Step 2 — KG entries

## Context

Step 1 locked the question and method: anchor gene set = Cyanorak `E.4 ∪ D.1.3` across 13 Cyanorak-annotated strains, expanded via eggnog Cyanobacteria-level ortholog groups to cover all 19 cyano strains. Step 2 builds the concrete KG entries this analysis will operate on:

- **Anchor gene set** — the 375 distinct locus_tags from the 13 Cyanorak strains that will seed the ortholog expansion.
- **Strain table** — 19 cyano strains with KG facts + Cyanorak source-CSV enrichment for clade/subcluster/pigment.
- **Background publications** — 4 N-perturbation studies in the KG (cited in paper.md Background; not analyzed here).

Step 2 also resolves the two outlier flags from step 1 (SS120 = 12 anchors; WH7803 = 27 anchors) — both turn out to be real biological gene losses, not annotation artifacts.

## What I did

```bash
uv run python 2_kg_selection/scripts/01_anchor_genes.py
uv run python 2_kg_selection/scripts/02_strain_table.py
uv run python 2_kg_selection/scripts/03_n_publications.py
uv run python 2_kg_selection/scripts/qc_outlier_coverage.py
```

All scripts use the Python API (`from multiomics_explorer import ...`) and write CSVs + log files to `2_kg_selection/data/`.

## Results

### Anchor gene set per strain (Cyanorak `E.4 ∪ D.1.3`)

`data/01_anchor_summary.csv`, 13 rows.

| Strain | Anchor rows | Unique locus_tags | Unique gene names |
|---|---|---|---|
| Prochlorococcus AS9601 (HLII) | 40 | 24 | 23 |
| Prochlorococcus MED4 (HLI) | 44 | 28 | 27 |
| Prochlorococcus MIT9301 (HLII) | 42 | 25 | 24 |
| Prochlorococcus MIT9303 (LLIV) | 45 | 27 | 26 |
| Prochlorococcus MIT9312 (HLII) | 40 | 24 | 23 |
| Prochlorococcus MIT9313 (LLIV) | 45 | 27 | 26 |
| Prochlorococcus NATL1A (LLI) | 46 | 28 | 26 |
| Prochlorococcus NATL2A (LLI) | 46 | 28 | 28 |
| Prochlorococcus SS120 (LLII) | **16** | **12** | **12** |
| Synechococcus WH8102 (5.1/III) | 71 | 43 | 40 |
| Synechococcus WH7803 (5.1/V) | **43** | **28** | **27** |
| Synechococcus CC9311 (5.1/I) | 71 | 43 | 40 |
| Synechococcus sp. BL107 (5.1/IV) | 63 | 38 | 37 |

- 612 anchor rows total across 13 strains.
- 375 unique locus_tags across the union.
- Two outliers (SS120, WH7803) investigated below — both biological, not annotation.

Full per-gene table: `data/01_anchor_genes.csv` (612 rows; columns `strain, organism_name_kg, locus_tag, gene_name, product, gene_category, term_id`).

### 19-strain table with Cyanorak source-CSV enrichment

`data/02_strain_table.csv`, 19 rows. Selected columns:

| Strain | Genus | Genes | KG clade | Cyanorak clade | Cyanorak subcluster / subclade / pigment | In Cyanorak picocyano CSV |
|---|---|---|---|---|---|---|
| Prochlorococcus AS9601 | Prochlorococcus | 1951 | HLII | HLII | HL / HLII / Lba | ✓ |
| Prochlorococcus MED4 | Prochlorococcus | 1976 | HLI | HLI | HL / HLI / Lba | ✓ |
| Prochlorococcus MIT0801 | Prochlorococcus | 2358 | LLI | LLI | LL / LLI / Hba | ✓ |
| Prochlorococcus MIT9301 | Prochlorococcus | 1935 | HLII | HLII | HL / HLII / Lba | ✓ |
| Prochlorococcus MIT9303 | Prochlorococcus | 3114 | LLIV | LLIV | LL / LLIV / Hba | ✓ |
| Prochlorococcus MIT9312 | Prochlorococcus | 1978 | HLII | HLII | HL / HLII / Lba | ✓ |
| Prochlorococcus MIT9313 | Prochlorococcus | 2948 | LLIV | LLIV | LL / LLIV / Hba | ✓ |
| Prochlorococcus NATL1A | Prochlorococcus | 2226 | **LLII** *(KG bug)* | **LLI** | LL / LLI / Hba | ✓ |
| Prochlorococcus NATL2A | Prochlorococcus | 2214 | **LLII** *(KG bug)* | **LLI** | LL / LLI / Hba | ✓ |
| Prochlorococcus RSP50 | Prochlorococcus | 1870 | HLI | — | — | ✗ |
| Prochlorococcus SS120 | Prochlorococcus | 1964 | LLII | LLII | LL / LLII / Hba | ✓ |
| Synechococcus WH8102 | Synechococcus | 2881 | **null** *(KG bug)* | III | 5.1 / IIIa / 3c | ✓ |
| Synechococcus WH7803 | Synechococcus | 2621 | **null** *(KG bug)* | V | 5.1 / V / 3a | ✓ |
| Synechococcus CC9311 | Synechococcus | 3052 | **null** *(KG bug)* | I | 5.1 / Ia / 3dA | ✓ |
| Synechococcus sp. BL107 | Synechococcus | 2595 | **null** *(KG bug)* | IV | 5.1 / IVa / 3dA | ✓ |
| Synechococcus PCC 7002 | Synechococcus | 3207 | null | — | — | ✗ |
| Synechococcus elongatus PCC 7942 | Synechococcus elongatus | 2717 | null | — | — | ✗ |
| Synechococcus elongatus UTEX 2973 | Synechococcus elongatus | 2714 | null | — | — | ✗ |
| Thermosynechococcus vestitus BP-1 | Thermosynechococcus | 2449 | null | — | — | ✗ |

Confirmed in script log: 2 strains with KG/Cyanorak clade mismatch (NATL1A, NATL2A), 4 strains with null KG clade where Cyanorak has it (the 4 marine Syn), 5 strains outside the Cyanorak picocyano CSVs (RSP50, PCC 7002, PCC 7942, UTEX 2973, T. vestitus BP-1). Both KG bugs were already logged in `gaps_and_friction.md` at step 1; this script reproduced them programmatically.

### 4 N-perturbation publications (Background context)

`data/03_n_publications.csv`, 4 rows.

| Year | First author | DOI | Organism(s) | Omics | Experiments |
|---|---|---|---|---|---|
| 2006 | Tolonen | 10.1038/msb4100087 | MED4, MIT9313 | MICROARRAY | 6 (time-course, N depletion) |
| 2017 | Read | 10.1038/ismej.2017.88 | MED4 | RNASEQ | 1 (N-deprived TSS mapping) |
| 2017 | Domínguez-Martín | 10.1128/mSystems.00008-17 | SS120 | PROTEOMICS | 1 (N-limited) |
| 2025 | Weissberg | 10.1101/2025.11.24.690089 | MED4, HOT1A3 | PROTEOMICS, RNASEQ | 10 (coculture, N recycling) |

These supply field-relevant Background for `paper.md` — they are *not* analyzed by this inventory analysis.

### Outlier QC — SS120 and WH7803 are real biological gene losses

`data/qc_outlier_coverage.csv`. For each outlier, list anchor genes missing relative to the peer consensus (gene name in ≥half of same-lineage peers), and check whether the gene exists in the outlier's proteome under a different annotation via `genes_by_function(search_text=gene_name, organism=outlier)`.

**SS120** (LLII Prochlorococcus): 14 anchor genes missing relative to consensus; **0 recovered under different annotation**. All 14 are truly absent. Composition:

| Gene | Peers with it | Outlier present? |
|---|---|---|
| ureA, ureB, ureC, ureD, ureE, ureF, ureG | 8/8 | no — entire urease operon absent |
| urtA, urtB, urtC, urtD, urtE | 8/8 | no — entire urea ABC transporter absent |
| nirA | 4/8 | no |
| focA | 4/8 | no |

The 12 SS120 anchor genes are: ntcA, glnB, glsF, glnA, amt1, cynABDS (cyanate), cysK, metB, metY, MetB-like, metC, carA. SS120 keeps ammonium uptake, cyanate utilization, central N assimilation, and Met biosynthesis — but loses urea uptake/hydrolysis entirely.

**WH7803** (marine Syn): 17 anchor genes missing relative to peers; **0 recovered under different annotation**. Composition:

| Gene | Peers with it | Outlier present? |
|---|---|---|
| ureA, ureB, ureC, ureD, ureE, ureF, ureG | 3/3 | no — entire urease operon absent |
| urtB, urtC, urtD, urtE | 3/3 | no — most of urea transport absent (urtA appears present) |
| amt2, cynA, cynB, cynD, merR, ntrX | 1/3 | no — variable across marine Syn |

WH7803's N-uptake repertoire (per anchor): NtcA, GlnB, GlsF, GlnA, AMT1, nitrate/nitrite pathway (nrtP×2, narB, narM, nirA), plus Mo cofactor biosynthesis. **Lacks urease + most urea transport** — matches its low coverage.

**Interpretation [interpretation].** Both outliers reflect *gene-content reduction* relative to their lineage peers — the kind of cross-strain variation this inventory analysis is meant to surface. The low Cyanorak anchor counts at step 1 were misleading as a "data quality concern" — they are the actual biology. Genome streamlining in SS120 (1964 genes, smallest in the LL group) and the absence of urea-cycle machinery in WH7803 are both well-documented in the literature for these strains, but here we verified them directly from the KG with no a-priori claim.

## Surprises

- **SS120 lacks the entire urea pathway** — both urtABCDE (uptake) and ureABCDEFG (hydrolysis). 12 of its 14 "missing" anchor genes are these two operons.
- **WH7803 also lacks urease** (and most urea transport) — same biological pattern as SS120 but in marine *Synechococcus*. Keeps nitrate/nitrite pathway and urtA. Convergent loss of urease across lineages.
- **PCC 7942 and UTEX 2973 have identical N-machinery present sets** — both 24 groups. Expected for closely-related *S. elongatus* lab strains; the inventory cleanly confirms it.
- **Thermosynechococcus BP-1 lacks the cysK/metB/metY Met-biosynthesis cluster** that's universal among picocyano in scope. Worth a Discussion note (step 6) — thermophile Met biosynthesis may use different enzymes.
- **8 of 9 orphan anchors are paralog copies whose gene-name is still represented elsewhere in the inventory.** The Cyanobacteria-level orthology splits lineage-specific paralogs (e.g., MED4's cynA orphan but other strains have cynA in a different group). The inventory captures gene names robustly but undercounts a few paralog copies.
- **Cyanorak NATL1A/NATL2A `Pigment type` = Hba, not Lba.** All HL Prochlorococcus in the Cyanorak Pro CSV are Lba; all LL Prochlorococcus are Hba. Adds a second piece of evidence for the NATL1A/2A loader bug being a systematic mismapping, not just a clade-field typo.

## Ortholog expansion (researcher-redirected mid-step extension)

After the four scripts above completed, the researcher requested pulling the
ortholog-bridge expansion forward into step 2, so the full gene-set inventory
exists before step 3 framing. Per just-in-time formalization, this is the
right call — frame the comparison once we can see what the data looks like.
(Methodology note logged in `gaps_and_friction.md`.)

```bash
uv run python 2_kg_selection/scripts/04_anchor_to_ortholog.py
uv run python 2_kg_selection/scripts/05_ortholog_inventory.py
uv run python 2_kg_selection/scripts/qc_ortholog_inventory.py
```

### Anchor → ortholog mapping (`04_anchor_to_ortholog.py`)

- Input: 375 unique anchor locus_tags.
- Output: 366 (97.6%) mapped to an eggnog Cyanobacteria-level group; **54 unique ortholog groups** result.
- 9 anchors have no Cyanobacteria-level group (and no Prochloraceae-level group either); only 1 of them (MED4 cynA / PMM0370) has a Bacteria-level group (nrtC). The other 8 are paralog copies whose lineage IS represented in the inventory via a different anchor's group (verified in `qc_ortholog_inventory.py`).
- Group breakdown: 35 cross-genus groups, 19 single-genus groups; median 8 members per group, max 21.

### 19-strain × 54-group inventory (`05_ortholog_inventory.py`)

Output: `05_inventory_matrix.csv` (wide, 19 × 54, cells = copy counts) + `05_inventory_members.csv` (long, 549 rows). Per-strain N(groups present) and total gene copies:

| Strain | N(groups present) | Total copies |
|---|---:|---:|
| Prochlorococcus AS9601 (HLII) | 24 | 24 |
| Prochlorococcus MED4 (HLI) | 27 | 27 |
| Prochlorococcus MIT0801 (LLI) | **28** | 28 |
| Prochlorococcus MIT9301 (HLII) | 25 | 25 |
| Prochlorococcus MIT9303 (LLIV) | 27 | 28 |
| Prochlorococcus MIT9312 (HLII) | 24 | 24 |
| Prochlorococcus MIT9313 (LLIV) | 26 | 26 |
| Prochlorococcus NATL1A (LLI) | 28 | 28 |
| Prochlorococcus NATL2A (LLI) | 28 | 29 |
| Prochlorococcus RSP50 (HLI) | **26** | 26 |
| Prochlorococcus SS120 (LLII) | **12** | 12 |
| Synechococcus WH8102 (5.1/III) | 40 | 43 |
| Synechococcus WH7803 (5.1/V) | 28 | 29 |
| Synechococcus CC9311 (5.1/I) | 40 | 41 |
| Synechococcus sp. BL107 (5.1/IV) | 38 | 39 |
| Synechococcus PCC 7002 | **35** | 39 |
| Synechococcus elongatus PCC 7942 | **24** | 26 |
| Synechococcus elongatus UTEX 2973 | **24** | 26 |
| Thermosynechococcus vestitus BP-1 | **28** | 29 |

**Bold = previously zero in pure Cyanorak (now covered via ortholog bridge).** The bridge added MIT0801 (28), RSP50 (26), PCC 7002 (35), PCC 7942 (24), UTEX 2973 (24), T. vestitus BP-1 (28) to the inventory.

Group ubiquity distribution (n strains where present): 7 groups universal (n=19), 2 in n=17, 1 in n=16, 8 in n=15, 3 in n=13, then a long tail down to 2 groups in n=1.

### QC results (`qc_ortholog_inventory.py`)

**Q1 — 7 universal (n=19) groups** = the canonical universally-conserved cyano N machinery:

| Group ID | Consensus name | Product |
|---|---|---|
| eggnog:1G07U@1117 | ntcA | global nitrogen regulatory protein |
| eggnog:1G255@1117 | glnA | glutamine synthetase, type I |
| eggnog:1G5QJ@1117 | glnB | nitrogen regulatory protein P-II |
| eggnog:1G0XM@1117 | glsF | ferredoxin-dependent glutamate synthase |
| eggnog:1G19V@1117 | carA | carbamoyl-phosphate synthase small subunit |
| eggnog:1G0FS@1117 | moeB | molybdopterin biosynthesis protein MoeB |
| eggnog:1G03T@1117 | metC | cystathionine beta-lyase |

All 7 cross-genus (Prochlorococcus + Synechococcus + Thermosynechococcus). NtcA + GlnB (regulation), GlnA + GlsF (assimilation), CarA (carbamoyl-P), MetC (Met biosynthesis), MoeB (Mo cofactor). This is what should be universally conserved; the inventory passes the sanity check.

**Q2/Q3/Q4 — Per-strain inventories for the outliers and previously-uncovered strains.** Full per-strain lists in the log. Highlights:

- **SS120** (12 present / 42 absent): keeps ntcA, glnA, glnB, glsF, carA, metC, metB, metY, cysK, moeB, amt1; **lacks all urea (urtABCDE, ureABCDEFG), all nitrate/nitrite (nrtP, narB, narM, nirA), all cyanate (cynABDS, cynH), most Mo cofactor, all paralog regulators (ntrB, ntrX)**. Minimum-machinery cyano in scope.
- **WH7803** (28 / 26): has full nitrate pathway, Mo cofactor, partial cyanate (cynS + cynH); **lacks urease + all urea transport except urtA**. Marine Syn that uses nitrate but not urea.
- **PCC 7002** (35 / 19): the richest non-marine cyano. Has nitrate, urea, urease, cyanate (partial), Mo cofactor, plus `glnN1` (alternative Gln synthetase) and 3 copies of amt1.
- **PCC 7942 = UTEX 2973** (24 each, identical present sets): both *S. elongatus*; both have nitrate pathway + cyanate (cynA, cynD, cynS) + glnN1; **both lack all urea pathway**. Identical biology — expected (they're closely-related lab strains).
- **T. vestitus BP-1** (28 / 26): has urea pathway, nitrate pathway, GlsF; **lacks cysK, metB, metY** (the Met-biosynthesis cluster that's universal among the picocyano). Thermosynechococcus shows different Met biosynthesis machinery.

**Q5 — Bridge integrity check:** all 549 member rows are in the 19-cyano scope. No Alteromonas / Marinobacter / Pseudomonas / etc. members leaked through the eggnog Cyanobacteria-level groups. Bridge is clean.

**Q6 — Paralogs:** 16 (strain × group) entries have multi-copy. Notable: PCC 7002 amt1 × 3 (the most extreme), WH8102 urtA × 2, multiple strains have cysK × 2 (CC9311, WH7803, WH8102, BL107), multiple have cynD × 2 (PCC 7002, WH8102, PCC 7942, UTEX 2973). Full table at `data/qc_paralogs.csv`.

**Q7 — Orphan anchors at higher tiers:** 0 of 9 orphans have a Prochloraceae-level group. Only 1 (PMM0370 cynA) has a Bacteria-level group (eggnog:COG0715@2 = nrtC). The other 8 (PMT2232 ureE; SYNW2462/SYNW2463 nrtP; sync_2895 moaB; sync_2280 amt2; sync_2840/sync_2903 cynH; BL107_06829 focA) have NO eggnog group at any tier.

**Per-orphan undercount audit (Q7b, post-hoc).** After the initial "orphans are paralog copies, don't matter" framing was challenged, I audited each orphan against the inventory: does any group in the inventory with the same `consensus_gene_name` actually include the orphan's strain as a member? Result:

| Orphan locus | Strain | Gene | Undercount before recovery |
|---|---|---|---|
| PMM0370 | MED4 | cynA | yes — MED4 absent in cynA group |
| PMT2232 | MIT9313 | ureE | yes — MIT9313 absent in both ureE groups |
| SYNW2462 + SYNW2463 | WH8102 | nrtP | yes — WH8102 absent in nrtP group (2 paralog copies unrepresented) |
| sync_2840 + sync_2903 | CC9311 | cynH | yes — CC9311 absent in cynH group (2 paralog copies unrepresented) |
| BL107_06829 | BL107 | focA | yes — BL107 absent in focA group |
| sync_2280 | CC9311 | amt2 | yes — amt2 was not in the inventory at all |
| sync_2895 | CC9311 | moaB | (already represented via another locus_tag, but kept as its own group for symmetry) |

**Recovery (`06_orphan_recovery.py`, researcher-redirected).** Earlier I framed these as a documented limitation to keep ortholog-tier purity. Researcher overrode: don't throw out the orphans. Each (strain, gene_name) orphan combination became one synthetic singleton group with `group_id = anchor_singleton:<gene>:<strain>`, `source = anchor_singleton`. 9 orphan locus_tags collapsed into 7 synthetic groups (WH8102 nrtP and CC9311 cynH each merge 2 paralog copies into one group with copy_count=2). The inventory matrix grew from 19 × 54 → **19 × 61**. Members CSV grew from 549 → 558 rows. Per-strain copy-count deltas after recovery: MED4 +1, MIT9313 +1, WH8102 +2, CC9311 +4, BL107 +1; net +9 copies — matches the orphan total. All other strains unchanged.

**Tier-mixing acknowledgment:** the inventory is now a mix of 54 eggnog Cyanobacteria-level groups (cross-strain orthology) and 7 anchor singleton groups (strain-specific, no orthology). The `source` column on `04_ortholog_groups.csv` distinguishes them. Downstream analysis should treat synthetic singletons as strain-specific "anchor presence" indicators, not as cross-strain comparable units.

## Decisions

**2026-05-19 — Two-tier strain stratification adopted.** All 19 strains are in scope for the ortholog-bridge inventory. Within the inventory output:
- 11 Cyanorak picocyano strains carry Cyanorak `Pigment` (Lba/Hba/3a/3c/3dA) as an additional annotation dimension.
- 8 strains outside Cyanorak picocyano (RSP50 + 4 non-marine Syn-family + T. vestitus) carry no pigment annotation; comparison will be presence/absence on the ortholog axis only.

**2026-05-19 — `clade_canonical` column added to `02_strain_table.csv`** = `cyanorak_clade` if present, else `clade_kg`. Downstream steps use `clade_canonical`. The KG `clade_kg` column is kept for traceability of the bugs.

**2026-05-19 — SS120 and WH7803 kept in scope.** Earlier they were flagged for coverage investigation; the investigation showed the low counts are real biology. No exclusion warranted.

**2026-05-19 — Ortholog tier locked at eggnog Cyanobacteria (specificity_rank=2).** The Cyanobacteria tier gives 366/375 anchor coverage with clean bridge integrity (Q5). Going to Prochloraceae would split lineage-specific paralogs further (smaller groups, less cross-comparability). Going to Bacteria would let in non-cyano members (Alteromonas, Marinobacter). Cyanobacteria is the right tier.

**2026-05-19 — Orphan anchors recovered as synthetic singleton groups (researcher-redirected reversal).** Earlier decision was to leave orphans out for ortholog-tier purity. Researcher overrode that — don't throw out data. Each (strain, gene_name) orphan combination became one synthetic singleton group via `06_orphan_recovery.py`. Inventory matrix grew 19 × 54 → 19 × 61. The synthetic-singleton convention is principled (a strain-specific anchor that has no Cyanobacteria-level eggnog group is still data); the trade-off is tier mixing, addressed by the `source` column on `04_ortholog_groups.csv` (`eggnog` vs `anchor_singleton`).

### Items deferred to later steps

- **Step 3 framing:** how to order rows (clade-aware) and columns (Cyanorak category breakdown); whether to use Cyanorak `Pigment` field as a row annotation; presence/absence vs copy-count encoding for the visualization.
- **Step 6 Discussion:** MIT9303 has ntcA × 2 — unusual Prochlorococcus duplication of the global N regulator, worth a brief interpretive note. PCC 7942 and UTEX 2973 have identical 24-group inventories — confirms the known near-isogeny of these *S. elongatus* lab strains.

## Advance rationale

KG entries enumerated, anchor set built, ortholog expansion performed, and 7-question QC suite passed: bridge integrity holds (no out-of-cyano leakage), the 7 universal groups are biologically sensible (canonical N machinery), the SS120/WH7803 outliers are explained as real biology, and orphan anchors are documented as a tolerated limitation. The 19 × 54 presence/copy-number matrix exists at `data/05_inventory_matrix.csv` and is the de-facto answer to the step-1 question; step 3 framing will decide how to present and order it (clade-aware ordering, category breakdown, heatmap visualization).

---

## Decide-gate checklist

- **Outputs produced.**
  - Scripts (run with `uv run python 2_kg_selection/scripts/<name>.py`): `01_anchor_genes.py`, `02_strain_table.py`, `03_n_publications.py`, `qc_outlier_coverage.py`, `04_anchor_to_ortholog.py`, `05_ortholog_inventory.py`, `qc_ortholog_inventory.py`, `06_orphan_recovery.py`.
  - Data: `01_anchor_genes.csv` (612 rows), `01_anchor_summary.csv` (13), `02_strain_table.csv` (19 × 17), `03_n_publications.csv` (4), `qc_outlier_coverage.csv` (31), `04_anchor_to_ortholog.csv` (366 rows), `04_anchors_without_group.csv` (9), `04_ortholog_groups.csv` (**61** post-recovery), `05_inventory_members.csv` (**558** post-recovery), `05_inventory_matrix.csv` (**19 × 61** post-recovery), `06_synthetic_groups.csv` (7), `qc_universal_groups.csv` (7), `qc_paralogs.csv` (16), `qc_orphan_anchors.csv` (9). Plus 8 .log files.
  - No figures (figures come in step 5).
- **Results presented.** Anchor summary, 19-strain enriched table, N-publications table, outlier QC, ortholog group inventory, per-strain group counts, ubiquity distribution, 7-question QC suite — all shown inline above with same numbers as the underlying CSVs.
- **QC gate (7 named checks).**
  - Q0a. Anchor extraction reproduces per-strain Cyanorak counts from step-1 grounding → match.
  - Q0b. 19-strain table reproduces both step-1-logged KG bugs (NATL1A/2A clade; null Syn clade) → confirmed programmatically.
  - Q0c. SS120 and WH7803 low anchor counts investigated → real biology (urea pathway loss), not annotation gaps.
  - Q1. 7 universal-coverage groups → all are canonical cyano N machinery (NtcA, GlnA, GlnB, GlsF, MetC, MoeB, CarA).
  - Q2/Q3/Q4. Per-strain group inventories for outliers + previously-uncovered strains → biologically interpretable; SS120/WH7803 absences match expectations; PCC 7942 = UTEX 2973 identical (expected).
  - Q5. Bridge integrity → all 549 member rows in-scope cyano; no Alteromonas/Marinobacter leakage.
  - Q6. Paralog inventory → 16 (strain × group) multi-copy entries identified for downstream copy-aware analysis.
  - Q7. Orphan anchors → 9 orphan loci, Q7b audit showed 6 strain × gene cells were undercounted; resolved by `06_orphan_recovery.py` adding 7 synthetic singleton groups (one per (strain, gene_name) orphan combination). Inventory matrix now 19 × 61.
- **Decisions made this step.** Five: two-tier stratification; `clade_canonical` convention; keep SS120/WH7803 in scope; lock ortholog tier at eggnog Cyanobacteria (for cross-strain comparison); **recover orphan anchors as synthetic singleton groups** (so no anchor gene is dropped from the inventory).
- **Advance rationale.** Anchor + ortholog inventory is built and QC'd. The 19 × 54 matrix is the inventory the step-1 question asks for; step 3 will decide *how* to present and interpret it.
