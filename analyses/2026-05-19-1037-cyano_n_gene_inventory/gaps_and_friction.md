# Gaps and friction log

Append-only. Captures methodology, KG data, and tooling friction encountered during this analysis. Distinct from decisions (which live in step `notebook.md` files).

## 2026-05-19 — Step 1 retrofitted to disk after dialogue

**What happened.** The step-1 brainstorming dialogue (locked question, gene-set candidate comparison, scope decisions) ran in chat before the scaffold was created. The methodology specifies that the empty scaffold (`paper.md` skeleton, `gaps_and_friction.md`, `.gitignore`, empty `1_question/notebook.md`) should be written first, then the dialogue progressively populates `1_question/notebook.md` live. Here, the scaffold was created after the dialogue concluded and approved.

**Workaround.** `1_question/notebook.md` was reconstructed from the chat transcript at decide-phase close. Capture is faithful (queries, counts, alternatives, decisions all present), but the live-typed-as-it-happens pattern was not followed.

**Impact on methodology.** Low. The dialogue content is fully captured; what was lost is the audit trail of when each clarifying exchange happened relative to each KG query. If this becomes a pattern across multiple analyses, consider either: (a) a pre-dialogue scaffold hook that fires when research-methodology is loaded, or (b) acceptance that step 1's notebook is allowed to be retrofitted at decide time, with the scaffold-before-dialogue rule relaxed.

## 2026-05-19 — KEGG `ko00910` "Nitrogen metabolism" is far narrower than its name suggests in Prochlorococcus

**What happened.** During step-1 grounding, KEGG pathway `ko00910` returned only 6 genes in MED4 — cynABDS (cyanate operon), glnA, glsF. This excludes NtcA, GlnB/P-II, urease, urea transport, AMT1 — all canonical N-machinery in Prochlorococcus per Cyanorak. KEGG's pathway boundaries are not the same as a cyanobacterial researcher's intuition for "nitrogen metabolism gene".

**Workaround.** Switched to Cyanorak roles `E.4 ∪ D.1.3` (28 genes in MED4 with the expected content). Cyanorak is curated specifically for marine picocyanobacteria and aligns better with the field's working definition.

**Impact on methodology.** Methodology gap: when scoping a gene set, KEGG pathway names can be misleading proxies. Documenting "compare three candidate definitions on a sample organism before locking" as a step-1 grounding pattern would prevent silent narrowing.

## 2026-05-19 — Anti-hallucination correction: Synechococcus habitat labels came from intrinsic knowledge, not KG

**What happened.** During step-1 grounding I labeled Synechococcus strains "marine", "coastal", "freshwater" and Thermosynechococcus as "thermophile" — but the KG `clade` field is `null` for all of them. The labels came from my intrinsic knowledge of these well-studied lab strains, not from a KG query. This violates Rule 1 (KG is sole data source) and is exactly the kind of slip the methodology aims to catch.

**Workaround.** Researcher flagged the inconsistency and pointed to the Cyanorak source CSVs (see KG-bug entry below), which carry the authoritative clade and subcluster info. Labels are being re-derived from those CSVs in step 2.

**Impact on methodology.** Reinforces an existing rule rather than revealing a new gap. Worth a process check: before claiming any organism property in step 1, verify it came from a KG field, not background knowledge.

## 2026-05-19 — KG bug: Synechococcus `clade` field is null; source data populates it

**What happened.** `list_organisms` returns `clade=null` for every Synechococcus and Thermosynechococcus strain in the KG. The Cyanorak source CSV at `/home/oweissberg/work/multiomics/multiomics_biocypher_kg/data/Cyanorak  Organism Table  synechococcus.csv` has the full Scanlan/Mazard SubCluster + Clade + SubClade + Pigment classification for the 4 marine Synechococcus in our scope:

| Strain | SubCluster | Clade | SubClade | Pigment |
|---|---|---|---|---|
| WH8102 | 5.1 | III | IIIa | 3c |
| WH7803 | 5.1 | V | V | 3a |
| CC9311 | 5.1 | I | Ia | 3dA |
| BL107 | 5.1 | IV | IVa | 3dA |

PCC 7002, PCC 7942, UTEX 2973, and Thermosynechococcus vestitus BP-1 are outside the Cyanorak picocyanobacteria scope and would need a different taxonomy source if their clade is needed downstream.

**Workaround.** For this analysis, use the Cyanorak source CSV as the authoritative clade assignment for the 4 marine Synechococcus; leave the 4 non-marine strains without a sub-cluster clade label and refer to them by genus+species+strain.

**Impact on KG.** Loader bug — the `Clade` (and likely `SubCluster`, `SubClade`, `Pigment type`) columns from the Cyanorak Synechococcus table are not propagating into `Organism.clade` in the graph for Synechococcus nodes. Same loader probably handles both Pro and Syn tables; needs to verify why one populates and the other doesn't. May also be that the Synechococcus `clade` value is being mapped to a different KG property (e.g., subcluster) and the API isn't exposing it. Worth investigating upstream.

## 2026-05-19 — Methodology compression: ortholog expansion + QC pulled into step 2 by researcher redirect

**What happened.** Step 2 was about to close after the initial 4 scripts (anchor extraction, strain table, N publications, outlier QC). The researcher requested pulling the ortholog-bridge expansion (script 04 + 05) and the cross-strain inventory QC (script qc_ortholog_inventory) into step 2 before commit, so the full gene-set inventory exists before step 3 framing. This is methodologically step-4/5 work (methods + analyze) pulled into step-2 (KG entries).

**Workaround.** Did it. Inventory matrix (19 strains × 54 groups) sits at `2_kg_selection/data/05_inventory_matrix.csv` at end of step 2. Step 3 (framing) becomes a presentation-and-ordering decision rather than a methodology decision.

**Impact on methodology.** Consistent with just-in-time formalization — *"look at the data before drafting the plan"*. The trade-off: step 2's commit is much larger than typical, and step 3 will be unusually short. Worth a methodology observation: for inventory-style analyses (no hypothesis test), the 6-step flow may naturally compress into 2-3 substantive steps. If this pattern recurs, consider an "inventory analysis" methodology variant. One occurrence is a note; revisit if seen again.

## 2026-05-19 — Orphan-anchor inventory undercount (eggnog Cyanobacteria-level lineage asymmetry)

**What happened.** Of the 375 Cyanorak `E.4 ∪ D.1.3` anchor locus_tags, 9 returned no eggnog Cyanobacteria-level ortholog group. A post-hoc audit (`qc_ortholog_inventory.py` Q7b) showed that 6 of these 9 caused real undercounts in the 19-strain × 54-group inventory matrix: the strain shows cell=0 for a gene it actually has, because the orphan's lineage-specific ortholog group is not represented at the eggnog Cyanobacteria tier.

Affected cells (5 strains × 6 genes):

| Strain | Gene | Orphan locus_tag(s) |
|---|---|---|
| Prochlorococcus MED4 | cynA | PMM0370 |
| Prochlorococcus MIT9313 | ureE | PMT2232 |
| Synechococcus WH8102 | nrtP | SYNW2462, SYNW2463 |
| Synechococcus CC9311 | cynH | sync_2840, sync_2903 |
| Synechococcus CC9311 | amt2 | sync_2280 (gene fully absent from inventory) |
| Synechococcus sp. BL107 | focA | BL107_06829 |

**Workaround.** Accept and document. Considered (but rejected): augment inventory with cyanorak-curated groups as a fallback for orphans (would mix tiers); per-cell patch (would mix data sources). For step 6, any interpretive claim that hinges on one of these 6 cells must be cross-checked against the original CSV.

**Impact on KG.** Suggests that the eggnog-Cyanobacteria-level orthology in the KG is not symmetric across cyano lineages — for some genes, only one of the lineage-specific paralog groups makes it to the Cyanobacteria tier (the other stays at a more specific tier or has no group). Worth a longer-term KG audit: which Cyanorak roles have eggnog Cyanobacteria-level coverage in both Prochlorococcus and Synechococcus lineages, and which are asymmetric? For this analysis, the inventory matrix is correct in 1020 of 1026 cells (99.4%) and incorrect in 6.

## 2026-05-19 — KG bug: NATL1A and NATL2A clade is wrong in KG (LLII instead of LLI)

**What happened.** `list_organisms` returns `clade=LLII` for both Prochlorococcus NATL1A and NATL2A. The Cyanorak source CSV (`Cyanorak  Organism Table  prochlorococcus.csv`) lists both as **LLI** — which matches the established literature (NATL1A and NATL2A are the canonical LLI reference strains). All 8 other Prochlorococcus strains in our scope match between the KG and the source CSV; only these two are wrong.

**Workaround.** Use the Cyanorak source CSV as the authoritative Prochlorococcus clade assignment. Corrected clade table for this analysis:

| Strain | KG `clade` | Corrected (Cyanorak CSV) |
|---|---|---|
| NATL1A | LLII | **LLI** |
| NATL2A | LLII | **LLI** |
| (all others) | (correct) | (matches) |

**Impact on KG.** Likely a parsing or mapping error in the organism loader — the LLI/LLII columns may have been swapped or misread for these two strains specifically. Worth a targeted fix upstream. Impact on field-facing tools that filter or group by clade is non-trivial: any user querying for "LLI Prochlorococcus" via the KG would currently miss the two best-characterized LLI strains.
