# Step 1 — Research question

## Context

User prompt: *"I want to add papers to the kg. The general plan is to use PaperBLAST: select interesting or high-rank genes, run PaperBLAST and identify candidate papers to add."*

This is a **KG-enrichment / curation** analysis rather than a hypothesis-driven one. It still runs through the 6-step flow, but the deliverable is a curated candidate-publication shortlist, not a scored result against a biological hypothesis.

Scope was settled through a clarifying dialogue. The forks resolved:

1. **End goal.** Among (gene→literature annotation layer / new omics datasets / curated reading list only), the researcher chose **new omics datasets** — find papers carrying expression/proteomics/metabolomics data to ingest as new Publication+Experiment+DE, the same shape as the existing 42 dataset-papers. PaperBLAST is the discovery route. Later refined to: **do both / just improve the KG / don't special-case by gene type** — so functional-characterization papers are also captured, not only omics datasets.
2. **Seed-gene basis.** **High-response genes per experiment**, breadth-first across all organism types/publications. Start with one round (~20–30 genes), expand to all experiments later. Heavy cross-experiment overlap expected → **caching + dedup are first-class requirements**.
3. **PaperBLAST access.** A probe found PaperBLAST is behind Cloudflare (HTTP 403, JS challenge). Resolved by a **hybrid workflow** (see Locked question). The researcher pushed back — correctly — on an earlier assumption that hypothetical/poorly-annotated seeds would yield few hits: PaperBLAST searches by *sequence similarity* against a large characterized-protein universe, so the query gene's own annotation state does not gate hits. The surviving concern is output-type (PaperBLAST returns characterization papers, not omics datasets), which is a triage step, not a seed problem.
4. **Deliverable size.** A shortlist of **~10–20+ candidate publications**.

## Locked question

> Build a deduplicated, prioritized shortlist of **~10–20+ candidate publications to add to the KG**, discovered by running PaperBLAST on a breadth-first set of high-response genes spanning all organism types, triaged for what each paper offers (omics dataset vs functional characterization) and filtered against the 42 publications already in the KG.

### Workflow (hybrid, three stages)

1. **Claude builds a prioritized gene table** (scripted, from the KG) — ~20–30 genes for round one, breadth-first across organism types (Prochlorococcus strains, Synechococcus, Alteromonas/heterotrophs), prioritized by being high-response across *many* experiments. Columns: `locus_tag`, `protein_id` (RefSeq, for reliable PaperBLAST lookup), `gene_name`, `product`, `organism`, `annotation_state`, priority signals, ready-to-paste query string. Sequence-null genes flagged/dropped.
2. **Researcher runs PaperBLAST offline** on that list and saves the result pages (the one non-scriptable link).
3. **Claude analyzes the downloads** (scripted) — parse hits → papers, dedup across genes, drop papers already in the KG, triage each survivor (which seed gene(s), top %identity, paper type, organism relevance) → ranked candidate-publication shortlist + builder-repo handoff note.

### In scope
- The prioritized gene list (stage 1).
- Analysis of whatever PaperBLAST results the researcher brings back (stage 3) → candidate-publication shortlist + handoff note.
- Caching and deduplication (cross-gene, and against the existing 42 KG publications).

### Out of scope
- The actual KG ingestion (lives in the builder repos `biocypher_kg` / `explorer`).
- The manual PaperBLAST step itself (researcher's, documented as a limitation per Rule 5).

### Deferred to step 2 (KG entries)
- Enumerate organisms and their experiments; confirm how broadly DE `rank_up`/`rank_by_effect` are populated vs needing a |log2FC|-based rank.
- Enumerate the existing 42 publications (the dedup target).
- Build the prioritized gene table (interactive — the researcher asked to work through this together).

### Deferred to step 3 (framing)
- Exact prioritization formula (cross-experiment frequency, per-organism quota, tie-breaks).
- Triage rubric for "candidate to add" (omics-dataset vs characterization; relevance threshold; %identity cutoffs).

## KG context

Grounding queries run during the dialogue. Counts and structural findings that shaped scope.

### How publications relate to the KG (load-bearing)
- `Publication` (42 nodes) carries metadata (title, authors, DOI, abstract, journal, year, study_type) and connects to the graph **only** via `Has_experiment` → `Experiment` → (`Changes_expression_of` → Gene, derived metrics, clustering, metabolite assays). **There is no direct Gene→Publication "characterized-in" edge.** A KG publication is a *dataset container*, not a literature annotation. (Source: `kg_schema`, `docs://guide/concepts`.)
- PaperBLAST returns the opposite: functional-characterization papers about sequence-similar homologs. Reconciling the two is the triage step in stage 3.

### Seed organism / example gene
- `resolve_gene('PMM0246')` → ntcA, "global nitrogen regulatory protein", Prochlorococcus MED4 (the example in the user's PaperBLAST URL).
- `gene_aa_sequence`: ntcA (PMM0246, 244 aa) and amt1 (PMM0263, 486 aa) have sequences; **PMM0838 has a null sequence** (`not_matched`).

### Experiment landscape (MED4 slice, orientation only)
- `list_experiments(organism='MED4', summary=True)`: 41 experiments matched; by omics MICROARRAY 19 / RNASEQ 11 / PROTEOMICS 8 / + vesicle & paired; by treatment light 10, carbon 8, nitrogen 8, viral/iron 3 each, coculture 2, etc. Whole KG has ~195 experiments across 37 organisms.
- DE rank fields: `experiment_Tolonen2006_ntcA_depletion_MED4` DE rows have **null** `rank_up`/`rank_by_effect` (one experiment checked; see `gaps_and_friction.md`).

### Access probe
- PaperBLAST web interface returns HTTP 403 / Cloudflare JS challenge to WebFetch and browser-UA curl. Drives the hybrid workflow. (See `gaps_and_friction.md`.)

## Decisions

**2026-05-30 — KG-enrichment analysis, hybrid PaperBLAST workflow.** Goal is to improve the KG by surfacing ingestible publications (omics datasets primarily, characterization papers also kept). PaperBLAST accessed manually by the researcher (Cloudflare block); Claude builds the gene list and analyzes downloads.

**2026-05-30 — Seed = high-response genes, breadth-first across organisms.** Round one ~20–30 genes; prioritize genes high-response across many experiments; per-organism quota for breadth; expand to all experiments later. Caching + dedup mandatory given overlap.

**2026-05-30 — Dedup against the existing 42 KG publications.** The shortlist must be genuinely new; PaperBLAST will re-surface papers already ingested.

## Surprises

- **Publications are dataset containers, not literature annotations** — no Gene→Publication edge. Shapes what "adding a PaperBLAST paper" means.
- **PaperBLAST is Cloudflare-gated** — no simple programmatic access.
- **Stored DE rank fields can be null** (Tolonen experiment) — "high-rank" may need a computed |log2FC| rank.
- **Some genes have null sequences** (PMM0838) — must filter seeds.
- **Annotation state does not gate PaperBLAST hits** (researcher correction) — PaperBLAST is sequence-similarity-based against a broad characterized-protein universe; the KG-internal ortholog-group size is not a predictor of external hits.

## Advance rationale

The deliverable (candidate-publication shortlist), the hybrid workflow, the seed strategy (breadth-first high-response, ~20–30 round one), and the dedup targets are all fixed. The remaining unknowns (exact prioritization formula, rank-field availability, triage rubric) are step-2/3 concerns that depend on data not yet pulled. Step 2 can begin — building the prioritized gene list, interactively with the researcher.

---

## Decide-gate checklist

- **Outputs produced.**
  - `.gitignore` (template)
  - `paper.md` (skeleton; Question populated)
  - `gaps_and_friction.md` (3 entries: Cloudflare block, null rank fields, null sequences)
  - `1_question/notebook.md` (this file)
  - No scripts/data/figures (step 1 is a conversation).
- **Results presented.** Grounding counts and findings shown inline above and in chat (publication model, MED4 experiment landscape, access probe, sequence/rank-field findings).
- **QC gate.**
  - Publication→Gene relationship checked against `kg_schema` → no direct edge; publications are dataset containers.
  - PaperBLAST access probed (WebFetch + curl) → Cloudflare 403 → hybrid workflow adopted.
  - DE rank fields inspected on one experiment → null → |log2FC| fallback flagged for step 2.
  - Seed sequence availability spot-checked → PMM0838 null → seed filtering required.
- **Decisions made this step.** Three (logged above): hybrid KG-enrichment workflow; breadth-first high-response seed (~20–30 round one); dedup against existing 42 publications.
- **Advance rationale.** Scope locked and approved by the researcher; remaining unknowns are data-dependent step-2/3 concerns. Step 2 begins by building the prioritized gene list interactively.
