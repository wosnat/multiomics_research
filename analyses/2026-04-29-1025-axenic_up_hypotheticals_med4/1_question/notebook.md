# Step 1 — Research question

## Context

User prompt (2026-04-29): *"I want to explore some of the upregulated genes in weissberg 2025 RNASEQ axenic. many of them have hypothetical function - I want to learn about them and their potential role in this experiment."*

This step locks the research question. The user came in with a per-gene contextualization goal scoped to one publication, one omics, one condition. The clarifying dialogue locked organism (MED4 only), goal shape (per-gene KG dossier as primary, pattern-finding as stretch), candidate-set definition (`significant_up` × `annotation_quality ≤ 1`), and a parallel **KG-enhancement sub-question** (log what the KG can't answer well; consolidate proposals at step 6).

## What I did

Used `superpowers:brainstorming` with the two methodology overrides applied (capture in this notebook, not `docs/superpowers/specs/`; advance to step 2, not `writing-plans`). Four clarifying questions, one at a time. Each question multiple-choice with a recommendation; the user added the KG-enhancement sub-question and the broader KG-facet list (cluster assignments, derived metrics, ortholog facets) as a mid-stream expansion.

KG queries issued during this step:

```
list_publications(author="Weissberg", verbose=True)
list_experiments(publication_doi=["10.1101/2025.11.24.690089"],
                 background_factors=["axenic"],
                 omics_type=["RNASEQ"],
                 verbose=True, limit=20)
```

→ established the Weissberg 2025 publication (1 match, DOI `10.1101/2025.11.24.690089`) and the 2 axenic RNA-seq experiments under it (MED4 single-contrast, HOT1A3 time-course). Locked on the MED4 axenic RNA-seq experiment as the sole data source.

## KG context

[KG] **Publication.** *Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for Prochlorococcus Prolonged Survival* — DOI `10.1101/2025.11.24.690089` (Weissberg, Aharonovich, Sher 2025); `list_publications(author="Weissberg")` returns exactly 1 result.

[KG] **In-scope experiment** (one):
- `experiment_id`: `10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic`
- Organism: *Prochlorococcus* MED4
- Background factors: `axenic`, `light`
- Omics: `RNASEQ`
- Treatment: `PRO99-lowN nutrient starvation` vs `PRO99-lowN exponential growth` (control)
- Statistical test: DESeq2
- Light / medium / temperature: continuous light, 20 µmol photons m⁻² s⁻¹, PRO99-lowN, 24 °C
- `is_time_course`: **false** (single contrast; one log2FC per gene)
- `table_scope`: `all_detected_genes`
- `gene_count`: 1849 detected
- `genes_by_status`: 405 `significant_up`, 472 `significant_down`, 972 `not_significant`
- Growth phase: `nutrient_limited`

[KG] **Out of scope at this analysis** (one): the parallel HOT1A3 axenic RNA-seq experiment (`…_hot1a3_rnaseq_axenic`, time-course d18 / d31 / d60+89, 778 sig_up across the experiment). Held as a possible follow-up analysis if MED4 dossier-building proves informative.

[KG] **`annotation_quality` semantics** (from `genes_by_function` schema): 0 = hypothetical, 1 = has description (no named product), 2 = named product, 3 = well-annotated. The candidate set uses `annotation_quality ≤ 1` (strictly poorly characterized). Set size — i.e., how many genes are simultaneously `significant_up` and `annotation_quality ≤ 1` — is unknown until step 2 quantifies it.

## Results

### Locked research question

> In the *Prochlorococcus* MED4 axenic RNA-seq experiment of Weissberg et al. 2025 (PRO99-lowN N-starvation vs PRO99-lowN exponential growth, single contrast, `nutrient_limited` phase), among genes called `significant_up` (DESeq2 padj < 0.05, log2FC > 0) that are annotated as hypothetical or poorly characterized (`annotation_quality ≤ 1`), what does the KG already tell us — per gene — that bears on each gene's potential role under N-starvation?

### Per-gene KG facets to assemble (primary deliverable)

For each candidate gene **and for each of its orthologs** (via homolog group):

| Facet | KG source (provisional — confirmed at step 2) |
|---|---|
| Annotation | `gene_overview` / `gene_details` (product, description, gene_category, BRITE) |
| Ontology | `gene_ontology_terms` (GO, KEGG, COG, etc.) |
| Genomic neighborhood / operon | `gene_details` (start/end/strand, neighbors) |
| Cluster assignments | `gene_clusters_by_gene` (any clustering analyses participated in) |
| Derived metrics | `gene_derived_metrics` (every `metric_type` × `value_kind` available) |
| Expression in other studies | `gene_response_profile` / `differential_expression_by_gene` (other organisms / treatments / TPs) |
| Homolog group + orthologs | `gene_homologs` then iterate facets above per ortholog |

Tool selection above is provisional — confirmed against actual API signatures at step 2 (per `feedback_check_api_before_cypher`).

### Stretch question (b)

> If patterns emerge across the dossier — shared ortholog groups, neighborhood proximity to known N-regulon members, co-occurrence in clusters, shared derived-metric profiles, or co-expression with known N-stress genes in other studies — document candidate functional bins for the upregulated hypotheticals.

Activated only if step 5 dossier review surfaces visible patterns; not pre-architected.

### KG-enhancement sub-question (meta)

> As the analysis proceeds, log every place where the KG fails to answer a question we want to ask about these genes, and what addition to the KG would help future analyses of this kind.

Mechanics: each blocked question gets a dated entry in `gaps_and_friction.md` (already required by methodology). Step 6 evaluation includes a dedicated **KG enhancement proposals** sub-section that consolidates the entries into a prioritized list (frequency × impact). The analysis itself does not run new bioinformatics — it harvests gap signal.

Researcher-supplied leads already on the table at step 1 (will become F-entries the moment the analysis tries the corresponding question and is blocked):
- Amino-acid sequences absent from the KG due to per-node size limits — needed for ad-hoc BLAST / family searches / structure prediction.
- Batch-bioinformatics layers: Pfam / InterPro domain hits, SignalP signal-peptide predictions, TMHMM transmembrane predictions, AlphaFold structure summaries.

### Out of scope (this analysis)

- HOT1A3 axenic RNA-seq.
- All coculture data (RNA-seq or proteomics).
- The proteome dimension on the MED4 axenic side.
- Re-doing differential-expression statistics — the publication's `significant_up` call is taken as given.
- Wet-lab follow-up on candidate hypotheticals.
- Running new bioinformatics ourselves (BLAST, domain searches, structure prediction). The KG-enhancement sub-question flags what would be useful; it does not produce it.

### Rejected alternatives (clarifying-dialogue trail)

| Question | Rejected option | Why rejected |
|---|---|---|
| Q1 — analysis identity | (b) branch off active discordance analysis | RNA-seq-only / no proteome dimension; doesn't need paired-FC machinery |
| Q2 — organism scope | (b) HOT1A3 only | first pass on the streamlined MED4 genome makes hypotheticals more interpretable; HOT1A3 follow-up possible |
| Q2 — organism scope | (c) both organisms | doubles the work for a first dossier; pattern-finding (b stretch) easier on one organism first |
| Q3 — end goal | (c) prioritized shortlist | scoring on shaky ground without a per-gene context dossier first |
| Q4.1 — "upregulated" | (ii) custom magnitude cut | re-doing DE statistics is out of scope; magnitude noted in dossier separately |
| Q4.1 — "upregulated" | (iii) all positive log2FC | publication's significance call is the anchor; sub-threshold genes are noise for "potential role" reasoning |
| Q4.2 — "hypothetical" | (I) `annotation_quality = 0` strict | DUF-style "has description but no named product" genes share the same "poorly understood, potentially interesting" character |

### Success criterion (end-of-analysis deliverable, for orientation)

A per-gene dossier table indexed by MED4 locus tag, with columns for each KG facet listed above, plus a parallel ortholog-facet section. Ranked or grouped by whatever organization makes the patterns most readable (decided at step 5). Plus, attached to the same analysis, a prioritized list of KG-enhancement proposals harvested from `gaps_and_friction.md`.

How the dossier is laid out, how orthologs are aggregated when there are many, and how patterns (if any) are surfaced are deliberately deferred to step 3 (framing) and step 4 (methods).

## Surprises

None for step 1 (clarifying-dialogue step).

## Decisions

**2026-04-29 — New analysis, not branch from active discordance analysis (option (a) at Q1).** The proteome-transcriptome discordance analysis (`analyses/2026-04-27-1638-…`) stays paused at end-of-step-2; this analysis runs independently. Same publication, different question: RNA-seq-only per-gene contextualization vs paired-FC discordance.

**2026-04-29 — MED4 axenic RNA-seq only (option (a) at Q2).** Single experiment, single contrast. HOT1A3 axenic RNA-seq is held as a possible follow-up analysis once the MED4 dossier-building approach is validated.

**2026-04-29 — Per-gene KG dossier as primary; pattern-finding as stretch (option (a) with (b) as stretch at Q3).** Primary deliverable is the per-gene assembly of all KG facets for each candidate gene and its orthologs. Stretch is activated only if patterns are visible after the dossier is built.

**2026-04-29 — Candidate set: `significant_up` × `annotation_quality ≤ 1` (option (i)+(II) at Q4).** Anchors to the publication's DESeq2 significance call (no re-statistics) and to a broad "poorly characterized" annotation cut. Set size to be quantified at step 2; if it falls outside ~10–100 the parameters get re-examined.

**2026-04-29 — Per-gene facets expanded to include cluster assignments and derived metrics (researcher-directed expansion).** After I proposed annotation / ontology / neighborhood / homologs / cross-study expression as the dossier facets, the user added cluster assignments and derived metrics, applied to both the gene itself and to each of its orthologs.

**2026-04-29 — KG-enhancement sub-question runs alongside the analysis (researcher-directed addition).** The meta sub-question — "what could be added to the KG to help with this kind of analysis?" — is answered by accumulating dated F-entries in `gaps_and_friction.md` as the analysis hits limitations, then consolidating them into a prioritized list at step 6. This analysis does not itself run new bioinformatics; it harvests gap signal. Researcher-supplied initial leads (AA sequences, batch bioinformatics layers) are documented in the sub-question section above and will become F-entries when actually hit.

## Decide-gate checklist

- **Outputs produced** — `1_question/notebook.md` (this file); `paper.md` Question section populated; `gaps_and_friction.md` header with the active KG-enhancement sub-question framing; scaffold (`paper.md` skeleton, `gaps_and_friction.md` header, `.gitignore`, `1_question/`).
- **Results presented** — locked question, per-gene KG facets table, stretch question, KG-enhancement sub-question framing, out-of-scope list, rejected-alternatives table, success criterion all shown inline above and in chat.
- **QC gate** — KG publication identity verified via `list_publications(author="Weissberg")` → 1 result, DOI `10.1101/2025.11.24.690089`, year 2025, study type "Transcriptomics and Proteomics (multi-omics)" — sourced from KG response, not memory. Experiment scope verified via `list_experiments(publication_doi=[DOI], background_factors=["axenic"], omics_type=["RNASEQ"], verbose=True)` → 2 results; selected `…_med4_rnaseq_axenic` (single-contrast, 1849 detected, 405 sig_up). No statistical computation in this step.
- **Decisions made this step** — five locked decisions above (analysis identity, organism scope, goal shape, candidate-set definition, sub-question structure).
- **Advance rationale** — research question fully scoped, all four clarifying questions answered, design approved by researcher; ready for step 2 (KG entries) — formal selection of the candidate-gene set (`significant_up` × `annotation_quality ≤ 1` from the locked experiment), with QC of the set size, the annotation_quality distribution, and a smoke-test of which per-gene KG facets actually return data on a few candidates.
