# Step 1 — Research question

## Context

User prompt (2026-04-27): *"I want to know about the proteome/transcriptome correlation/differences. is it related to specific pathways? specific protein sizes?"*

This step locks the research question. The user came in with a generic discordance question and two pre-specified candidate axes (pathway, size). Through five clarifying questions the dataset, contrast, and additional classification axes were progressively scoped down to a single publication and up to a broad multi-axis decomposition of the discordance signal.

Critical constraint clarified during dialogue: the KG holds **no abundance / read-count data** — only differential-expression fold-changes and a small set of derived metrics. This rules out any analysis that needs absolute mRNA or protein levels (e.g., translation-efficiency / protein-per-transcript ratios). The locked question is purely a **fold-change-vs-fold-change discordance** analysis at matched timepoints.

## What I did

Used `superpowers:brainstorming` with the two methodology overrides applied (capture in this notebook, not `docs/superpowers/specs/`; advance to step 2, not `writing-plans`). Five clarifying questions, one at a time. Each question multiple-choice with a recommendation; the user selected one option per question and added one mid-stream correction (Q5 follow-up — see Decisions).

KG queries issued during this step:

```
list_experiments(summary=True)
list_derived_metrics(summary=True)
list_publications(limit=50)
list_experiments(publication_doi=["10.1101/2025.11.24.690089"], limit=20, verbose=True)
```

→ established the 4 within-publication paired transcriptome+proteome studies in the KG (Waldbauer 2012, Weissberg 2025, Kratzl 2024, Ma 2022) and the 10 experiments + matched-timepoint structure of Weissberg 2025.

## Results

### Locked research question

> In *Prochlorococcus* MED4 and *Alteromonas macleodii* HOT1A3 under PRO99-lowN nitrogen starvation (axenic and in coculture, Weissberg et al. 2025, DOI `10.1101/2025.11.24.690089`), where mRNA and protein responses disagree per gene at matched timepoints, what gene-level features explain the disagreement?

The discordance question is decomposed across multiple gene-level classification axes (see below). Step 3 framing prioritizes among these axes. Step 1 locks *what we are asking*, not *how we will answer it*.

### Scope — paired data in scope

| Organism | Condition | Paired RNA-seq + proteomics timepoints |
|---|---|---|
| *Prochlorococcus* MED4 | N-starvation, coculture | d18, 31, 60, 89 |
| *Prochlorococcus* MED4 | N-starvation, axenic | d14 (RNA-seq matches first proteome point; later proteome timepoints excluded) |
| *Alteromonas* HOT1A3 | N-starvation, coculture | d18, 31, 60, 89 |
| *Alteromonas* HOT1A3 | N-starvation, axenic | d18, 31 |

Per-organism analysis (separate locus-tag spaces); cross-organism comparison runs at the pathway / ortholog-group level.

### Out of scope

- The two single-timepoint coculture-vs-axenic RNA-seq contrasts in Weissberg 2025 (no matching proteomics).
- MED4 axenic d31, d89 proteome-only data — excluded because there is no matching RNA-seq, and "no extractable RNA" cannot be treated as "mRNA biologically absent" (logged in `../gaps_and_friction.md` F2).
- All other publications. Held in reserve for step-6 cross-condition replication: Waldbauer 2012 (DOI `10.1371/journal.pone.0043432`, MED4 diel, PAIRED_RNASEQ_PROTEOME with `peak_time_*`, `protein_transcript_lag_h`, `damping_ratio` derived metrics already in the KG).

### Classification axes for discordant genes

Listed in priority order as proposed at step 1; step-3 framing will narrow to the strongest axes for the actual hypothesis, others become exploratory secondary questions or get harvested as caveats at step 6.

1. **Pathway / function** — KEGG, KEGG BRITE pathway/metabolism trees, COG, Cyanorak (Pro only), GO Biological Process, Pfam.
2. **Protein size** — AA length and/or molecular_mass (`Polypeptide.sequence_length`, `Polypeptide.molecular_mass`).
3. **Protein architecture / localization** — primarily GO Cellular Component (the populated localization source in the KG, accessed via `Gene_located_in_cellular_component`) and KEGG BRITE location-implying trees (Transporters br02000, Secretion system br02044, Bacterial motility br02035, Photosynthesis proteins br00194, Ribosome br03011, Cell envelope), accessed via `Gene_has_kegg_ko` → `Kegg_term_in_brite_category`. Supplementary where annotation exists: `Polypeptide.signal_peptide`, `Polypeptide.transmembrane_regions` (sparse coverage expected).
4. **Direction-of-discordance category** — partition genes into mRNA-up/protein-flat, mRNA-flat/protein-up, opposite-direction (and concordant as the comparator) before applying axes 1–3, 5–11. Different categories likely have different biological causes.
5. **Ortholog conservation breadth** — `Gene.closest_ortholog_group_size`, `Gene.closest_ortholog_genera` count.
6. **Cross-condition consistency** — same gene, axenic vs coculture: discordance present in both = property of the gene; discordance only in coculture = condition-specific regulation.
7. **Cross-organism conservation of discordance** — ortholog-pair join MED4 ↔ HOT1A3 via `Gene_in_ortholog_group`; do orthologs show the same discordance signature?
8. **Temporal pattern** — early-only (d18), persistent, late-only (d89), and trajectory shape across the timecourse. Early discordance suggests translation-rate / lag differences; late suggests turnover / stability differences.
9. **Operon / genomic context** — derivable from `Gene.start` / `Gene.end` / `Gene.strand`; per-organism operon prediction at step 4. Genes in the same operon should have correlated mRNA but can have decorrelated protein.
10. **Hydrophobicity / GRAVY** — computed from `Polypeptide.sequence_length`-derived features (sequence not directly stored on the node; needs extraction at step 4 if MS-detection bias dominates). Not in KG.
11. **Annotation quality** — `Polypeptide.is_reviewed` (Swiss-Prot vs TrEMBL), `Polypeptide.annotation_score`, `Gene.function_description == "hypothetical"` flag. Uncharacterized proteins may show artefactual discordance from poor annotation.

### Deferred to step 4 methods (per just-in-time formalization)

- Pathway annotation source priority (likely COG + KEGG as cross-organism backbone, supplemented by Cyanorak for MED4 only).
- Discordance metric (signed Δ, residual-from-regression, quadrant classification, etc.) — driven by what the data look like.
- Significance gating (all detected vs significant in either omic vs union).
- Per-timepoint vs trajectory-summary handling.

### Rejected alternatives (clarifying-dialogue trail)

| Question | Rejected | Why rejected |
|---|---|---|
| Q1 — kind of correlation/difference | (b) steady-state level discordance | needs absolute abundance, not in KG (only fold-changes + derived metrics) |
| Q1 — kind of correlation/difference | (c) temporal lag from raw trajectories | same reason — would need raw abundance over time |
| Q3 — scope across publications | (b) all 4 within-publication paired studies | adds non-marine model organisms (E. coli, P. putida, Synechococcus elongatus) — different biology, more confounds |
| Q3 — scope across publications | (c) all 4 + across-publication matched comparisons | mixes platforms, strains, study designs — unmanageable confounds |
| Q4 — contrast within Weissberg 2025 | (a) all four N-starvation paired contrasts × all timepoints | fine but no need to cap there — locked scope grew to (a) + Waldbauer reserve |
| Q4 — contrast within Weissberg 2025 | (c) deepest single contrast (MED4 coculture only) | wastes the symmetry between MED4 and HOT1A3 |
| Q5 — how to handle MED4 axenic d31, d89 proteome-only | (b) use as "extreme discordance" sub-analysis | researcher correction: "no extractable RNA" ≠ "mRNA absent"; cannot be framed as biological persistence |
| Q5 — how to handle MED4 axenic d31, d89 proteome-only | (c) persistent-protein-first reframing | same correction applies — measurement failure ≠ biological absence |

### Success criterion (end-of-analysis deliverable, for orientation)

For each (organism × condition × timepoint) paired contrast — 11 paired observations across the 4 organism-condition rows — produce a per-gene discordance score with category labels and pathway / size / architecture / conservation / annotation-quality / etc. annotations. Across the 11 contrasts, identify which classification axes carry the strongest signal: do certain pathways, sizes, localizations, or conservation classes show systematically large discordance? Are those signals consistent across timepoints, conditions, and organisms?

Step 5 produces scored tables, scatter plots (Δprotein vs ΔmRNA per timepoint), pathway / size / architecture stratifications, and cross-condition / cross-organism consistency checks. Step 6 harvests caveats — including methodological confounds (MS detection bias, low-abundance protein under-sampling) — and finalizes the paper.

How the discordance score is computed, how categories are defined, and how each axis is operationalized are deliberately deferred to step 3 (framing) and step 4 (methods).

## Surprises

None for step 1 (clarifying-dialogue step).

The closest thing to a surprise was the discovery that the KG already carries 6 paired-omics derived metrics (`peak_time_transcript_h`, `peak_time_protein_h`, `protein_transcript_lag_h`, `diel_amplitude_transcript_log2`, `diel_amplitude_protein_log2`, `damping_ratio`) — but all from a single diel paired experiment (Waldbauer 2012), not from Weissberg 2025. These are reserved for step-6 cross-condition replication.

## Decisions

**2026-04-27 — Scope locked to Weissberg 2025 (option (a) at Q3).** Within-publication paired RNA-seq + proteomics analysis. The 3 other publications with within-publication paired data (Waldbauer 2012, Kratzl 2024, Ma 2022) are out of scope; Waldbauer is held in reserve for step-6 replication, the two non-marine studies (Kratzl, Ma) are dropped because their organisms (E. coli, P. putida, Synechococcus elongatus) are too biologically distant.

**2026-04-27 — Both organisms in scope, all paired timepoints (option (b) at Q4 expanded).** MED4 (5 paired timepoints across coculture d18/31/60/89 + axenic d14) and HOT1A3 (6 paired timepoints across coculture d18/31/60/89 + axenic d18/31). Per-organism analysis with cross-organism comparison at the pathway / ortholog-group level.

**2026-04-27 — Matched-FC discordance only (option (a) at Q5, with researcher mid-stream correction).** Discordance is computed as per-gene paired (ΔmRNA, Δprotein) at matched timepoints. The MED4 axenic d31, d89 proteome-only data is excluded from the analysis because there is no matching RNA-seq, and "no extractable RNA" ≠ "mRNA biologically absent" (researcher-supplied correction; logged in `../gaps_and_friction.md` F2 and as a user-level feedback memory).

**2026-04-27 — Eleven classification axes adopted (researcher-directed expansion).** After I proposed three priority axes (architecture-bundle, direction-of-discordance, cross-condition consistency) the user requested all axes I had listed be retained. The locked axes list (1–11 above) is therefore deliberately broad; just-in-time formalization shifts to step 3 (framing) where priority among axes will be set.

**2026-04-27 — KEGG BRITE added to axis 3 (researcher-prompted).** BRITE location-implying trees (Transporters, Secretion system, Bacterial motility, Photosynthesis proteins, Ribosome, Cell envelope) are explicitly listed under axis 3 alongside GO Cellular Component. BRITE pathway/metabolism trees additionally feed axis 1.

## Decide-gate checklist

- **Outputs produced** — `1_question/notebook.md` (this file); `paper.md` Question section populated; `gaps_and_friction.md` entries F1 (KG metadata gap — MED4 axenic RNA-seq has no timepoint label) and F2 (anti-hallucination correction — measurement failure ≠ biological absence) logged; scaffold (`paper.md` skeleton, `gaps_and_friction.md` header, `.gitignore`, `1_question/`).
- **Results presented** — locked question, paired-data scope table, out-of-scope list, eleven classification axes, rejected-alternatives table, and success criterion all shown inline above and in chat.
- **QC gate** — KG publication identity verified via `list_publications(limit=50)` → DOI `10.1101/2025.11.24.690089` matched; experiment structure verified via `list_experiments(publication_doi=[DOI], verbose=True)` → 10 experiments, organisms, timepoints, table_scope all sourced from the KG response, not memory. Within-publication paired-omics filter via the omics_types list on each publication (see `list_publications`) → 4 publications matched. No statistical computation in this step.
- **Decisions made this step** — five locked decisions above (publication scope, organism scope, discordance flavor, axes list, BRITE inclusion).
- **Advance rationale** — research question fully scoped, all five clarifying questions answered, design approved by researcher; ready for step 2 (KG entries) — formal selection of the 8 paired experiments × matched timepoints from the KG, with QC of detection-coverage and gene-counts per pair.
