# Step 1 — Research question

## Context

User prompt (2026-04-27): *"I am working on Weissberg 2025, there are 2 main conditions: axenic and coculture. I want to know if Prochlorococcus is stressed in each of them."*

This step locks the research question. The publication was resolved via `list_publications(author="Weissberg")` to a single match: Weissberg, Aharonovich & Sher 2025, *bioRxiv*, "Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for Prochlorococcus Prolonged Survival" (DOI `10.1101/2025.11.24.690089`). The user is the first author. The KG indexes 10 experiments for this study, two organisms (*Prochlorococcus* MED4, *Alteromonas macleodii* HOT1A3), two omics (RNA-seq, proteomics), and three growth phases (`exponential`, `nutrient_limited`, `death`).

The user's published narrative is that axenic *Prochlorococcus* perishes under extreme N-deprivation while coculture cultures persist 90+ days, driven by *Alteromonas* recycling N back to *Pro*. The user's question reframes that: rather than asking *what mechanism saves coculture*, ask *what stress signatures does each condition exhibit on its own trajectory?*

## What I did

Used `superpowers:brainstorming` with the two methodology overrides applied (capture in this notebook, not `docs/superpowers/specs/`; advance to step 2, not `writing-plans`). Five clarifying questions, one at a time. Each question multiple-choice with a recommendation; user selected one option per question.

KG queries issued during this step:

```
list_publications(author="Weissberg", verbose=True)
```

→ resolved single matching publication (DOI `10.1101/2025.11.24.690089`), 10 experiments, organisms, omics types, growth phases.

## Results

### Locked research question

> Within each of axenic and coculture *Prochlorococcus* MED4 cultures from Weissberg 2025 (DOI `10.1101/2025.11.24.690089`), how do five stress axes — N-starvation, oxidative, proteotoxic, photosynthetic, and cell-death/late-stationary — evolve from an exponential-phase baseline through nutrient_limited and death-phase timepoints, scored independently in transcriptome and proteome?

Two parallel within-condition trajectories per axis per omics. No cross-condition statistical contrast (data limitation, see `../gaps_and_friction.md`).

### Scope

**In scope**

| Dimension | Value |
|---|---|
| Publication | Weissberg 2025 only (DOI `10.1101/2025.11.24.690089`) |
| Organism | *Prochlorococcus* MED4 only |
| Conditions | axenic and coculture, treated as two independent within-condition analyses |
| Omics | RNA-seq and proteomics, scored independently per axis |
| Stress axes (5) | N-starvation, oxidative, proteotoxic, photosynthetic, cell-death/late-stationary |
| TPs | exponential = baseline; nutrient_limited + death TPs are scored against baseline |

**Out of scope** — cross-condition contrast at matched TPs (data limitation); *Alteromonas* HOT1A3 data (off-question); external "healthy" reference from other studies (cross-study confounder); the "light" background factor (third experimental dimension in the KG, deferred); carbon/energy stress, envelope stress, isolated translation/ribosome stress (just-in-time formalization, deferred).

### Rejected alternatives (clarifying-dialogue trail)

| Question | Rejected | Why rejected |
|---|---|---|
| Q1 — meaning of "stressed" | (a) single general stress score | collapses axis-specific differences |
| Q1 — meaning of "stressed" | (b) N-starvation only | ignores oxidative-detox claim from abstract |
| Q1 — meaning of "stressed" | (c) oxidative only | same — too narrow |
| Q1 — meaning of "stressed" | (d) generalized "things going badly" | agnostic, masks mechanism |
| Q3 — comparison structure | (b) across-condition contrast at matched TPs | data does not support (user-reported) |
| Q3 — comparison structure | (c) both trajectory and across-condition contrast | (b) component infeasible |
| Q3 — comparison structure | (d) external healthy reference | cross-study confounder |
| Q4 — omics scope | (a) RNA-only / (b) protein-only | loses multi-omics framing |
| Q4 — omics scope | (d) combined RNA+protein single score | hides cross-omics divergence (informative) |
| Q5 — TPs | (a) all TPs including exponential | trivial baseline; effort wasted |
| Q5 — TPs | (b) score nutrient_limited + death without baseline | loses well-defined reference |

### Success criterion (end-of-analysis deliverable, for orientation)

For each (condition × axis × omics) cell — **2 × 5 × 2 = 20 trajectories** — a stress score curve from exponential baseline (= 0) through nutrient_limited and death TPs. Per trajectory, verdict on whether stress accumulates, when, and whether it plateaus or escalates. Per (condition × axis), cross-omics concordance check between RNA and protein. Cross-condition narrative is by visual inspection of curve shapes — explicitly **not** a statistical contrast.

How the stress score is computed (gene-set definition, scoring method, statistical treatment) is deliberately deferred to step 3 (framing) and step 4 (methods). Step 1 locks *what we are asking*, not *how we will answer it*.

## Surprises

None for step 1 (clarifying-dialogue step).

## Decisions

**2026-04-27 — Five-axis panel adopted.** User selected option (e) at Q1: stress measured as multi-axis panel rather than a single scalar. Five axes locked: N-starvation, oxidative, proteotoxic, photosynthetic, cell-death/late-stationary. Carbon/energy, cell-envelope, and isolated translation axes deferred (just-in-time formalization).

**2026-04-27 — Within-condition trajectory only (no cross-condition contrast).** User selected option (a) at Q3 with the explicit caveat *"would like b but we don't have the data for it."* The analysis answers "does stress accumulate within axenic? within coculture?" but not "is one more stressed than the other at matched TPs." Cross-condition reading is by visual inspection of curve shapes, not statistical. Logged in `../gaps_and_friction.md`.

**2026-04-27 — Both omics scored independently.** User selected option (c) at Q4. RNA and protein get separate stress trajectories per axis; cross-omics concordance becomes a finding, not a methodological requirement.

**2026-04-27 — Exponential as within-condition baseline.** User selected option (c) at Q5. Exponential phase is the reference (stress score ≈ 0 by construction); nutrient_limited and death TPs are scored against it.

## Decide-gate checklist

- **Outputs produced** — `1_question/notebook.md` (this file); `paper.md` Question section populated; `gaps_and_friction.md` entry F1 logged; scaffold (`paper.md` skeleton, `gaps_and_friction.md` header, `.gitignore`, `1_question/`).
- **Results presented** — locked question, scope (in/out), rejected-alternatives table, success criterion all shown inline above and in chat.
- **QC gate** — KG publication resolved via `list_publications(author="Weissberg")` → single match, DOI verified `10.1101/2025.11.24.690089`. Author / organism / omics / growth-phases all sourced from the KG response, not memory. No statistical computation in this step.
- **Decisions made this step** — four locked decisions above (panel, comparison structure, omics scope, baseline).
- **Advance rationale** — research question fully scoped, all five clarifying questions answered, design approved by researcher; ready for step 2 (KG entries) — identify the specific experiments, organisms, TPs, and omics rows that the question requires.
