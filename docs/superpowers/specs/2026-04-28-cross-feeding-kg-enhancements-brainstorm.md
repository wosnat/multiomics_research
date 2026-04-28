# Cross-feeding (Prochlorococcus / Alteromonas) — KG enhancement brainstorm

**Date:** 2026-04-28
**Status:** brainstorm only — no implementation, no scoping decisions
**Question:** If we want to identify cross-feeding between *Prochlorococcus* and *Alteromonas*, what is missing from the KG and could be added? Across: papers, schema, bioinformatic tools, organisms.

## Snapshot of current KG state (at time of brainstorm)

- 32 organisms; relevant pairs already present: Pro MED4 / NATL2A / MIT9312 / MIT9313 / SS120 paired in experiments with Alt MIT1002 / EZ55 / HOT1A3, plus an Alt (MarRef v6) reference proteome.
- 38 publications, 172 experiments; 77 are coculture, 12 publications are coculture-typed.
- Functional annotation today: eggNOG-derived KEGG, GO, COG, Pfam, EC, TIGR, Cyanorak, BRITE.
- DerivedMetric layer for column-level summaries (rhythmicity, diel amplitude, etc.).
- Paired RNA/protein only for Waldbauer 2012 MED4 diel; one PAIRED_RNASEQ_PROTEOME omics_type.

## What "cross-feeding" needs that the KG doesn't have

Cross-feeding is "a metabolite produced by one organism is consumed by another." None of those primitives exist as first-class nodes today:

- No **Metabolite / Compound** node anywhere
- No **Reaction** node (no KEGG R-numbers, no MetaCyc, no ModelSEED)
- No **Transporter** with substrate + directionality (import/export); transporter info is buried in protein-family strings
- No **secreted/excreted** flag on proteins (no signal-peptide property, no extracellular localization edge)
- No edge linking "gene A in MED4 produces X" with "gene B in HOT1A3 consumes X"

Everything has to be inferred indirectly through orthologs + functional descriptions, which is exactly why eggNOG-only annotation is the bottleneck.

## Papers worth adding

The corpus has the transcriptome backbone. The chemistry side is thin:

- **Metabolomics / exometabolome** for Pro–Alt or Pro–heterotroph pairs (Becker et al.; Braakman 2018 "Metabolic evolution and the self-organization of ecosystems"; Roth-Rosenberg DOM exchange; Szul et al. on Pro exudate).
- **Sher 2011** — "Response of *Prochlorococcus* ecotypes to co-cultured diverse marine bacteria" — foundational pairwise paper.
- **Christie-Oleza 2015** exudate composition (the 2017 + 2018 follow-ups are already in; this completes the set).
- **Stable-isotope probing / NanoSIMS** studies tracing carbon/nitrogen flow between Pro and Alt.
- **Capovilla / chitin** is in, but the rest of the Braakman/Cordero "metabolic coupling" series is not — those have explicit cross-feeding hypotheses to ground-truth against.
- Moreno-Cabezuelo 2023 is in as proteomics-only; its **metabolomics half is not**. Same applies to Kratzl 2024 and Ma 2022.

## Bioinformatic tools beyond eggNOG

eggNOG gives orthologs and category-level function. Cross-feeding needs metabolic context:

- **Genome-scale metabolic models** per strain (CarveMe, gapseq, ModelSEED) → produces Reactions and Compounds for free.
- **SMETANA / MICOM / COMETS** → community-level metabolic exchange predictions; literally outputs predicted cross-fed compounds.
- **TCDB / DeepTMHMM / SubCellPhonet** → real transporter substrate + topology, not just Pfam strings.
- **antiSMASH** → secondary metabolites and BGCs (siderophores, vitamins, infochemicals).
- **SignalP / SecretomeP / PSORTb** → secreted vs cytoplasmic, gates "can this even be exchanged."
- **DeepEC / CLEAN** → fills EC numbers eggNOG misses, especially where Alt coverage thins.

## Organisms worth adding or upgrading

The pairings are decent (5 Pro × 3 Alt). Gaps:

- **Alteromonas (treatment) entry has 0 genes** — placeholder. Replacing it with the actual EZ55/HOT1A3/MIT1002 mapping when papers used "Alteromonas sp." would unlock more pairings.
- **SAR11 / Pelagibacter** (HTCC1062, HTCC7211) — the other dominant heterotroph in oligotrophic gyres; co-occurs with Pro more than Alt does in the wild.
- **Roseobacter clade beyond Ruegeria** (e.g., *Phaeobacter*, *Sulfitobacter*) — known DMSP cross-feeders with cyanobacteria.
- **Marinobacter** — present in MarRef v6 as a reference proteome, but not as a profiled organism with experiments.
- **Synthetic-consortia** (3+ partner) papers — would let pairwise vs community effects be separated.

## Schema changes (rough order of payoff)

1. **Metabolite + Reaction nodes**, with `Gene -[catalyzes]-> Reaction -[produces|consumes]-> Metabolite`. Single biggest unlock — everything below leans on this.
2. **Transporter** node (or properties on Gene/Protein) with `substrate` + `direction` (import/export) + `localization` (inner/outer membrane, secreted).
3. **Cross-feeding hypothesis** edge: `Gene_A -[predicted_provides]-> Gene_B` annotated with the shared metabolite and evidence source (model prediction vs measured exometabolite vs isotope tracing).
4. **Exometabolome / spent-medium** as a first-class compartment alongside existing `vesicle` / `exoproteome` / `whole_cell` (the slot exists; data needs to flow in).
5. **Coculture phenotype** on `Tests_coculture_with` (growth-rate change, yield, dependence) — turns coculture from a context label into a measured outcome.

## Highest-leverage moves

- **#1 (metabolites + reactions)** — without it everything else is a workaround.
- **Lowest-effort move:** add the missing metabolomics data from papers already in the KG, since their other halves are already ingested.

---

*Captured from an interactive brainstorm. Not a design spec; not committed-to scope. Use as a starting point if/when cross-feeding moves to an active analysis.*
