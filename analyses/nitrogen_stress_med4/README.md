# Nitrogen stress vs death in Prochlorococcus MED4

## Question

Can we distinguish nitrogen-specific stress from general nutrient stress, general stress, and dying in MED4's transcriptional response?

## Current understanding

MED4's response to N-deprivation proceeds in 5 phases (based on Tolonen 2006 time course, cross-referenced against all MED4 stress experiments in KG):

1. **N-specific alarm** (0-3h): 6 genes — N-scavenging (cynA, urtA, glnA) + 3 hypotheticals. 4/6 never respond to any other stress.
2. **Metabolic restructuring** (6h): 25 genes — NtcA regulon activates, ATP synthase operon adjusts, ribosomal proteins shift. Heavily shared with carbon stress.
3. **General stress cascade** (12h): 38 genes — HLI proteins, photosynthesis changes, Calvin cycle. Light/viral/salt stress overlaps dominate.
4. **Damage control** (24h): 11 genes — chaperones (htpG, groEL). Protein quality control.
5. **Photosystem collapse / death** (48h): 28 genes — PSI/PSII shutdown, strongest HLI blast, DNA damage (lexA/recA).

Of 108 total N-responsive genes: **21 N-specific, 45 nutrient-general, 42 general-stress.**

**The Weissberg day 14 axenic data resembles phase 5 (dying), not phase 1-2 (N-stress).** The coculture survival state should be compared to phases 1-2.

N-deprivation rapidly becomes a carbon/energy crisis — 45/108 genes are shared with carbon stress. This reframes the paper: what looks like "nitrogen response" is largely "nutrient starvation leading to energy collapse."

## Key files

| File | Description |
|------|-------------|
| [methods.md](methods.md) | Data sources, gene identification, caveats |
| [gaps_and_friction.md](gaps_and_friction.md) | KG/MCP/skill issues found during analysis |
| [exploration/](exploration/) | Session logs with detailed findings |
| [data/](data/) | Extracted data files (CSVs) |
| [scripts/](scripts/) | Analysis and extraction scripts |
| [results/exploration/](results/exploration/) | Quick figures from exploration |
| [results/final/](results/final/) | Publication-quality figures |

## Exploration log

- [2026-03-30: Orientation and first responders](exploration/2026-03-30-orientation-and-first-responders.md) — KG orientation, Tolonen 3h first responders, cross-stress specificity, RNA vs protein divergence
- [2026-03-30: Timepoint classification](exploration/2026-03-30-timepoint-classification.md) — All 108 Tolonen genes classified as N-specific/nutrient/general by cross-stress comparison

## Data sources

8 nitrogen experiments + 22 other-stress experiments for MED4, from 11 publications across microarray, RNA-seq, and proteomics. See [methods.md](methods.md) for details.
