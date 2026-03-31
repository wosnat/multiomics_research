# 2026-03-30: Orientation and N-stress first responders

## Session goal

Explore nitrogen stress in MED4 across the KG. Distinguish N-specific from general stress and dying.

## Orientation

### Nitrogen experiments for MED4

8 experiments from 3 publications:

| Publication | Year | Platform | Design |
|---|---|---|---|
| Tolonen et al. (10.1038/msb4100087) | 2006 | Microarray | N-deprivation time course (0, 3, 6, 12, 24, 48h), urea as sole N, cyanate as sole N |
| Read et al. (10.1038/ismej.2017.88) | 2017 | RNA-seq | N-depleted vs N-replete (3h, 12h, 24h). Filtered: top 50% by expression. |
| Weissberg et al. (10.1101/2025.11.24.690089) | 2025 | RNA-seq + Proteomics | Long-term N-starvation (14-89 days), axenic vs coculture with HOT1A3 |

### All MED4 experiments (for cross-stress comparison)

30 experiments across 8 treatment types:
- nitrogen_stress (8), light_stress (8), carbon_stress (4), coculture (3), iron_stress (3), viral (3), phosphorus_stress (1), salt_stress (1)

## Top DE genes under N-starvation (Weissberg axenic RNA-seq, day 14)

### Upregulated (602 sig entries, 466 genes)

Top categories: Unknown (163), Stress response (59), Coenzyme metabolism (41), Carbohydrate metabolism (28)

Top genes by |log2FC|:
- PMM1404 (hli) — high light inducible protein, log2FC 9.6
- PMM0689 (hli) — log2FC 7.7
- PMM1828 (hypothetical) — log2FC 6.7
- PMM0374 (tatA) — twin-arginine translocase, log2FC 6.4
- PMM1385 (hli) — log2FC 6.0
- PMM0370 (cynA) — cyanate transporter, log2FC 5.2
- PMM0336 (PTOX) — plastoquinol terminal oxidase, log2FC 5.0

**Note:** HLI proteins dominate the top upregulated genes. These are general stress markers, not N-specific.

### Downregulated (640 sig entries, 499 genes)

Top categories: Unknown (126), Translation (58), Photosynthesis (50), Stress response (49)

Top genes by |log2FC|:
- PMM1992 (hypothetical) — log2FC -7.3
- PMM0227 (sat) — sulfate adenylyltransferase, log2FC -6.4
- PMM0468 (psaJ) — PSI subunit IX, log2FC -6.0
- PMM0627 (pcb) — chlorophyll-binding protein, log2FC -5.8
- PMM0469 (psaF) — PSI reaction centre, log2FC -5.5
- PMM0785 (prkB) — phosphoribulokinase, log2FC -5.5
- PMM0552 (csoS2) — carboxysome shell, log2FC -5.4
- PMM1439 (atpE) — ATP synthase epsilon, log2FC -5.3

**Note:** Massive photosynthesis shutdown — PSI, PSII, carboxysome, ATP synthase, Calvin cycle all collapsed. But this is day 14 — cells are near death. Is this N-stress or dying?

## N-metabolism genes across all treatments

Queried 17 known N-metabolism genes (amt1, cynABDS, glnA, glnB, glsF, ntcA, ureABC, urtABCDE) across all MED4 experiments.

### Mirror pattern: N-starvation vs coculture

| Gene | Function | N-stress | Coculture |
|------|----------|----------|-----------|
| ntcA | Global N regulator | UP (3.4-4.5) | DOWN (-11.5) |
| glnA | Glutamine synthetase | UP (3.4-5.1) | DOWN (-2.8) |
| amt1 | Ammonium transporter | UP (2.3-3.9) | — |
| cynABDS | Cyanate transport/hydrolysis | UP (2.5-5.7) | DOWN (-5 to -13) |
| urtABCDE | Urea transport | UP (2.3-5.4) | DOWN (-3 to -13) |
| ureABC | Urease | UP (2.3-3.1) | DOWN (-4.6) |

**Interpretation:** Alteromonas provides enough N that MED4 doesn't need to scavenge. Consistent with nitrogen recycling thesis.

### RNA vs protein divergence in coculture (living cells, 90 days)

| Gene | Protein trend (day 18→89) | RNA trend |
|------|---------------------------|-----------|
| ntcA | UP throughout (1.8→1.2) | Trending down, sig at day 89 (-2.1) |
| glnA | UP throughout (2.6→2.2) | DOWN (-1.7 to -2.2) |
| glnB | UP, increasing (1.4→3.3) | DOWN (-2.0 to -2.6) |
| cynABD | UP (2.5→3.9) | Mostly no RNA signal |
| cynS | UP, increasing (1.7→4.7) | Modest RNA up at day 89 (1.6) |
| urtA | UP throughout (1.7→2.5) | Down at day 60 (-2.0) |
| ureABC | UP throughout | RNA not significant |

**Interpretation:** Cells made these proteins early and they persist. Don't need to keep making mRNA. Post-transcriptional survival strategy.

## Tolonen time course dynamics

| Timepoint | Sig genes | Up | Down | Pattern |
|-----------|-----------|-----|------|---------|
| 0h | 13 | 5 | 8 | Baseline noise |
| 3h | 6 | 6 | 0 | Pure N-specific |
| 6h | 236 | 113 | 123 | Major response wave |
| 12h | 307 | 145 | 162 | Peak response |
| 24h | 216 | 75 | 141 | Skewed DOWN — shutdown |
| 48h | 90 | 48 | 42 | Collapse |

## 3h first responders — the N-specific signature

Only 6 genes in MED4's genome respond within 3h:

| Gene | Locus tag | Product | log2FC at 3h |
|------|-----------|---------|-------------|
| cynA | PMM0370 | Cyanate transporter | 4.31 |
| urtA | PMM0970 | Urea transporter | 3.91 |
| glnA | PMM0920 | Glutamine synthetase | 3.76 |
| — | PMM0958 | Hypothetical (DUF1830) | 4.82 |
| — | PMM1462 | Hypothetical (adjacent to glnB) | 4.73 |
| — | PMM0687 | Hypothetical | 2.68 |

### Cross-stress specificity of 3h first responders

| Gene | N-stress | Iron | Phosphorus | Carbon | Light | Viral |
|------|----------|------|-----------|--------|-------|-------|
| cynA | UP always | — | — | — | — | — |
| PMM0958 | UP always | — | — | — | — | — |
| PMM1462 | UP always | — | — | — | — | — |
| PMM0687 | UP always | — | — | — | — | — |
| urtA | UP | — | DOWN | mixed | DOWN | — |
| glnA | UP | — | — | — | DOWN | — |

**Conclusion:** 4 of 6 first responders (cynA, PMM0958, PMM1462, PMM0687) are exclusively nitrogen-specific across all MED4 experiments in the KG. The remaining 2 (urtA, glnA) are primarily N-specific but have minor responses to other stresses.

## Working hypothesis

N-deprivation response proceeds in phases:
1. **N-specific** (0-3h): N-scavenging genes only (6 genes)
2. **Nutrient stress** (6-12h): Broader metabolic restructuring
3. **General stress** (12-24h): Photosynthesis shutdown, HLI activation
4. **Collapse/death** (24-48h): Declining gene count

Coculture rescues cells at phase 1-2, allowing a stable survival state.

## Next steps

- Classify 6h/12h/24h/48h genes by N-specificity (cross-reference with other stresses)
- Compare Tolonen 48h to Weissberg axenic day 14
- Determine whether HLI activation is N-specific or general stress
- Check what's unique to coculture survival state
