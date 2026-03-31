# 2026-03-30: Tolonen timepoint classification

## Approach

Took all 108 significantly DE genes from the Tolonen 2006 N-deprivation time course (0-48h, microarray).
For each gene, queried all MED4 experiments in the KG (via Cypher) to find which other stresses it responds to.

Classification rules:
- **N-specific**: only responds to nitrogen_stress (and coculture, which is N-related)
- **Nutrient-stress**: also responds to iron, phosphorus, and/or carbon stress (but NOT light, viral, salt)
- **General-stress**: also responds to light, viral, or salt stress

Note: this analysis was done with jq/python one-offs parsing oversized MCP results — should be redone as a proper script using the Python API (see gaps_and_friction.md).

## Summary

| Timepoint | N-specific | Nutrient | General | Total |
|-----------|-----------|----------|---------|-------|
| 3h | 2 | 2 | 2 | 6 |
| 6h | 8 | 13 | 4 | 25 |
| 12h | 7 | 13 | 18 | 38 |
| 24h | 1 | 4 | 6 | 11 |
| 48h | 3 | 13 | 12 | 28 |
| **Total** | **21** | **45** | **42** | **108** |

## Key observations

### Phase 1: N-specific alarm (3h)
Only 6 genes. Two are purely N-specific (PMM0687, PMM1462 hypotheticals). The classic N-scavenging genes (cynA, urtA, glnA) also respond to other nutrient stresses but with much weaker effects.

### Phase 2: NtcA regulon + carbon/energy restructuring (6h)
25 new genes. N-specific wave includes **ntcA**, **glnB** (PII), **cynB**, **tatA**, and an HLI (PMM1390). Nutrient-stress wave dominated by **ATP synthase operon** (atpA/C/D/F/G/H) and **ribosomal proteins**, all shared with carbon stress. N-deprivation forces energy rebalancing.

### Phase 3: General stress cascade (12h)
38 new genes — largest wave. General-stress genes dominate (18/38): **HLI proteins**, **photosynthesis** (psaM), **Calvin cycle** (fbaA), more ribosomal proteins, **sigma factor** rpoD8.

### Phase 4: Chaperones and damage control (24h)
Only 11 new genes. **htpG** (N-specific chaperone), **groL2/groL1** (chaperonins), **EF-Tu**. Protein quality control — consistent with cellular damage.

### Phase 5: Photosystem collapse (48h)
28 new genes. Dominated by **photosystem shutdown**: psaJ/F/L/I/B/D (PSI), psbO/T/B/K (PSII), pcb, prkB. Most shared with carbon stress. General-stress additions are ALL HLI proteins (the strongest responders). Also **pstS** (phosphate transporter), **lexA/recA** (SOS/DNA damage).

## Biological interpretation

1. **0-3h: N-specific alarm** — Pure nitrogen scavenging. "Find more N."
2. **6h: Metabolic restructuring** — NtcA regulon activates. ATP synthase and ribosomes adjust. "Rebalance energy for low-N."
3. **12h: General stress** ��� Light-harvesting stress proteins, Calvin cycle. "The cell is broadly stressed."
4. **24h: Damage control** — Chaperones. "Proteins are misfolding."
5. **48h: Shutdown** — Photosystem collapse, HLI blast, DNA damage response. "The cell is dying."

The carbon_stress overlap is pervasive (45/108 genes). N-deprivation rapidly becomes a C/energy crisis.

## Implications for the paper rewrite

The Weissberg day 14 axenic RNA-seq (cells committed to death) most resembles Tolonen 48h: photosystem collapse, HLI activation, general stress. This is a dying/shutdown signature, not primarily N-stress.

The coculture survival state should be compared to Tolonen 3-6h: do coculture cells maintain the early N-specific response without progressing to shutdown?
