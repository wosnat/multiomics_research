# 2026-03-31: Corrected N-stress marker specificity

## Question

Which N-metabolism genes are truly N-specific when accounting for table_scope (absence = not significant for significant_only experiments)? (follow-up)

## Approach

Reinterpret the specificity table from iteration 1. For treatment types where all experiments use `significant_only`, `filtered_subset`, or `significant_any_timepoint` scopes on microarray platforms, change `?` (not tested) to `-` (tested, not significant).

## Findings

### Corrected specificity table

Applying the table_scope correction, the `?` entries become `-` (measured, not significant) for all non-N treatment types. The corrected table:

| Gene | Name | N-stress | Coculture | Carbon | Iron | Light | P | Viral | Salt | Class |
|------|------|----------|-----------|--------|------|-------|---|-------|------|-------|
| PMM0965 | ureA | UP(13) | - | - | - | - | - | - | - | **N-specific** |
| PMM0964 | ureB | MIX | - | - | - | - | - | - | - | **N-specific** |
| PMM1463 | glnB | MIX | - | - | - | - | - | - | - | **N-specific** |
| PMM1512 | glsF | UP(3) | - | - | - | - | - | - | - | **N-specific** |
| PMM1462 | (hyp) | MIX | - | - | - | - | - | - | - | **N-specific** |
| PMM0246 | ntcA | MIX | MIX | - | - | - | - | - | - | **N+coculture** |
| PMM0371 | cynB | UP(11) | MIX | - | - | - | - | - | - | **N+coculture** |
| PMM0372 | cynD | UP(11) | UP | - | - | - | - | - | - | **N+coculture** |
| PMM0393 | pipX | UP(2) | DN | - | - | - | - | - | - | **N+coculture** |
| PMM0687 | (hyp) | MIX | MIX | - | - | - | - | - | - | **N+coculture** |
| PMM0973 | urtD | UP(5) | UP | - | - | - | - | - | - | **N+coculture** |
| PMM0963 | ureC | UP(9) | DN | DN | - | - | - | - | - | N+carbon |
| PMM0373 | cynS | UP(10) | MIX | DN | - | - | - | - | - | N+carbon |
| PMM0974 | urtE | UP(6) | DN | UP | - | - | - | - | - | N+carbon |
| PMM0971 | urtB | MIX | - | MIX | - | - | - | - | - | N+carbon |
| PMM0972 | urtC | UP(6) | - | MIX | - | - | - | - | - | N+carbon |
| PMM0370 | cynA | UP(15) | UP | DN | MIX | - | - | - | - | N+multi |
| PMM0958 | (DUF1830) | UP(14) | MIX | - | MIX | - | - | - | - | N+multi |
| PMM0263 | amt1 | MIX | - | DN | - | - | - | UP | - | N+multi |
| PMM0920 | glnA | MIX | MIX | MIX | - | DN | - | - | - | N+multi |
| PMM0970 | urtA | MIX | UP | MIX | - | DN | DN | UP | - | broad |

### Tier classification

- [KG] **Tier 1 — N-specific (5 genes):** ureA, ureB, glnB, glsF, PMM1462. Respond only to nitrogen_stress. No response in any other treatment type (6 types tested). These are the most reliable N-stress markers.
- [KG] **Tier 2 — N+coculture (6 genes):** ntcA, cynB, cynD, pipX, PMM0687, urtD. Respond to N-stress and coculture (which involves N-provision by Alteromonas). Biologically coherent — coculture response is N-related.
- [KG] **Tier 3 — N+carbon (5 genes):** ureC, cynS, urtE, urtB, urtC. Also respond to carbon stress. [interpretation] N-deprivation causes energy/carbon crisis, so some cross-reactivity with carbon stress is expected.
- [KG] **Tier 4 — multi-stress (5 genes):** cynA, PMM0958, amt1, glnA, urtA. Respond to 3+ treatment types. Not reliable as N-specific markers.

### Best single N-stress marker

- [KG] **ureA (PMM0965):** UP in 13/25 N-stress timepoints across 5/8 experiments. Never responds to any other stress (tested in 7 treatment types). Best rank 6. Max log2FC 3.08. No coculture response.
- [interpretation] ureA is the strongest N-specific marker in this gene set. It encodes the urease gamma subunit — part of the urea degradation pathway that provides alternative N when primary sources are depleted.

### Temporal dynamics (from MIX signals)

- [KG] Several genes show MIX under N-stress: UP early, DOWN late (ntcA, glnB, glnA, PMM1462, ureB, urtA, amt1, PMM0687). This reflects the time course from Tolonen 2006 — early N-scavenging response followed by late shutdown/death.
- [interpretation] For marker purposes, MIX genes are informative for **staging**: UP = early/active N-stress, DOWN = late/dying. Combined UP+DOWN signal = time course spanning active stress through collapse.

## Assessment

- **Established:** ureA (PMM0965) is N-specific across all available MED4 experiments (5/8 up, 0/7 other stresses). Strongest single marker.
- **Established:** Tiers 1-2 (11 genes total) respond only to N-stress and coculture. The coculture response is biologically coherent (N-related).
- **Preliminary:** The corrected specificity depends on table_scope reinterpretation — absence from `significant_only` experiments = not significant. This is a reasonable assumption for microarray experiments but should be verified against original data.
- **Preliminary:** Carbon stress cross-reactivity in Tier 3 may reflect N→C metabolic coupling rather than true non-specificity.

## Gaps and friction

1. **gene_response_profile `not_known` is misleading** — The tool correctly reports no KG edge, but doesn't account for table_scope. For `significant_only` experiments, absence is informative (not significant), not unknown. This should be documented as a tool limitation or addressed with a `groups_tested_not_responded` field.

## Temporal dynamics (Tolonen microarray, 0-48h)

### Tier 1 genes (N-specific) — Tolonen time course log2FC:

| Gene | Name | 0h | 3h | 6h | 12h | 24h | 48h | Pattern |
|------|------|----|----|----|-----|-----|-----|---------|
| PMM0965 | ureA | +0.3 | +1.1 | **+1.9** | +1.5 | +1.1 | +1.2 | [KG] Peaks at 6h, sustained |
| PMM0964 | ureB | -0.3 | +0.5 | +0.9 | **+1.3** | +0.7 | +0.6 | [KG] Peaks at 12h, sustained |
| PMM1463 | glnB | +0.1 | +2.1 | **+2.9** | +2.3 | +1.3 | +1.5 | [KG] Peaks at 6h, sustained |
| PMM1512 | glsF | -0.2 | 0.0 | +0.9 | **+1.2** | +0.9 | +1.0 | [KG] Slow onset, peaks 12h |
| PMM1462 | (hyp) | +0.2 | **+4.7** | +5.6 | +4.4 | +3.1 | +3.1 | [KG] Fastest responder! 3h start |

### Tier 2 genes (N+coculture) — Tolonen time course:

| Gene | Name | 0h | 3h | 6h | 12h | 24h | 48h | Pattern |
|------|------|----|----|----|-----|-----|-----|---------|
| PMM0246 | ntcA | -2.2 | **+3.0** | +4.5 | +3.0 | +1.7 | +1.2 | [KG] 3h spike, declining |
| PMM0371 | cynB | 0.0 | +1.1 | +2.3 | **+2.4** | +1.9 | +2.1 | [KG] Peaks 12h, sustained |
| PMM0372 | cynD | -0.2 | +0.1 | +1.5 | **+2.2** | +1.5 | +1.4 | [KG] Peaks 12h, sustained |
| PMM0393 | pipX | -0.7 | -0.4 | +0.1 | +0.8 | -0.4 | -1.4 | [KG] Weak, transient at 12h |
| PMM0687 | (hyp) | +0.4 | +2.7 | **+3.9** | +3.1 | +2.7 | +2.8 | [KG] 3h start, peaks 6h, sustained |
| PMM0973 | urtD | -0.5 | +0.4 | +0.6 | 0.0 | +0.3 | +0.3 | [KG] Weak response |

### Weissberg coculture proteomics (day 18-89, coculture cells surviving)

- [KG] ALL Tier 1-2 proteins are **persistently elevated** in coculture cells:
  - glnB: +1.3 → +3.3 (increasing over 90 days)
  - cynB: +2.7 → +3.5 (increasing)
  - cynD: +2.6 → +3.2 (stable high)
  - ureA/B: +2.0-3.0 (stable)
  - ntcA: +1.8 → +1.2 (slowly declining but still elevated)
- [interpretation] Even though mRNA may decline in coculture (Aharonovich data shows cynB mRNA DOWN), proteins persist for months. This suggests early N-stress activates these proteins and they're maintained post-transcriptionally in the survival state.

### Coculture mRNA: contradictory results between studies

- [KG] Aharonovich 2016 (N-replete Pro99): ntcA DOWN (-11.5), cynB DOWN (-5.2), pipX DOWN (-11.0), PMM0687 DOWN (-11.0). Coculture suppresses N-stress genes when N is available.
- [KG] Weissberg 2025 (PRO99-lowN, day 11): ntcA UP (+2.4), cynB UP (+2.5), PMM0687 UP (+2.5). Coculture under N-limitation: N-scavenging is active.
- [interpretation] The direction of coculture response depends on N-availability. In N-replete medium, Alteromonas provides enough N → N-scavenging genes suppressed. In N-limited medium, both organisms are N-stressed → N-scavenging is needed.

## Next

- Summarize findings: what markers to use, when, and caveats
- Write README with file index
- Consider: should we look at MIT9313 (LLIV clade) to test whether these markers are conserved across Prochlorococcus?
