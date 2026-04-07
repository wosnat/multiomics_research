# KG rules

Extends SKILL.md Rule 1 (KG is the sole data source) with
operational details for scoping, gap handling, and API choice.

## Common KG gaps
- Missing organisms or strains not yet in the KG
- Annotations not loaded (e.g., some GO terms, regulatory elements)
- Expression data not available for a condition/study
- No adjusted p-values in the original dataset
- Missing metadata (timepoints, replicates, normalization method)
- No physiological data (growth curves, Fv/Fm, cell counts)

## Separating KG from literature

Clearly separate KG-derived findings from literature context:
- Use a distinct format: "**From the KG:** ... **From the
  literature (not verified in this KG):** ..."
- Never cite a specific paper from intrinsic knowledge as if
  it were a KG reference
- If the user needs literature support, say so — this system is
  a KG explorer, not a literature review tool

## Scoping: what to verify

Before diving into analysis, confirm the KG has what you need:
- Are the needed organisms present? (`list_organisms`)
- Are there experiments for the conditions of interest?
  (`list_experiments`)
- Are relevant publications available? (`list_publications`)
- Is the gene set complete? (check `total_matching` vs `returned`)
- Document any gaps: "The KG does not contain [X] data for [Y]"

## MCP vs Python API

| Use MCP when | Use Python API when |
|---|---|
| Browsing and orienting | Extracting full datasets |
| Checking gene counts and summaries | Running enrichment or statistics |
| Navigating between tools | Producing tables and figures |
| Quick lookups (< 50 results) | Any computation on gene lists |
| Working in chat mode | Result set > MCP limit |
| Cross-experiment overview (`gene_response_profile`) | Reshaping or aggregating data |
| | Building response matrices (`response_matrix()`) |
| | Comparing gene sets (`gene_set_compare()`) |

### When to switch from MCP to utilities

If you need to reshape, aggregate, or compute on MCP results, use
the analysis utilities first. Write scripts only for project-specific
logic not covered by the utilities. Do not attempt to work around
MCP limits in chat with one-off parsing.
