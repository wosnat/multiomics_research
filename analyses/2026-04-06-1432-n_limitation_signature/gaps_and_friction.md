# Gaps and friction points

## KG data bugs

## KG gaps

- **Read 2017 filtered_subset scope:** Read 2017 DE data is filtered to top
  50% of genes by expression level (table_scope: filtered_subset). There is
  no easy programmatic way to know which genes were excluded from the filter
  vs. genuinely not DE. This affects extended signature classification:
  tolonen_only_read_absent genes could be truly absent from the study or
  just below the expression threshold. [gap]

## MCP friction

## Skill/methodology friction

## Process retrospective

### What worked

- **sig_utils reusable package:** Building a small `sig_utils` package with
  standardized IO (`load_signature_csv`, `RESULTS_DIR`) made scripts
  composable and results reproducible across multiple analysis steps.

- **Toy data verification before real data:** Testing score computation on
  synthetic data with known expected outputs caught bugs in the rank score
  metric before applying it to real KG extracts.

### What didn't work

### Proposed changes

**To the skill:**

**To the MCP/KG:**

- **Add table_scope info to DE result rows:** DE result records should
  include the `table_scope` field in each row (not just as a filter option),
  so that downstream scripts can programmatically distinguish
  all_detected_genes records from filtered_subset records without requiring
  a separate metadata query. This would eliminate the Read 2017
  filtered_subset ambiguity in the extended signature classification.
