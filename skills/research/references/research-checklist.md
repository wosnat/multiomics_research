# Research checklist

Step-by-step protocol for answering a biological research question
using the multiomics KG.

---

## Before starting

- [ ] Restate the research question precisely
- [ ] Identify what kind of question it is:
  - **Lookup** — gene identity, function, annotations → MCP only
  - **Survey** — one gene across conditions, or one condition across
    genes → MCP + maybe extraction
  - **Comparative** — two conditions, two organisms, pathway
    comparison → extraction + statistics required
  - **Discovery** — enrichment, clustering, unbiased screen →
    full extraction + computational pipeline required
- [ ] For comparative/discovery: create the analysis directory
  structure (see artifacts guide)

---

## Orientation (run once)

### Scope

- [ ] `list_organisms` — confirm target organisms are in the KG
- [ ] `list_publications` — find relevant studies
- [ ] `list_experiments(summary=true)` — understand available
  conditions, time points, omics types
- [ ] `list_experiments` with filters — get specific experiment IDs

### Gene identification

- [ ] Identify genes of interest using appropriate discovery tools:
  - By name/keyword: `resolve_gene`, `genes_by_function`
  - By pathway/ontology: `search_ontology` → `genes_by_ontology`
  - By homology: `search_homolog_groups` → `genes_by_homolog_group`
  - By prior knowledge: start with `resolve_gene` for known names
- [ ] Resolve all gene names to locus tags
- [ ] Check for paralogs — multiple locus tags per gene name?
  - If yes: `gene_overview` on all locus tags, treat separately
- [ ] `gene_overview` — confirm annotation types, expression data
  availability, ortholog coverage

### Orientation outputs

- [ ] Create `methods.md` stub (research question, data scope, gene
  set)
- [ ] Write first exploration log
  (`exploration/YYYY-MM-DD-orientation.md`)
- [ ] Initialize `gaps_and_friction.md`

### Gate: scope check

- [ ] Are the needed organisms present?
- [ ] Are there experiments for the conditions of interest?
- [ ] Are all genes resolved to locus tags? Paralogs identified?
- [ ] Is the gene set complete? (check `total_matching` vs
  `returned`)
- [ ] Document any gaps: "The KG does not contain [X] data for [Y]"

---

## Research loop (per iteration)

Run this checklist for each iteration of the research loop.

### Question

- [ ] State a testable claim or comparison (not open-ended)
- [ ] Mark type: `hypothesis`, `exploratory`, or `follow-up`

### Explore

- [ ] Use MCP for browsing, Python scripts for computation
- [ ] For cross-experiment comparison: use `gene_response_profile`
  or extract via API and aggregate in a script
- [ ] Save any ad hoc computation to `scripts/explore_*.py`
- [ ] If comparing across platforms: note caveat, compare direction
  and rank not magnitude

### Log

- [ ] Write exploration log entry
  (`exploration/YYYY-MM-DD-{topic}.md`) **during** the iteration
- [ ] Tag every finding: `[KG]`, `[interpretation]`, or `[gap]`

### Assess

- [ ] Classify findings: `established`, `preliminary`, or
  `speculative`
- [ ] Self-check: did I use knowledge not from the KG? Is it tagged?
- [ ] Append any gaps/friction to `gaps_and_friction.md`
- [ ] Update `methods.md` if established findings changed scope
- [ ] Decide: next question, or conclude?

---

## Synthesis (after the loop)

- [ ] Finalize `methods.md` — coherent, publication-ready
- [ ] Write or update `README.md` — summary, findings, file index,
  exploration log links
- [ ] Produce publication artifacts in `results/`
- [ ] Populate `references.md` or `references.bib`

---

## Review (self or peer)

- [ ] Every number traces to a named KG query or script output
- [ ] Findings are tagged with source (`[KG]`, `[interpretation]`,
  `[gap]`) in exploration logs
- [ ] No gene name used without its locus tag
- [ ] Paralogs and orthologs handled correctly
- [ ] Truncation metadata respected (not counting from limited rows)
- [ ] Findings classified (established / preliminary / speculative)
- [ ] Statistical methods appropriate and documented
- [ ] KG gaps and limitations stated in exploration logs and
  `gaps_and_friction.md`
- [ ] All artifact files exist and are consistent with the narrative
- [ ] Exploration logs form a followable chain (`## Next` links)
- [ ] Claims match the data — no over-interpretation
