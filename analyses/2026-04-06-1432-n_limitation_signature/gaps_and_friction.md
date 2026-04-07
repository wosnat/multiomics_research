# Gaps and friction points

## KG data bugs

1. **Axenic RNA-seq timepoint label:** The Weissberg axenic RNA-seq experiment has `timepoint="single"` and `timepoint_hours=NaN` instead of a meaningful label (the experiment is day 11 starvation vs exponential). This caused it to be silently excluded from trajectory plots until manually patched. Other single-timepoint experiments likely have the same issue. **Action:** Consider labeling non-time-course experiments with a descriptive timepoint (e.g., "day 11") rather than null/single.

## KG gaps

1. **Read 2017 filtered_subset scope:** Read 2017 DE data is filtered to top 50% of genes by expression level (table_scope: filtered_subset). There is no easy programmatic way to know which genes were excluded from the filter vs. genuinely not DE. This affects extended signature classification: `tolonen_only_read_absent` genes could be truly absent from the study or just below the expression threshold.

2. **No gene product in signature-building path without verbose DE:** The `differential_expression_by_gene` default (non-verbose) output lacks `product` and `gene_category`. Since these are essential for interpreting signatures, every extraction script must remember to set `verbose=True`. If forgotten, unnamed genes remain opaque until a separate `gene_overview` call.

## MCP friction

1. **`limit=None` not obvious:** The Python API supports `limit=None` to get all results, but this isn't the default and isn't prominently documented. First instinct is to set an arbitrary high limit (e.g., 10000). The MCP tool docs mention it but it's easy to miss.

2. **No bulk "genes in experiment" query:** To build a signature, we needed all significant DE genes for an experiment — but also needed to know which genes were in the dataset but NOT significant (for the `tolonen_only_read_ns` classification). This required two calls: one with `significant_only=True` for the signature, and the full dataset to get the `study_b_all_locus_tags` set. A single call with a flag like `include_background_tags=True` would simplify this.

3. **Cross-experiment scoring requires manual timepoint grouping:** The scoring script had to manually loop over `de_df.groupby("timepoint")` because DE results come as a flat table. A utility like `score_gene_set_per_timepoint()` that takes a gene set + DE data + metric function would reduce boilerplate. This is a candidate for `multiomics_explorer.analysis`.

## Skill/methodology friction

1. **Research methodology skill loaded too late.** The brainstorming and spec were mostly complete before checking the skill. Artifacts structure, source tagging, and logging requirements were retrofitted. The skill should be loaded BEFORE brainstorming, not as a review step after.

2. **Spec formulas were unclear to the researcher.** The mean signed log2FC formula was written in shorthand ("reference direction as expected sign") that was clear to the implementer but not to the domain expert. The spec should include worked examples for every formula — this was only added after the user asked "what is reference direction?", "why log2fc and not rank?", "so we put abs fc?"

3. **Plan contained untested code.** ~800 lines of Python were written in the plan document without running any of it. The `intersect_de_lists` merge was broken (column collision), peak selection used log2fc instead of rank, study names were hardcoded. All caught by the user during review. Plans should contain pseudocode or function signatures, not implementation code — or the code should be tested before the plan is finalized.

4. **Subagent reviews were skipped.** The subagent-driven-development skill requires spec compliance review + code quality review after each task. These were skipped for speed. No bugs were caught, but this was luck — the gitignore issue and the axenic timepoint issue would have been caught by review.

## Process retrospective

### What worked

1. **Brainstorming was thorough and collaborative.** One question at a time, multiple choice format kept alignment. User pushback on formulas, rank vs log2fc, and study naming made the spec much better than the first draft.

2. **sig_utils reusable package.** Separating metrics from scripts paid off immediately — scoring core and extended signatures was trivial. The same functions will work for Approaches B/C with different gene sets.

3. **Toy data verification before real data.** User's call to slow down ("I think we are running too fast"). Built confidence in sig_utils code before committing to real KG extraction. No bugs found, but the verification itself was the value — hand-calculated expected values for every metric.

4. **Subagent dispatch for mechanical tasks.** Tasks 5-12 went smoothly. Parallel dispatch for independent tasks (9+10) saved time. The extraction and plotting scripts were well-specified enough that subagents needed no clarification.

5. **Incremental spec review.** User reviewed each section of the spec before moving on. Every review round surfaced real issues (rank vs log2fc, direction handling, study naming, logging requirements). This is slower but produces a much better spec.

### What didn't work

1. **Research methodology skill loaded as afterthought.** Artifacts, source tagging, logging all retrofitted after spec was mostly done.

2. **Formulas unclear to researcher.** Multiple rounds of "what does this mean?" questions that should have been prevented by including worked examples from the start.

3. **Plan code was untested.** Broken merge logic, wrong peak selection, hardcoded names — all found during review, none found by tests (because there were no tests until the user asked for them).

4. **Verbose DE data underused.** Extracted with `verbose=True` (product, gene_category in every row) but didn't join into signature CSVs until the very end. Did a redundant `gene_overview` + `gene_ontology_terms` lookup for information already in the data.

5. **Subagent reviews skipped.** Traded quality for speed. Got away with it this time.

6. **Axenic timepoint issue discovered late.** The "single" timepoint label caused the axenic RNA-seq point to be invisible in all trajectory plots. Found and fixed only after all 12 tasks were complete.

### Proposed changes

**To the skill (research-methodology):**

1. **Hard gate: load research-methodology skill before brainstorming.** Not after the spec, not as a review step. Before. The artifacts structure and logging requirements shape the entire analysis design.

2. **Require worked examples for every formula in specs.** If a spec includes a mathematical formula, it must include at least one worked example with concrete numbers showing input → computation → output. "Reference direction as expected sign" is not sufficient.

3. **Add a "join metadata early" rule.** When extracting DE data with `verbose=True`, immediately join `product` and `gene_category` into any gene-level output files (signatures, summaries). Don't defer annotation to a separate step.

4. **Add a "verify single-timepoint experiments" checkpoint.** After extraction, explicitly check for experiments where `timepoint` is null/single and confirm they're handled correctly in downstream scripts.

**To the MCP/KG:**

1. **Add `table_scope` to DE result rows.** Each DE row should carry its experiment's `table_scope` so scripts can distinguish `all_detected_genes` from `filtered_subset` without a metadata join.

2. **Label single-timepoint experiments meaningfully.** Instead of `timepoint=null` for non-time-course experiments, use the experimental timepoint (e.g., "day 11") from the study metadata. The current "single" label is a data loss.

3. **Add a `background_gene_count` field to DE results.** When querying significant genes, knowing the total background size (how many genes were tested) is essential for enrichment analysis and for computing normalized ranks. Currently requires a separate query.

4. **Consider a `score_gene_set()` analysis utility.** Takes a gene set (locus tags + directions) + experiment ID, returns per-timepoint activation scores. This is the core operation we built manually in `score_signature.py` and it's general enough to be a first-class utility in `multiomics_explorer.analysis`.

**To the process (superpowers workflow):**

1. **Don't write implementation code in plan documents.** Write function signatures, expected I/O, and pseudocode. Implementation happens in actual files where it can be tested. Plans with untested code create a false sense of completeness.

2. **Toy data verification as a mandatory step.** When building a reusable utility package (like sig_utils), verify with hand-calculated toy data before proceeding. This should be a plan task, not an ad-hoc intervention.

3. **Don't skip subagent reviews for "simple" tasks.** The tasks that seem simple are exactly the ones where silent bugs hide. At minimum, run spec compliance review for tasks that produce data outputs.

**To the research methodology skill — documentation requirements:**

Analysis directories should include these documents (in addition to the existing artifacts guide requirements):

1. **Data manifest** (`data/DATA_MANIFEST.md`): For each CSV, record row count, gene count, timepoints, which script produced it, and a one-line description. Column schema for recurring formats. Without this, files in `data/` are opaque after a few days.

2. **Results manifest** (`results/RESULTS_MANIFEST.md`): Same for results — what each file contains, which script produced it, figure descriptions.

3. **Decision log** (`decisions.md`): Design decisions with rationale. Not "what was done" (that's methods.md) but "why it was done this way and what alternatives were considered." Decisions made during brainstorming exist only in chat history, which is ephemeral. Critical for the researcher to understand their own analysis after a gap, and for reviewers/collaborators.

4. **Caveats for interpretation** (`caveats.md`): Distinct from gaps_and_friction (which is about process/tooling). Lists things a reader of the RESULTS needs to know before drawing conclusions — platform coverage differences, single-timepoint limitations, circularity in self-scoring, etc. This is the "fine print" that should accompany any publication figure or claim.

5. **Full scoring tables in exploration logs.** Not just prose summaries — include the actual numbers in markdown tables so a reviewer can scan all conditions at a glance without opening CSVs.
