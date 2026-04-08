# Gaps and friction points

## KG gaps

1. **No well-powered non-N stress experiments for MED4.** Phosphate (34 genes), iron (112 genes), salt (28 genes), viral (25-45 genes) — all too sparse for signature-based scoring. Carbon experiments have no non-significant rows. The only usable non-N control (Steglich high white light) has 198 genes. Limits our ability to test N-specificity of the signature.

2. **Read 2017 filtered_subset scope.** Top 50% of genes by expression — 857 vs 1,697 genes in Tolonen. 73 signature-relevant genes absent from Read may be expressed but below the threshold. No way to recover them without the full dataset.

3. **Experiment metadata quality is excellent.** Light conditions, temperature, medium, statistical test — all present in `list_experiments(verbose=True)`. Enabled the diel concern resolution without leaving the KG. Diel cluster data (Zinser 2009) was valuable for investigation even though it didn't confound the analysis.

## MCP/API friction

1. **`list_experiments` doesn't return matching_genes.** Had to call `differential_expression_by_gene(summary=True)` per experiment to get gene counts — N+1 queries during scoping. Adding `matching_genes` to `list_experiments` results would save significant time.

2. **`resolve_gene` parameter name.** The parameter is `identifier`, not `search_term`. Error message was clear but the name isn't intuitive for first-time use.

3. **`experiments_to_dataframe()` use case unclear.** Available but not used — we went directly to `to_dataframe()` for experiment listing. The when-to-use-which could be documented more clearly.

## MCP/API — what worked well

1. **`to_dataframe()` worked reliably.** No flattening issues throughout the analysis.

2. **`verbose=True` puts table_scope on every row.** Exactly what we needed — no fallback to envelope required.

3. **table_scope field makes dataset limitations transparent.** Read's `filtered_subset` limitation is visible in every row, enabling correct downstream handling.

## Methodology insights

1. **Rank normalization matters.** Normalizing by total_genes (v1 approach) compresses all ranks to ~0.95+, making rank score ≈ hit rate. Normalizing by n_significant_in_direction gives meaningful spread. Discovered during interactive review — would have been missed without the walkthrough.

2. **Early timepoint exclusion is symmetric.** Both Tolonen 0h/3h and Read 3h show transient responses. v1 excluded only Tolonen's early timepoints. Excluding Read 3h fixed the rbcL misclassification.

3. **Proteomics gene coverage bias is systematic.** Missing 42 genes aren't random — they're small/membrane proteins (hli, photosystem subunits, ribosomal) with better-than-average ranks. This biases proteomics scoring toward larger, more detectable proteins.

4. **RNA/protein discordance is genuine DE difference, not coverage artifact.** Tested by scoring RNA-seq with only the 147 proteomics-detectable genes. RNA-seq axenic still scored 0.548 vs proteomics 0.066 for the same gene set.

## Process retrospective

### What worked

1. **Methodology skill loaded before brainstorming.** Artifacts structure, notebook discipline, and source tagging shaped the design from the start. v1's biggest problem fixed.

2. **Interactive exploration caught real issues.** rbcL discordance, Read 3h transient, rank normalization compression, proteomics coverage bias, RNA/protein discordance test — all discovered through walking through specific genes with the researcher. None of these would have been caught by running the pipeline end-to-end.

3. **Researcher's questions drove the best discoveries.** "Is it the missing genes or the DE?" → restricted-gene-set test confirming genuine discordance. "The score is dominated by hit rate" → normalization fix. "What about diel?" → light conditions investigation. "What is the rank distribution of missing proteomics genes?" → systematic coverage bias finding.

4. **Toy-data verification before real data.** 31 tests caught the normalization change cleanly. When we changed the formula, tests told us immediately what broke and why.

5. **Single metric decision.** One metric understood thoroughly > three black boxes. When the metric had a problem (normalization), we caught it because we only had one thing to reason about.

6. **Controls as a design principle.** Classifying every experiment upfront (reference/negative/target) gave us something to validate the metric against before trusting it on targets.

7. **sig_utils separation.** Clean boundary between reusable methodology and analysis-specific scripts. The extraction utility avoided duplicated API logic. Signature and scoring modules are genuinely packageable.

8. **Three scoring tiers.** Top (38) > Core (189) > Extended (531) consistently — confirms the signal concentrates in the strongest markers and extended genes dilute.

### What didn't work

1. **Notebook fell behind after Step 3.** Steps 4-6 have no notebook entries. The richest exploration (proteomics coverage bias, rank normalization discussion, RNA/protein discordance test, photosynthesis gene walkthrough) exists only in the chat conversation. A future reader has methods.md and decisions.md but not the messy investigative trail.

2. **Script 01 never written.** The discovery/scoping step was done entirely via interactive MCP queries and inline Python. The scoping CSV exists but the process isn't reproducible from a script. If someone wants to rescore with updated KG data, they can't re-run Step 1.

3. **Log for Step 1 missing.** No `logs/01_discover_experiments.log`.

4. **Manifests written retroactively.** DATA_MANIFEST was written at the end, not after each extraction step as the methodology specifies.

5. **Steglich light control too sparse.** Should have checked gene counts before including. The 198-gene filtered_subset matched only 43 of 189 signature genes — unusable for signature-based scoring. The score (0.400) is misleading.

6. **`__pycache__` directories in the analysis tree.** Should be gitignored.

7. **`--explore` output not captured in logs.** Scripts have `--explore` flags that print marker gene traces to stdout, but this output went to chat, not to log files. The logs captured during scoring are sparse (summary lines only).

8. **Plotting cosmetics.** Trajectory x-axis in hours (hard to read as days), cramped x-labels in tier comparison. Affects publication readiness.

### Proposed changes

**To the research-methodology skill:**

1. **Hard gate: notebook entry committed before next step begins.** The skill says "capture key data points and conclusions in the notebook" but doesn't enforce when. We did it for Steps 1-3 and then stopped under momentum. A hard gate ("no step proceeds without a committed notebook entry for the previous step") would prevent this.

2. **Guidance for interactive vs scripted steps.** Discovery/scoping steps are naturally interactive — you browse, discuss, classify. Forcing into a script feels wrong. But leaving unscripted means unreproducible. Proposed pattern: "interactive steps produce a frozen output file (CSV) that downstream scripts consume, plus a notebook entry documenting the reasoning. Optionally, write a script that reproduces the frozen output."

3. **Log verbosity standard.** Scripts write to `logs/` but the skill doesn't say what detail level. The `--explore` output (marker gene traces, QC diagnostics) should go into the log file, not just stdout. Proposed: "logs capture everything the researcher needs to verify the step without rerunning it."

4. **Manifest update timing.** "Updated with every extraction step" is correct but clashes with interactive flow. Practical rule: "manifest updated in the same commit that adds the new data file."

5. **Methodological revision pattern.** When the metric/formula needs to change mid-analysis (as happened with rank normalization), the skill should describe the rerun workflow: update the utility, rerun tests, rescore everything, add a rerun entry to the notebook explaining what changed and why. We did this correctly but ad hoc.

6. **Working example for formula specs.** The statistical rigor reference says "every formula must include worked examples." This was followed in the brainstorming (rank score worked example) but the spec's formula description was still somewhat abstract. The toy test data ended up being the real worked examples — perhaps the skill should say "toy test = worked example."

**To the MCP/KG:**

1. **Add `matching_genes` to `list_experiments` results.** Saves N+1 queries during scoping. The information is already computed — just not exposed.

2. **Consider a `docs://analysis/signature_scoring` guide.** The pattern we built (extract → summarize per gene → intersect references → apply signature → score → permute) is general enough to document as a tutorial. Would help future users doing similar cross-study analysis.

3. **Rename `resolve_gene` parameter to `search_term` or `query`.** `identifier` is technically correct but not what a user types first.

**To the process (superpowers workflow):**

1. **Notebook entries should block next step.** The executing-plans skill should check for a notebook entry commit before marking a step complete. This is the gate that slipped.

2. **Interactive discovery steps need a pattern.** The plan had a script for Step 1, but we skipped it because interactive exploration was more natural. The plan should distinguish "interactive steps (produce output CSV + notebook entry)" from "scripted steps (produce output + log + notebook entry)."
