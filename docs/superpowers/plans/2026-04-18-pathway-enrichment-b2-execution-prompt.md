# B2 execution — prompt for a new Claude Code session

Copy the block below into a fresh Claude Code chat at the `multiomics_research` repo root (`/home/osnat/github/multiomics_research`). Do not summarize or edit the block — paste verbatim.

---

## Prompt to paste

You are continuing work on the `multiomics_research` project. Your task is to execute the pathway enrichment B2 analysis end-to-end following the written plan. This is a research execution session, not a planning session — the spec and plan are locked.

**Primary repo:** `/home/osnat/github/multiomics_research`
**Sibling repo (MCP server source):** `/home/osnat/github/multiomics_explorer` (editable install via `uv.sources`)

### Required reading, in order

1. **Load the `research-methodology` skill** via the Skill tool. This is non-negotiable — it defines the step-protocol (do → show → explore → decide with two commits per step and three hard gates), notebook discipline, artifact layout, and the eight rules (KG-only data source, locus tags over gene names, source tagging, etc.). Per CLAUDE.md: load this BEFORE touching code.
2. **Read the spec:** [`docs/superpowers/specs/2026-04-18-pathway-enrichment-b2-design.md`](specs/2026-04-18-pathway-enrichment-b2-design.md). Pay attention to §5.0 (step-protocol pace), §7 (methodology M1/M2/M3), §8 (risks), §10 (meta-deliverables).
3. **Read the plan:** [`docs/superpowers/plans/2026-04-18-pathway-enrichment-b2-plan.md`](2026-04-18-pathway-enrichment-b2-plan.md). Tasks 0–14 are your execution list. Every task is a bite-sized step; do not collapse do/show/explore/decide into a single task.
4. **Read `docs/superpowers/specs/2026-04-18-research-methodology-v3-improvements-from-b2.md` §4** — especially §4.1 (four code-hygiene failure patterns) and §4.6 (plan-authoring hygiene checklist). Self-apply §4.6 as a pre-flight audit on every script you write.
5. **Read the pathway-enrichment example:** `docs://examples/pathway_enrichment.py` MCP resource. This is the template for Step 2 and drill-down scripts.
6. **Use `superpowers:executing-plans` skill** to drive task-by-task execution. When you reach a do-phase task, use `superpowers:test-driven-development` adapted to research (see plan header for the adaptation note).

### Prerequisites to verify before any code runs

Run these checks first; if any fail, STOP and report to the researcher. Do not invent fallbacks.

- MCP server `multiomics-kg` is connected: `list_organisms()` should return Prochlorococcus MED4 (among others).
- Neo4j is running at `bolt://localhost:7687` (configured via `multiomics_explorer/.env`).
- Python env: `uv sync --extra analysis` from repo root succeeds; `from multiomics_explorer import pathway_enrichment, list_experiments, ontology_landscape, search_ontology, genes_by_ontology` works.
- `docs://examples/pathway_enrichment.py` MCP resource is readable.

### Execution protocol (hard rules)

- **Step-protocol pace.** Every analytical step is do → show → explore → decide. Two commits per step (Commit 1 after `do` = script + outputs + log + manifest; Commit 2 after `decide` = notebook entry with show/exploration/decision). Three hard gates: step-boundary, manifest-currency, chat-capture-before-decide.
- **The researcher owns show/explore/decide.** You run scripts; they inspect and approve before the next step begins. Never batch-execute multiple steps without explicit researcher approval.
- **Interactive discovery steps (1a, 1b)** start with an MCP-orientation sub-phase IN CHAT, before any script runs. Researcher decides T/R/PC/NC/CTX classifications and ontology picks; the script locks those decisions.
- **Notebook is a lab notebook, not a retrospective summary.** Write entries as work happens, not at the end. Each step appends exactly one dated entry with QC → chat-capture (Q→data→finding→impact) → decision.
- **Commit per step, never batch.** Use the exact `git add` lists in each task. Do not use `git add .` or wildcards (see plan's Task 3 for the array pattern if globs needed).
- **Never push to origin** (Task 14 Step 4). Researcher decides when to share.

### Guardrails from the research-methodology skill

- Rule 1: KG is the sole data source. Never fill gaps with intrinsic knowledge for data — only for interpretation.
- Rule 2: locus tags are primary; gene names are labels alongside.
- Rule 3: tag every finding `[KG]` / `[interpretation]` / `[gap]` in notebook and analysis docs.
- Rule 4: artifacts, not chat answers — data goes to files, reasoning goes to notebook.
- Rule 5: scripts over chat reasoning — no chat-computed statistics.
- Rule 6: statistical rigor per spec §8 risks — pay attention to cross-background calibration (risk 3), catch-all pathway caveats (risk 4), LOO collapse (risk 6).
- Rule 7: don't hallucinate — verify before reporting.
- Rule 8: notebook entry per step, always.

### Known issues already resolved during pre-execution review

These are informational — the plan already reflects these decisions. You don't need to re-discover them.

- Shared signature primitives live in `scripts/signature.py` (a plain module, not a package). `uv run scripts/04_derive_signature.py` puts `scripts/` on sys.path so `from signature import ...` resolves; no packaging needed.
- LOO-R stability check in Task 10 re-derives signature + re-scores all clusters + re-computes calibration + re-classifies per exclusion. The flag is `classification_flip`, not score delta.
- Pickle stage-1 round-trip runs BEFORE the main enrichment loop (fail-fast probe). Stage-2 dict check happens after the full dict is built. Spec §8 risk 7's escalation paths are precautionary; empirical round-trip passed on 2026-04-19.
- `$ANALYSIS_DIR` is persisted to `.analysis_dir` in Task 0. Every later task's first two lines must be `cd /home/osnat/github/multiomics_research && export ANALYSIS_DIR="$(cat .analysis_dir)"`.
- `list_experiments` lacks an `experiment_ids` parameter; Task 1 uses `list_experiments(verbose=True, limit=None)` + local filter (API gap is logged in api_coverage.md pre-seed).

### First move

After loading the skill and reading the spec + plan, tell the researcher:

> "Ready to execute B2. I've loaded the research-methodology skill, read the spec + plan + v3 meta-doc §4, and understand the step-protocol. Should I run prerequisite checks and proceed to Task 0 (scaffolding)?"

Wait for the researcher's go-ahead before Task 0. Do not skip this confirmation.

---

## If something goes wrong during execution

- If you hit a blocker not covered by the plan, stop and report to the researcher. Do not improvise around it.
- If a spec §8 risk fires, follow the contingency in order. Do not skip the empirical probe when one is specified.
- If you catch yourself collapsing do/show/explore/decide into a single task, stop — that's a spec §5.0 violation.
- Record all friction (KG data bugs, MCP friction, skill drift, API surface issues) to `analyses/<dir>/gaps_and_friction.md` as it happens, not retrospectively. The file is pre-seeded in Task 12 Step 5 with the baseline from pre-execution review — append the execution-time delta only.
