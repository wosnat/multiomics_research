# CLAUDE.md

## Project Overview

Research skills, evaluation framework, and usage logging for the multiomics KG MCP server.

This repo is the **consumer-facing** companion to [`multiomics_explorer`](https://github.com/osnat/multiomics_explorer), which builds the MCP tools. This repo:
- Packages research skills as a Claude Code plugin for researchers
- Evaluates how well Claude uses the MCP tools (trigger accuracy, chain quality, research output)
- Logs and analyzes real usage patterns to guide tool improvements

## Architecture

```
multiomics_research/            # This repo (consumer)
  skills/research-methodology/  # Domain rules for KG research (reference skill)
  skills/recipes/               # On-demand analysis protocols (one skill per method)
  hooks/                        # Usage logging hooks
  evals/                        # Evaluation cases and runners
  scripts/                      # Analysis and eval scripts
  benchmarks/                   # Eval results (gitignored)

multiomics_explorer/            # Sibling repo (builder)
  mcp_server/                   # MCP tools (18 tools)
  api/                          # Python API
  kg/                           # Query builders
  tests/                        # Tool correctness tests
```

## Dependencies

### MCP Server

This repo depends on the `multiomics-kg` MCP server from `multiomics_explorer`. The server must be running or accessible.

**Local development** (current):
```bash
# MCP server runs via uv from the sibling repo
uv run --directory /home/osnat/github/multiomics_explorer multiomics-kg-mcp
```

**Remote** (future): will point to a hosted MCP endpoint.

### Neo4j

The MCP server requires Neo4j at `bolt://localhost:7687` (or configured via `multiomics_explorer/.env`).

## Plugin Structure

This repo is structured as a Claude Code plugin. During development, load it with:
```bash
claude --plugin-dir /home/osnat/github/multiomics_research
```

For researchers, it will be installable from a marketplace (private git repo).

## Research methodology

**Load the `research-methodology` skill BEFORE invoking
`superpowers:brainstorming` for step 1 of an analysis.** The skill
contains KG usage rules, gene identity rules, anti-hallucination
patterns, scripts-over-chat-reasoning, and the 6-step research flow.
Loading after step 1 is committed means retrofitting.

### The 6-step flow

Every research analysis advances through 6 steps:

1. **Research question** — user prompt + clarifying questions →
   locked question (uses `superpowers:brainstorming` with overrides)
2. **KG entries** — relevant publications, experiments, organisms,
   data types
3. **Analysis framing** — selection + framing (hypothesis, controls,
   expected outcome) in KG terms
4. **Methods** — ad-hoc Python module using one item from step 3
   as a driving example
5. **Analyze** — run the method; produce scored outputs, figures,
   tables
6. **Evaluate** — assess against framing; harvest caveats; finalize
   paper

Steps 1–3 form the research proposal (locked at end of step 3).
Steps 4–6 execute against it.

### Intra-step rhythm: do → show → explore → decide

Every step advances through **do → show → explore → decide**. The
**decide** phase produces a minimal notebook.md checklist and pauses
for explicit researcher approval before committing. One commit per
step, at decide close. See
`skills/research-methodology/references/step-protocol.md` for commit
timing, the decide-gate checklist, and hard gates.

### Just-in-time formalization

Terms, predictions, metrics, stability checks, decisions, and
caveats enter the analysis **only when the data demands them**.
Nothing is enumerated in advance "just in case." No upfront
taxonomies of hypotheses, predictions, or thresholds. If you find
yourself listing things the analysis might need before the data has
arrived, stop.

### What is NOT used in this workflow

- **No `superpowers:writing-plans`** — the 6-step flow replaces the
  monolithic plan; each step's intent lives in its own `notebook.md`
- **No `superpowers:executing-plans`** — the decide-gate per step
  replaces task-level checkpoint reviews
- **No per-analysis `spec.md` / `plan.md` / `decisions.md` /
  `hypotheses.md`** — content lives in per-step `notebook.md` +
  single analysis-root `paper.md`
- **No cross-file labels** (H1 / P3 / D8 / T4) — labels are
  document-scoped and must be paired with a short readable name

On-demand tools that remain available: `superpowers:brainstorming`
(step 1), `superpowers:verification-before-completion`,
`superpowers:systematic-debugging`, `superpowers:requesting-code-review`.

## Evaluation Framework

Three evaluation layers:

| Layer | What it measures | How |
|-------|-----------------|-----|
| **Usage logging** | Real tool usage patterns | PostToolUse/PreToolUse hooks → JSONL |
| **Tool chain evals** | Does Claude pick the right tool sequence? | Research questions → expected chains → `claude -p` |
| **Skill quality evals** | Does the research skill improve output? | A/B benchmark: with-skill vs without-skill |

### Running evals

```bash
# Analyze logged usage data
python scripts/analyze_usage.py ~/.claude/logs/multiomics-kg-usage.jsonl

# Run chain evaluation (requires Neo4j + MCP server)
python scripts/run_chain_eval.py evals/chain_evals.yaml

# Run trigger evaluation (requires claude CLI)
python scripts/run_trigger_eval.py evals/trigger_evals.json --skill-path skills/research-methodology
```

### Logs

Usage logs are written to `~/.claude/logs/multiomics-kg-usage.jsonl` by the PostToolUse hook. Each line is a JSON object with tool name, parameters, timing, and response metadata.

## Key Files

| File | Purpose |
|---|---|
| `.claude-plugin/plugin.json` | Plugin manifest |
| `skills/research-methodology/SKILL.md` | Research methodology rules (reference skill) |
| `skills/recipes/` | On-demand analysis protocols |
| `hooks/hooks.json` | Usage logging hook configuration |
| `hooks/log-mcp-usage.sh` | Hook script that appends to JSONL |
| `evals/chain_evals.yaml` | Tool chain evaluation cases |
| `evals/trigger_evals.json` | Skill trigger evaluation cases |
| `evals/research_questions.json` | End-to-end research quality cases |
| `scripts/analyze_usage.py` | Usage log analysis |
| `scripts/run_chain_eval.py` | Chain evaluation runner |
| `scripts/run_trigger_eval.py` | Trigger evaluation runner |
