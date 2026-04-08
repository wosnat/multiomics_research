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

**Load the `research-methodology` skill BEFORE brainstorming.** It
contains the rules for KG usage, gene identity, artifact structure,
notebook discipline, and anti-hallucination. These rules shape the
analysis design — loading after the spec means retrofitting.

For specific analysis types, invoke the corresponding recipe skill
(e.g., `enrichment`, `response-matrix`, `conservation`).

Use the superpowers workflow for research: brainstorm the question,
write a plan, execute with checkpoints. The methodology skill
provides the domain rules; superpowers provides the process
discipline.

### Naming conventions

- **Spec steps ≠ plan tasks.** Spec steps describe the analytical
  pipeline (Step 1, Step 2, ...). Plan tasks describe
  implementation work and may split one spec step into multiple
  tasks. Specs must not reference files by task number. Plans must
  cross-reference the spec step in each task description (e.g.,
  "Task 5: Build signature utilities [Step 3]"). This prevents
  confusion when the plan introduces scaffolding or verification
  tasks that don't correspond to a spec step.

### Process overrides

- **Don't skip subagent reviews** for tasks that produce data
  outputs. At minimum, run spec compliance review. The tasks that
  seem simple are where silent bugs hide.

- **Notebook-commit gates in plans.** When writing implementation
  plans for research analyses, every plan must include a
  notebook-commit gate between data-producing steps. The gate is:
  "commit notebook entry for Step N before beginning Step N+1."
  This ensures the executing agent treats the notebook as a
  blocking dependency, not a nice-to-have.

- **Interactive discovery steps.** Rule 5 (scripts over chat
  reasoning) has a carve-out for interactive discovery/scoping
  steps — see the research-methodology skill for the
  frozen-output + notebook-entry pattern.

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
