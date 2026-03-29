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
  skills/research/              # Research methodology skill
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
python scripts/run_trigger_eval.py evals/trigger_evals.json --skill-path skills/research
```

### Logs

Usage logs are written to `~/.claude/logs/multiomics-kg-usage.jsonl` by the PostToolUse hook. Each line is a JSON object with tool name, parameters, timing, and response metadata.

## Key Files

| File | Purpose |
|---|---|
| `.claude-plugin/plugin.json` | Plugin manifest |
| `skills/research/SKILL.md` | Research methodology skill |
| `hooks/hooks.json` | Usage logging hook configuration |
| `hooks/log-mcp-usage.sh` | Hook script that appends to JSONL |
| `evals/chain_evals.yaml` | Tool chain evaluation cases |
| `evals/trigger_evals.json` | Skill trigger evaluation cases |
| `evals/research_questions.json` | End-to-end research quality cases |
| `scripts/analyze_usage.py` | Usage log analysis |
| `scripts/run_chain_eval.py` | Chain evaluation runner |
| `scripts/run_trigger_eval.py` | Trigger evaluation runner |
