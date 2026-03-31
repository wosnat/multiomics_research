# multiomics_research

Research methodology, evaluation framework, and usage logging for multi-omics analysis of marine cyanobacteria and heterotrophic bacteria using an AI-powered knowledge graph.

## Project Goals

1. **Publish academic papers** — both the AI/LLM research framework and biological findings from multi-omics analysis of *Prochlorococcus*, *Alteromonas*, and *Synechococcus*.
2. **Publish research tools** — a Claude Code plugin, MCP server, and Neo4j knowledge graph for community use.

## Multi-Repo Architecture

This project spans three repositories:

| Repo | Role | Key outputs |
|------|------|-------------|
| [multiomics_biocypher_kg](https://github.com/wosnat/multiomics_biocypher_kg) | **KG Builder** — constructs a Neo4j knowledge graph from multi-omics datasets using BioCypher | Neo4j database, Docker image |
| [multiomics_explorer](https://github.com/wosnat/multiomics_explorer) | **KG Tools** — MCP server (18 tools) and Python API for querying the KG | `multiomics-kg-mcp` server, `multiomics_explorer` Python package |
| [multiomics_research](https://github.com/wosnat/multiomics_research) (this repo) | **Research Framework** — skills, evaluations, and logging for LLM-driven analysis | Claude Code plugin, eval suite, usage analytics |

```
multiomics_biocypher_kg/        # Build the KG
  ├── config/                   #   BioCypher schema & adapter configs
  ├── multiomics_kg/            #   Adapters, loaders, schema
  └── docker/                   #   Neo4j Docker deployment

multiomics_explorer/            # Query the KG
  ├── multiomics_explorer/      #   Python API + MCP server (18 tools)
  ├── tests/                    #   Tool correctness tests
  └── evaluation_data/          #   Test fixtures

multiomics_research/            # Use the KG for research (this repo)
  ├── skills/research/          #   Research methodology skill
  ├── evals/                    #   Evaluation cases
  ├── scripts/                  #   Eval runners & analysis
  ├── hooks/                    #   Usage logging
  └── notebooks/                #   Analysis notebooks
```

## Knowledge Graph Contents

The KG integrates transcriptomics, proteomics, and genomic data for marine microbes:

- **15 organisms** — *Prochlorococcus* (7 strains: MED4, MIT9312, MIT9313, NATL2A, AS9601, MIT9301, NATL1A), *Alteromonas macleodii* (3 strains: EZ55, HOT1A3, MIT1002), *Synechococcus* (WH8102, CC9311)
- **76 differential expression experiments** from 21 publications
- **11 treatment types** — coculture, nitrogen/carbon/iron/phosphorus/light/plastic/salt stress, darkness, viral infection, growth state
- **3 omics platforms** — RNA-seq (48 experiments), microarray (24), proteomics (4)

## Research Skill

The core of this repo is a Claude Code research skill ([skills/research/SKILL.md](skills/research/SKILL.md)) that guides systematic multi-omics analysis. It enforces:

- **KG as sole data source** — never rely on LLM training data for gene facts
- **Locus tags over gene names** — prevents paralog conflation (e.g., `katA` maps to 2 distinct catalase genes per Alteromonas strain with opposite expression patterns)
- **Reproducibility** — analyses produce scripts, data files, and methods documentation, not just chat answers
- **Statistical rigor** — no cross-study p-value comparison, proper multiple testing correction, scope qualification

### Skill References

| File | Purpose |
|------|---------|
| [SKILL.md](skills/research/SKILL.md) | Core methodology: 3-phase workflow (orientation → gene work → expression) |
| [research-checklist.md](skills/research/references/research-checklist.md) | Step-by-step protocol with checkboxes |
| [artifacts-guide.md](skills/research/references/artifacts-guide.md) | Standardized directory structure for reproducible analyses |
| [anti-hallucination.md](skills/research/references/anti-hallucination.md) | Concrete LLM failure modes with prevention strategies |

## Evaluation Framework

Three layers measure how well the LLM uses the KG tools:

| Layer | What it measures | Eval cases |
|-------|-----------------|------------|
| **Usage logging** | Real tool usage patterns via hooks | Continuous (JSONL) |
| **Tool chain evals** | Does the LLM pick correct tool sequences? | [6 cases](evals/chain_evals.yaml) |
| **Trigger evals** | Does the research skill activate appropriately? | [20 cases](evals/trigger_evals.json) (10 should trigger, 10 should not) |
| **Research quality evals** | Does the skill improve analysis output? | [5 cases](evals/research_questions.json) |

## Setup

### Prerequisites

- Python >= 3.11
- [uv](https://docs.astral.sh/uv/) package manager
- [jq](https://jqlang.github.io/jq/) (for usage logging hooks)
- Neo4j running at `bolt://localhost:7687` (see [multiomics_biocypher_kg](https://github.com/wosnat/multiomics_biocypher_kg) for Docker setup)
- [multiomics_explorer](https://github.com/wosnat/multiomics_explorer) cloned as a sibling directory

### Install

```bash
# Clone all three repos
git clone https://github.com/wosnat/multiomics_biocypher_kg.git
git clone https://github.com/wosnat/multiomics_explorer.git
git clone https://github.com/wosnat/multiomics_research.git

# Install this repo's dependencies
cd multiomics_research
uv sync --extra analysis
```

### Verify MCP Server

```bash
# Start Neo4j (if not already running)
cd ../multiomics_biocypher_kg
docker compose up -d

# Test MCP server
uv run --directory ../multiomics_explorer multiomics-kg-mcp
```

### Load as Claude Code Plugin

```bash
claude --plugin-dir /path/to/multiomics_research
```

The plugin provides the research skill and usage logging hooks automatically.

## Running Evaluations

```bash
# Analyze logged usage data
python scripts/analyze_usage.py ~/.claude/logs/multiomics-kg-usage.jsonl

# Run tool chain evaluation (requires Neo4j + MCP server)
python scripts/run_chain_eval.py evals/chain_evals.yaml

# Run skill trigger evaluation (requires claude CLI)
python scripts/run_trigger_eval.py evals/trigger_evals.json --skill-path skills/research
```

## Usage Logging

The plugin's PostToolUse hook logs every MCP tool call to `~/.claude/logs/multiomics-kg-usage.jsonl`. Each entry records:

- Tool name and input parameters
- Response metadata (truncation, result counts, errors)
- Session ID and timestamp

Analyze patterns with:

```bash
python scripts/analyze_usage.py ~/.claude/logs/multiomics-kg-usage.jsonl --output report.md
```

## Project Status

- [x] KG built and deployed (Neo4j Docker)
- [x] MCP server with 18 tools
- [x] Research skill with anti-hallucination safeguards
- [x] Evaluation framework scaffolded (chain, trigger, quality)
- [x] Usage logging hooks
- [ ] Baseline eval results
- [ ] First biological analysis (driving example)
- [ ] Paper draft

## License

TBD
