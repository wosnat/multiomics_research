# Pathway Enrichment B1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **REQUIRED SKILL:** Load the `research-methodology` skill BEFORE starting any task. It contains non-negotiable rules for KG research, artifact structure, and notebook discipline.
>
> **Process overrides from CLAUDE.md:**
> - Don't skip subagent reviews for tasks that produce data outputs
> - Notebook-commit gates: commit notebook entry for each step before beginning the next
> - Interactive discovery steps: frozen output + notebook entry, not scripts

**Goal:** Build a reusable `enrich_utils` package for ontology-aware pathway enrichment analysis, validate it with toy and real KG data, then run the full pathway enrichment analysis on Weissberg 2025 MED4 data across all v2 experiments.

**Architecture:** Two-phase pipeline. Phase 0 builds and validates `enrich_utils` (4 modules: extraction, hierarchy, survey, enrichment). Phase 0c is a potential MCP pause point — if hierarchy extraction requires tools that don't exist yet, we stop, file an MCP enhancement request, and wait. Phase 1 runs the interactive analysis using the validated utilities.

**Tech Stack:** Python 3.13, pandas, scipy (fisher_exact), statsmodels (multipletests), matplotlib/seaborn (heatmaps). KG access via `multiomics_explorer` Python API and `run_cypher`. Run scripts with `uv run` from `multiomics_research` root.

**Spec:** `docs/superpowers/specs/2026-04-09-pathway-enrichment-b1-design.md`

**Predecessor data:** `analyses/2026-04-08-1038-n_limitation_signature_v2/data/de_*.csv`

---

## File Structure

### Analysis directory (created in Task 1)

```
analyses/YYYY-MM-DD-HHMM-pathway_enrichment_b1/
├── enrich_utils/           # Reusable methodology package
│   ├── __init__.py
│   ├── enrichment.py       # Fisher's exact, FDR, pathway coverage
│   ├── hierarchy.py        # Roll-up: gene × leaf → gene × level
│   ├── survey.py           # Coverage, term stats, ranking
│   ├── extraction.py       # KG data extraction (only KG-touching module)
│   ├── io.py               # Load/save helpers
│   └── tests/
│       ├── test_enrichment.py
│       ├── test_survey.py
│       └── test_extraction.py
├── scripts/                # Phase 1 pipeline scripts
│   ├── 01_extract_annotations.py
│   ├── 02_survey_landscape.py
│   ├── 03_define_pathways.py
│   ├── 04_run_enrichment.py
│   └── 05_plot_results.py
├── data/                   # Staged KG extracts
├── results/                # Analysis outputs
├── logs/                   # Per-step diagnostic logs
├── exploration/            # Research notebook
├── superpowers/            # Spec, plan, brainstorm-log copies
├── README.md
├── methods.md
├── decisions.md
├── caveats.md
├── gaps_and_friction.md
├── mcp_tool_requirements.md
└── references.md
```

---

## Phase 0: Build and Validate `enrich_utils`

### Task 1: Scaffold analysis directory and package structure [Spec: directory structure]

**Files:**
- Create: `analyses/YYYY-MM-DD-HHMM-pathway_enrichment_b1/` (full tree)
- Create: `enrich_utils/__init__.py`
- Create: All analysis doc stubs

- [ ] **Step 1: Create analysis directory with timestamp**

```bash
ANALYSIS_DIR="analyses/$(date +%Y-%m-%d-%H%M)-pathway_enrichment_b1"
mkdir -p "$ANALYSIS_DIR"/{enrich_utils/tests,scripts,data,results,logs,exploration,superpowers}
echo "Created $ANALYSIS_DIR"
```

- [ ] **Step 2: Create package __init__.py**

```python
# enrich_utils/__init__.py
"""Pathway enrichment utilities for multi-omics KG research.

Four modules:
- extraction: KG data extraction (annotations, hierarchies)
- hierarchy: roll-up from leaf to chosen level (pure DataFrame)
- survey: annotation landscape profiling and ranking (pure DataFrame)
- enrichment: Fisher's exact, FDR, pathway coverage (pure DataFrame)
"""
```

- [ ] **Step 3: Create analysis document stubs**

Create each of these files with a title header and empty sections per the artifacts guide in the research-methodology skill:
- `README.md` — with analysis title, question, links to notebook and spec
- `methods.md` — sections: Research question, Data scope, Gene selection, Statistical methods, Results summary, Limitations
- `decisions.md` — "Decision Log" header
- `caveats.md` — "Caveats for Interpretation" header
- `gaps_and_friction.md` — sections: KG data bugs, KG gaps, MCP friction, Skill/methodology friction, Process retrospective
- `mcp_tool_requirements.md` — "MCP Tool Requirements" header, sections: Queries needed, Proposed tool API
- `references.md` — "References and Citations" header
- `data/DATA_MANIFEST.md` — "Data Manifest" header
- `results/RESULTS_MANIFEST.md` — "Results Manifest" header
- `exploration/YYYY-MM-DD-notebook.md` — "Research Notebook: Pathway Enrichment B1" header

- [ ] **Step 4: Copy spec, plan, and brainstorm-log to superpowers/**

```bash
cp docs/superpowers/specs/2026-04-09-pathway-enrichment-b1-design.md "$ANALYSIS_DIR/superpowers/spec.md"
cp docs/superpowers/specs/2026-04-09-pathway-enrichment-b1-brainstorm-log.md "$ANALYSIS_DIR/superpowers/brainstorm-log.md"
cp docs/superpowers/plans/2026-04-09-pathway-enrichment-b1-plan.md "$ANALYSIS_DIR/superpowers/plan.md"
```

- [ ] **Step 5: Commit scaffold**

```bash
git add "$ANALYSIS_DIR"
git commit -m "scaffold: pathway enrichment B1 analysis directory"
```

---

### Task 2: Implement enrichment.py — Fisher's exact + FDR [Spec: Step 0a]

Core enrichment computation. Pure DataFrame in/out. No KG dependency.

**Files:**
- Create: `enrich_utils/enrichment.py`
- Create: `enrich_utils/tests/test_enrichment.py`

- [ ] **Step 1: Write test_enrichment.py — strong enrichment (up)**

```python
# enrich_utils/tests/test_enrichment.py
"""Toy-data tests for enrichment functions.

Each test uses hand-calculated expected values from the spec worked examples.
"""
import pandas as pd
import pytest
from scipy.stats import fisher_exact


def _make_de_df(genes_up, genes_down, genes_ns, locus_tags_all):
    """Helper: build a DE DataFrame with expression_status column."""
    rows = []
    for lt in locus_tags_all:
        if lt in genes_up:
            rows.append({"locus_tag": lt, "expression_status": "significant_up"})
        elif lt in genes_down:
            rows.append({"locus_tag": lt, "expression_status": "significant_down"})
        else:
            rows.append({"locus_tag": lt, "expression_status": "not_significant"})
    return pd.DataFrame(rows)


def _make_pathway_defs(pathways: dict[str, list[str]]) -> pd.DataFrame:
    """Helper: build pathway definitions DataFrame.

    pathways: {"pathway_name": [locus_tag, ...], ...}
    """
    rows = []
    for name, genes in pathways.items():
        rows.append({
            "pathway_id": name.lower().replace(" ", "_"),
            "pathway_name": name,
            "locus_tags": set(genes),
            "gene_count": len(genes),
        })
    return pd.DataFrame(rows)


class TestFisherEnrichment:
    """Test cases 1-7 from spec step 0a."""

    def test_strong_enrichment_up(self):
        """Case 1: 10-gene pathway, 8 DE up out of 10 in universe, 30 DE up total, 100 genes."""
        from enrich_utils.enrichment import run_enrichment

        all_genes = [f"G{i:03d}" for i in range(100)]
        pathway_genes = all_genes[:10]
        de_up = all_genes[:8] + all_genes[10:32]  # 8 in pathway + 22 outside = 30 total

        de_df = _make_de_df(
            genes_up=set(de_up), genes_down=set(), genes_ns=set(),
            locus_tags_all=all_genes,
        )
        pathway_defs = _make_pathway_defs({"TestPathway": pathway_genes})

        result = run_enrichment(de_df, pathway_defs, set(all_genes), "all_detected_genes")

        row_up = result[
            (result["pathway_id"] == "testpathway") & (result["direction"] == "up")
        ].iloc[0]

        # 2x2: [[8, 22], [2, 68]]
        expected_or, expected_p = fisher_exact([[8, 22], [2, 68]], alternative="greater")
        assert row_up["a"] == 8
        assert row_up["b"] == 22
        assert row_up["c"] == 2
        assert row_up["d"] == 68
        assert abs(row_up["odds_ratio"] - expected_or) < 0.001
        assert abs(row_up["p_value"] - expected_p) < 0.001
        assert row_up["test_type"] == "vs_genome"
        assert row_up["pathway_coverage"] == 1.0  # all 10 in universe

    def test_no_enrichment(self):
        """Case 2: 10-gene pathway, 1 DE up out of 10, 30 DE up total, 100 genes."""
        from enrich_utils.enrichment import run_enrichment

        all_genes = [f"G{i:03d}" for i in range(100)]
        pathway_genes = all_genes[:10]
        de_up = [all_genes[0]] + all_genes[10:39]  # 1 in pathway + 29 outside = 30 total

        de_df = _make_de_df(
            genes_up=set(de_up), genes_down=set(), genes_ns=set(),
            locus_tags_all=all_genes,
        )
        pathway_defs = _make_pathway_defs({"TestPathway": pathway_genes})

        result = run_enrichment(de_df, pathway_defs, set(all_genes), "all_detected_genes")
        row_up = result[
            (result["pathway_id"] == "testpathway") & (result["direction"] == "up")
        ].iloc[0]

        assert row_up["a"] == 1
        assert row_up["p_value"] > 0.05  # not significant

    def test_down_enrichment(self):
        """Case 3: same as case 1 but down direction."""
        from enrich_utils.enrichment import run_enrichment

        all_genes = [f"G{i:03d}" for i in range(100)]
        pathway_genes = all_genes[:10]
        de_down = all_genes[:8] + all_genes[10:32]  # 8 in pathway + 22 outside = 30

        de_df = _make_de_df(
            genes_up=set(), genes_down=set(de_down), genes_ns=set(),
            locus_tags_all=all_genes,
        )
        pathway_defs = _make_pathway_defs({"TestPathway": pathway_genes})

        result = run_enrichment(de_df, pathway_defs, set(all_genes), "all_detected_genes")
        row_down = result[
            (result["pathway_id"] == "testpathway") & (result["direction"] == "down")
        ].iloc[0]

        expected_or, expected_p = fisher_exact([[8, 22], [2, 68]], alternative="greater")
        assert row_down["a"] == 8
        assert abs(row_down["odds_ratio"] - expected_or) < 0.001

    def test_combined_test(self):
        """Case 4: pathway with 5 up + 5 down out of 10 members."""
        from enrich_utils.enrichment import run_enrichment

        all_genes = [f"G{i:03d}" for i in range(100)]
        pathway_genes = all_genes[:10]
        de_up = all_genes[:5] + all_genes[20:35]    # 5 in pathway + 15 outside = 20
        de_down = all_genes[5:10] + all_genes[35:50]  # 5 in pathway + 15 outside = 20

        de_df = _make_de_df(
            genes_up=set(de_up), genes_down=set(de_down), genes_ns=set(),
            locus_tags_all=all_genes,
        )
        pathway_defs = _make_pathway_defs({"TestPathway": pathway_genes})

        result = run_enrichment(de_df, pathway_defs, set(all_genes), "all_detected_genes")

        row_combined = result[
            (result["pathway_id"] == "testpathway") & (result["direction"] == "combined")
        ].iloc[0]
        row_up = result[
            (result["pathway_id"] == "testpathway") & (result["direction"] == "up")
        ].iloc[0]

        # Combined: 10 DE in pathway out of 40 DE total
        assert row_combined["a"] == 10
        assert row_combined["a"] + row_combined["b"] == 40
        # Up: only 5
        assert row_up["a"] == 5

    def test_empty_pathway_in_universe(self):
        """Case 5: pathway has 10 genes but 0 in this experiment's universe."""
        from enrich_utils.enrichment import run_enrichment

        universe_genes = [f"G{i:03d}" for i in range(100, 200)]  # disjoint from pathway
        pathway_genes = [f"G{i:03d}" for i in range(10)]

        de_df = _make_de_df(
            genes_up=set(universe_genes[:30]), genes_down=set(), genes_ns=set(),
            locus_tags_all=universe_genes,
        )
        pathway_defs = _make_pathway_defs({"TestPathway": pathway_genes})

        result = run_enrichment(de_df, pathway_defs, set(universe_genes), "all_detected_genes")
        row_up = result[
            (result["pathway_id"] == "testpathway") & (result["direction"] == "up")
        ].iloc[0]

        assert row_up["pathway_coverage"] == 0.0
        assert pd.isna(row_up["p_value"])  # test skipped

    def test_underpowered_pathway(self):
        """Case 6: pathway has 10 genes but only 3 in universe."""
        from enrich_utils.enrichment import run_enrichment

        all_genes = [f"G{i:03d}" for i in range(100)]
        pathway_genes = all_genes[:3] + [f"X{i:03d}" for i in range(7)]  # 3 in universe, 7 outside

        de_df = _make_de_df(
            genes_up=set(all_genes[:30]), genes_down=set(), genes_ns=set(),
            locus_tags_all=all_genes,
        )
        pathway_defs = _make_pathway_defs({"TestPathway": pathway_genes})

        result = run_enrichment(de_df, pathway_defs, set(all_genes), "all_detected_genes")
        row_up = result[
            (result["pathway_id"] == "testpathway") & (result["direction"] == "up")
        ].iloc[0]

        assert row_up["pathway_coverage"] == 0.3
        assert row_up["underpowered"] == True  # flagged

    def test_all_pathway_genes_de(self):
        """Case 7: all 10 pathway genes are DE up. Maximum enrichment."""
        from enrich_utils.enrichment import run_enrichment

        all_genes = [f"G{i:03d}" for i in range(100)]
        pathway_genes = all_genes[:10]
        de_up = all_genes[:10] + all_genes[20:40]  # all 10 in pathway + 20 outside = 30

        de_df = _make_de_df(
            genes_up=set(de_up), genes_down=set(), genes_ns=set(),
            locus_tags_all=all_genes,
        )
        pathway_defs = _make_pathway_defs({"TestPathway": pathway_genes})

        result = run_enrichment(de_df, pathway_defs, set(all_genes), "all_detected_genes")
        row_up = result[
            (result["pathway_id"] == "testpathway") & (result["direction"] == "up")
        ].iloc[0]

        assert row_up["a"] == 10
        assert row_up["c"] == 0  # none left in pathway that aren't DE
        assert row_up["p_value"] < 0.001


class TestFDRCorrection:
    """Case 8: BH correction across 20 pathways."""

    def test_bh_correction(self):
        from enrich_utils.enrichment import run_enrichment

        all_genes = [f"G{i:03d}" for i in range(500)]
        de_up = set(all_genes[:100])  # 100 DE up out of 500

        # 20 pathways: 2 strongly enriched, 18 not
        pathways = {}
        for i in range(20):
            if i < 2:
                # Enriched: 15 of 20 genes are DE up
                genes = all_genes[i*20:(i*20)+15] + [f"NEG{i}_{j}" for j in range(5)]
                # Hack: put neg genes inside universe
                genes = all_genes[i*20:(i*20)+20]  # first 15 are DE, last 5 not
            else:
                # Not enriched: 4 of 20 genes are DE up (= 100/500 * 20 = 4 expected)
                genes = all_genes[200 + i*20:200 + (i*20)+20]
            pathways[f"Pathway_{i:02d}"] = genes

        de_df = _make_de_df(genes_up=de_up, genes_down=set(), genes_ns=set(), locus_tags_all=all_genes)
        pathway_defs = _make_pathway_defs(pathways)

        result = run_enrichment(de_df, pathway_defs, set(all_genes), "all_detected_genes")
        up_results = result[result["direction"] == "up"].copy()

        # padj should exist and be >= p_value for all rows
        assert (up_results["padj"] >= up_results["p_value"]).all()
        # At least the 2 enriched pathways should have low padj
        enriched = up_results[up_results["pathway_id"].isin(["pathway_00", "pathway_01"])]
        assert (enriched["padj"] < 0.05).all()


class TestTestType:
    """Case 9: table_scope maps to test_type."""

    @pytest.mark.parametrize("scope,expected_type", [
        ("all_detected_genes", "vs_genome"),
        ("filtered_subset", "vs_filtered_genome"),
        ("significant_any_timepoint", "vs_all_responsive"),
        ("significant_only", "descriptive_only"),
    ])
    def test_table_scope_mapping(self, scope, expected_type):
        from enrich_utils.enrichment import run_enrichment

        all_genes = [f"G{i:03d}" for i in range(50)]
        de_df = _make_de_df(genes_up=set(all_genes[:10]), genes_down=set(), genes_ns=set(), locus_tags_all=all_genes)
        pathway_defs = _make_pathway_defs({"P1": all_genes[:5]})

        result = run_enrichment(de_df, pathway_defs, set(all_genes), scope)
        assert (result["test_type"] == expected_type).all()
```

- [ ] **Step 2: Run tests — verify they fail**

```bash
cd analyses/YYYY-MM-DD-HHMM-pathway_enrichment_b1
uv run pytest enrich_utils/tests/test_enrichment.py -v
```
Expected: ImportError — `enrich_utils.enrichment` doesn't exist yet.

- [ ] **Step 3: Implement enrichment.py**

```python
# enrich_utils/enrichment.py
"""Pathway enrichment: Fisher's exact test with FDR correction.

Pure DataFrame-in, DataFrame-out. Never calls the KG.
"""
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


_TABLE_SCOPE_TO_TEST_TYPE = {
    "all_detected_genes": "vs_genome",
    "filtered_subset": "vs_filtered_genome",
    "significant_any_timepoint": "vs_all_responsive",
    "significant_only": "descriptive_only",
}


def _fisher_one_pathway(
    de_genes: set[str],
    pathway_genes_in_universe: set[str],
    n_universe: int,
) -> dict:
    """Run one-sided Fisher's exact for a single pathway × direction.

    Returns dict with a, b, c, d, odds_ratio, p_value, fold_enrichment, expected.
    Returns NaN for p_value/odds_ratio if pathway has 0 genes in universe.
    """
    a = len(de_genes & pathway_genes_in_universe)
    c = len(pathway_genes_in_universe - de_genes)
    b = len(de_genes) - a
    d = n_universe - a - b - c

    if len(pathway_genes_in_universe) == 0:
        return {
            "a": 0, "b": len(de_genes), "c": 0, "d": n_universe - len(de_genes),
            "odds_ratio": np.nan, "p_value": np.nan,
            "fold_enrichment": np.nan, "expected": np.nan,
        }

    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
    expected = len(de_genes) * len(pathway_genes_in_universe) / n_universe
    fold_enrichment = a / expected if expected > 0 else np.nan

    return {
        "a": a, "b": b, "c": c, "d": d,
        "odds_ratio": odds_ratio, "p_value": p_value,
        "fold_enrichment": fold_enrichment, "expected": expected,
    }


def run_enrichment(
    de_df: pd.DataFrame,
    pathway_defs: pd.DataFrame,
    gene_universe: set[str],
    table_scope: str,
) -> pd.DataFrame:
    """Run Fisher's exact per pathway × direction (up, down, combined).

    Parameters
    ----------
    de_df : DataFrame with columns: locus_tag, expression_status
        expression_status in {significant_up, significant_down, not_significant}
    pathway_defs : DataFrame with columns: pathway_id, pathway_name, locus_tags (set), gene_count
    gene_universe : set of locus tags in this experiment × timepoint
    table_scope : str, maps to test_type column

    Returns
    -------
    DataFrame with one row per pathway × direction.
    """
    test_type = _TABLE_SCOPE_TO_TEST_TYPE.get(table_scope, "unknown")
    n_universe = len(gene_universe)

    genes_up = set(de_df.loc[de_df["expression_status"] == "significant_up", "locus_tag"])
    genes_down = set(de_df.loc[de_df["expression_status"] == "significant_down", "locus_tag"])
    genes_de = genes_up | genes_down

    rows = []
    for _, pw in pathway_defs.iterrows():
        pw_genes_in_universe = pw["locus_tags"] & gene_universe
        pw_coverage = len(pw_genes_in_universe) / pw["gene_count"] if pw["gene_count"] > 0 else 0.0
        underpowered = pw_coverage < 0.50

        for direction, de_set in [("up", genes_up), ("down", genes_down), ("combined", genes_de)]:
            stats = _fisher_one_pathway(de_set, pw_genes_in_universe, n_universe)
            rows.append({
                "pathway_id": pw["pathway_id"],
                "pathway_name": pw["pathway_name"],
                "direction": direction,
                **stats,
                "test_type": test_type,
                "pathway_coverage": pw_coverage,
                "n_pathway_genes_in_universe": len(pw_genes_in_universe),
                "n_pathway_genes_total": pw["gene_count"],
                "underpowered": underpowered,
            })

    result = pd.DataFrame(rows)

    # BH correction per direction group
    if len(result) > 0:
        padj_col = np.full(len(result), np.nan)
        for direction in ["up", "down", "combined"]:
            mask = (result["direction"] == direction) & result["p_value"].notna()
            if mask.sum() > 0:
                _, padj, _, _ = multipletests(
                    result.loc[mask, "p_value"], method="fdr_bh"
                )
                padj_col[mask.values.nonzero()[0]] = padj
        result["padj"] = padj_col

    return result


def run_enrichment_all_timepoints(
    de_df: pd.DataFrame,
    pathway_defs: pd.DataFrame,
    table_scope: str,
) -> pd.DataFrame:
    """Run enrichment per timepoint, concatenate results.

    de_df must have a 'timepoint' column. Single-timepoint experiments
    have one unique value (or NaN treated as a single group).
    """
    if "timepoint" not in de_df.columns or de_df["timepoint"].isna().all():
        gene_universe = set(de_df["locus_tag"])
        result = run_enrichment(de_df, pathway_defs, gene_universe, table_scope)
        result["timepoint"] = None
        return result

    frames = []
    for tp, tp_df in de_df.groupby("timepoint"):
        gene_universe = set(tp_df["locus_tag"])
        result = run_enrichment(tp_df, pathway_defs, gene_universe, table_scope)
        result["timepoint"] = tp
        frames.append(result)

    return pd.concat(frames, ignore_index=True)
```

- [ ] **Step 4: Run tests — verify they pass**

```bash
uv run pytest enrich_utils/tests/test_enrichment.py -v
```
Expected: All pass.

- [ ] **Step 5: Commit**

```bash
git add enrich_utils/enrichment.py enrich_utils/tests/test_enrichment.py
git commit -m "feat: add enrichment.py with Fisher's exact + FDR [Step 0a]"
```

---

### Task 3: Implement hierarchy.py — roll-up logic [Spec: Step 0b]

Pure DataFrame operations: propagate gene annotations up a hierarchy DAG.

**Files:**
- Create: `enrich_utils/hierarchy.py`
- Test in: `enrich_utils/tests/test_survey.py` (roll-up tests)

- [ ] **Step 1: Write roll-up tests in test_survey.py**

```python
# enrich_utils/tests/test_survey.py
"""Toy-data tests for hierarchy roll-up, survey, and pathway definitions.

Synthetic 3-level hierarchy with DAG structure and multi-annotation genes.
"""
import pandas as pd
import pytest


def _make_hierarchy():
    """Synthetic 3-level hierarchy (DAG, not tree).

    Level 0 (roots): R1, R2
    Level 1 (mid):   M1→R1, M2→R1, M3→R2, M4→R2
    Level 2 (leaves): L1→M1, L2→M1, L3→M2, L4→M2, L5→M3,
                      L6→M3, L6→M4 (DAG: L6 has two parents)
    """
    edges = pd.DataFrame([
        {"child_id": "M1", "parent_id": "R1", "child_level": 1, "parent_level": 0},
        {"child_id": "M2", "parent_id": "R1", "child_level": 1, "parent_level": 0},
        {"child_id": "M3", "parent_id": "R2", "child_level": 1, "parent_level": 0},
        {"child_id": "M4", "parent_id": "R2", "child_level": 1, "parent_level": 0},
        {"child_id": "L1", "parent_id": "M1", "child_level": 2, "parent_level": 1},
        {"child_id": "L2", "parent_id": "M1", "child_level": 2, "parent_level": 1},
        {"child_id": "L3", "parent_id": "M2", "child_level": 2, "parent_level": 1},
        {"child_id": "L4", "parent_id": "M2", "child_level": 2, "parent_level": 1},
        {"child_id": "L5", "parent_id": "M3", "child_level": 2, "parent_level": 1},
        {"child_id": "L6", "parent_id": "M3", "child_level": 2, "parent_level": 1},
        {"child_id": "L6", "parent_id": "M4", "child_level": 2, "parent_level": 1},  # DAG edge
    ])
    term_names = {
        "R1": "Root One", "R2": "Root Two",
        "M1": "Mid One", "M2": "Mid Two", "M3": "Mid Three", "M4": "Mid Four",
        "L1": "Leaf One", "L2": "Leaf Two", "L3": "Leaf Three",
        "L4": "Leaf Four", "L5": "Leaf Five", "L6": "Leaf Six",
    }
    return edges, term_names


def _make_annotations():
    """Genes annotated at various levels.

    G01 → L1 (leaf)
    G02 → L1, L3 (two leaves, different branches → both R1 children)
    G03 → M2 (mid-level directly, not via leaf)
    G04 → L5 (leaf under R2)
    G05 → L6 (leaf under M3 AND M4 — DAG convergence → one count in R2)
    G06 → R1 (root directly)
    G07 → (no annotation)
    G08 → L1, M1 (both leaf and its parent — dedup)
    """
    return pd.DataFrame([
        {"locus_tag": "G01", "term_id": "L1"},
        {"locus_tag": "G02", "term_id": "L1"},
        {"locus_tag": "G02", "term_id": "L3"},
        {"locus_tag": "G03", "term_id": "M2"},
        {"locus_tag": "G04", "term_id": "L5"},
        {"locus_tag": "G05", "term_id": "L6"},
        {"locus_tag": "G06", "term_id": "R1"},
        {"locus_tag": "G08", "term_id": "L1"},
        {"locus_tag": "G08", "term_id": "M1"},
    ])


class TestRollUp:

    def test_leaf_rolls_up_to_parent_and_root(self):
        """G01 annotated to L1 → appears in M1 and R1."""
        from enrich_utils.hierarchy import roll_up_to_level
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()

        rolled_1 = roll_up_to_level(annotations, hierarchy, target_level=1)
        rolled_0 = roll_up_to_level(annotations, hierarchy, target_level=0)

        g01_mid = set(rolled_1[rolled_1["locus_tag"] == "G01"]["term_id"])
        g01_root = set(rolled_0[rolled_0["locus_tag"] == "G01"]["term_id"])
        assert "M1" in g01_mid
        assert "R1" in g01_root

    def test_mid_level_annotation_propagates_up_not_down(self):
        """G03 annotated to M2 → appears in R1 at level 0, but NOT in L3/L4 at level 2."""
        from enrich_utils.hierarchy import roll_up_to_level
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()

        rolled_0 = roll_up_to_level(annotations, hierarchy, target_level=0)
        rolled_2 = roll_up_to_level(annotations, hierarchy, target_level=2)

        g03_roots = set(rolled_0[rolled_0["locus_tag"] == "G03"]["term_id"])
        g03_leaves = set(rolled_2[rolled_2["locus_tag"] == "G03"]["term_id"])
        assert "R1" in g03_roots
        assert "L3" not in g03_leaves  # no downward propagation
        assert "L4" not in g03_leaves

    def test_multi_branch_annotation(self):
        """G02 in L1 (→M1→R1) and L3 (→M2→R1) → appears in both M1 and M2 at level 1."""
        from enrich_utils.hierarchy import roll_up_to_level
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()

        rolled_1 = roll_up_to_level(annotations, hierarchy, target_level=1)
        g02_mid = set(rolled_1[rolled_1["locus_tag"] == "G02"]["term_id"])
        assert "M1" in g02_mid
        assert "M2" in g02_mid

    def test_dag_convergence_dedup(self):
        """G05 in L6 (→M3, →M4, both →R2) → counted once in R2."""
        from enrich_utils.hierarchy import roll_up_to_level
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()

        rolled_0 = roll_up_to_level(annotations, hierarchy, target_level=0)
        g05_roots = rolled_0[rolled_0["locus_tag"] == "G05"]
        r2_count = len(g05_roots[g05_roots["term_id"] == "R2"])
        assert r2_count == 1  # deduplicated

    def test_leaf_and_parent_dedup(self):
        """G08 in L1 AND M1 → counted once in M1 at level 1."""
        from enrich_utils.hierarchy import roll_up_to_level
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()

        rolled_1 = roll_up_to_level(annotations, hierarchy, target_level=1)
        g08_mid = rolled_1[rolled_1["locus_tag"] == "G08"]
        m1_count = len(g08_mid[g08_mid["term_id"] == "M1"])
        assert m1_count == 1

    def test_root_annotation_stays_at_root(self):
        """G06 annotated directly to R1 → appears at level 0, not at level 1 or 2."""
        from enrich_utils.hierarchy import roll_up_to_level
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()

        rolled_0 = roll_up_to_level(annotations, hierarchy, target_level=0)
        rolled_1 = roll_up_to_level(annotations, hierarchy, target_level=1)

        g06_roots = set(rolled_0[rolled_0["locus_tag"] == "G06"]["term_id"])
        g06_mid = set(rolled_1[rolled_1["locus_tag"] == "G06"]["term_id"])
        assert "R1" in g06_roots
        assert len(g06_mid) == 0  # no downward propagation

    def test_unannotated_gene_absent(self):
        """G07 has no annotation → absent from all levels."""
        from enrich_utils.hierarchy import roll_up_to_level
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()

        for level in [0, 1, 2]:
            rolled = roll_up_to_level(annotations, hierarchy, target_level=level)
            assert "G07" not in rolled["locus_tag"].values
```

- [ ] **Step 2: Run tests — verify they fail**

```bash
uv run pytest enrich_utils/tests/test_survey.py::TestRollUp -v
```
Expected: ImportError.

- [ ] **Step 3: Implement hierarchy.py**

```python
# enrich_utils/hierarchy.py
"""Hierarchy roll-up: propagate gene annotations upward through a DAG.

Pure DataFrame operations. Takes annotations (gene × term) and hierarchy
(child → parent edges with levels), produces gene × term at a target level.

Annotations propagate UP only (leaf → parent → root). Never down.
Deduplication: a gene appears once per term at the target level,
even if it reaches that term through multiple paths (DAG convergence).
"""
import pandas as pd


def _build_ancestor_map(
    hierarchy_df: pd.DataFrame,
    target_level: int,
) -> dict[str, set[str]]:
    """Build a map: term_id → set of ancestor term_ids at target_level.

    Traverses the hierarchy upward from each term. A term at the target level
    maps to itself. A term below target level maps to its ancestors at that level.
    A term above target level maps to empty set (no downward propagation).
    """
    # Build child→parents adjacency
    child_to_parents: dict[str, set[str]] = {}
    term_levels: dict[str, int] = {}

    for _, row in hierarchy_df.iterrows():
        child_id = row["child_id"]
        parent_id = row["parent_id"]
        child_to_parents.setdefault(child_id, set()).add(parent_id)
        term_levels[child_id] = row["child_level"]
        term_levels[parent_id] = row["parent_level"]

    # For each term, find all ancestors at target_level via BFS upward
    ancestor_cache: dict[str, set[str]] = {}

    def _ancestors_at_level(term_id: str) -> set[str]:
        if term_id in ancestor_cache:
            return ancestor_cache[term_id]

        level = term_levels.get(term_id)
        if level is None:
            ancestor_cache[term_id] = set()
            return set()

        if level == target_level:
            ancestor_cache[term_id] = {term_id}
            return {term_id}

        if level < target_level:
            # Term is above target level — no downward propagation
            ancestor_cache[term_id] = set()
            return set()

        # level > target_level — propagate up
        result = set()
        for parent in child_to_parents.get(term_id, set()):
            result |= _ancestors_at_level(parent)

        ancestor_cache[term_id] = result
        return result

    # Precompute for all known terms
    for term_id in list(term_levels.keys()):
        _ancestors_at_level(term_id)

    return ancestor_cache


def roll_up_to_level(
    annotations_df: pd.DataFrame,
    hierarchy_df: pd.DataFrame,
    target_level: int,
) -> pd.DataFrame:
    """Roll up gene annotations to a target hierarchy level.

    Parameters
    ----------
    annotations_df : DataFrame with columns: locus_tag, term_id
        Gene × term annotations (can be at any level).
    hierarchy_df : DataFrame with columns: child_id, parent_id, child_level, parent_level
    target_level : int, the hierarchy level to roll up to (0 = root)

    Returns
    -------
    DataFrame with columns: locus_tag, term_id
        Deduplicated: each gene appears at most once per target-level term.
    """
    ancestor_map = _build_ancestor_map(hierarchy_df, target_level)

    rows = []
    for _, row in annotations_df.iterrows():
        locus_tag = row["locus_tag"]
        term_id = row["term_id"]
        ancestors = ancestor_map.get(term_id, set())
        for anc in ancestors:
            rows.append({"locus_tag": locus_tag, "term_id": anc})

    if not rows:
        return pd.DataFrame(columns=["locus_tag", "term_id"])

    result = pd.DataFrame(rows).drop_duplicates()
    return result
```

- [ ] **Step 4: Run tests — verify they pass**

```bash
uv run pytest enrich_utils/tests/test_survey.py::TestRollUp -v
```
Expected: All pass.

- [ ] **Step 5: Commit**

```bash
git add enrich_utils/hierarchy.py enrich_utils/tests/test_survey.py
git commit -m "feat: add hierarchy.py with DAG-aware roll-up [Step 0b]"
```

---

### Task 4: Implement survey.py — profiling and pathway definitions [Spec: Step 0b]

Coverage stats, term-size distributions, ontology ranking, and pathway definition building.

**Files:**
- Create: `enrich_utils/survey.py`
- Add to: `enrich_utils/tests/test_survey.py`

- [ ] **Step 1: Write survey and pathway definition tests**

Add to `test_survey.py`:

```python
class TestSurveyOntology:

    def test_coverage(self):
        """Coverage = fraction of gene_universe with ≥1 annotation."""
        from enrich_utils.survey import survey_ontology
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()
        gene_universe = {"G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08"}

        profile = survey_ontology(annotations, hierarchy, gene_universe)
        # G07 has no annotation → 7/8 = 0.875
        assert profile["coverage"] == pytest.approx(7 / 8)

    def test_term_count_per_level(self):
        """Count terms with ≥1 gene at each hierarchy level."""
        from enrich_utils.survey import survey_ontology
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()
        gene_universe = {"G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08"}

        profile = survey_ontology(annotations, hierarchy, gene_universe)
        # Level 0: R1 (G01,G02,G03,G06,G08), R2 (G04,G05) → 2 terms
        assert profile["per_level"][0]["n_terms_with_genes"] == 2

    def test_term_size_distribution(self):
        """Term-size distribution at level 0."""
        from enrich_utils.survey import survey_ontology
        hierarchy, _ = _make_hierarchy()
        annotations = _make_annotations()
        gene_universe = {"G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08"}

        profile = survey_ontology(annotations, hierarchy, gene_universe)
        level_0 = profile["per_level"][0]
        # R1 has 5 genes, R2 has 2 genes
        assert level_0["min_genes"] == 2
        assert level_0["max_genes"] == 5


class TestBuildPathwayDefinitions:

    def test_build_at_level_1(self):
        """Pathway definitions at level 1 from toy hierarchy."""
        from enrich_utils.survey import build_pathway_definitions
        hierarchy, term_names = _make_hierarchy()
        annotations = _make_annotations()

        pw_defs = build_pathway_definitions(annotations, hierarchy, level=1, min_genes=1, term_names=term_names)

        # M1 should have: G01 (via L1), G02 (via L1), G08 (via L1 and M1 direct)
        m1 = pw_defs[pw_defs["pathway_id"] == "M1"].iloc[0]
        assert "G01" in m1["locus_tags"]
        assert "G02" in m1["locus_tags"]
        assert "G08" in m1["locus_tags"]

    def test_min_genes_filter(self):
        """Pathways with fewer than min_genes are excluded."""
        from enrich_utils.survey import build_pathway_definitions
        hierarchy, term_names = _make_hierarchy()
        annotations = _make_annotations()

        pw_defs = build_pathway_definitions(annotations, hierarchy, level=1, min_genes=3, term_names=term_names)

        # Only M1 has ≥3 genes (G01, G02, G08)
        assert len(pw_defs) >= 1
        assert "M1" in pw_defs["pathway_id"].values


class TestScopePathways:

    def test_coverage_per_universe(self):
        """scope_pathways_to_universe computes correct coverage fractions."""
        from enrich_utils.survey import build_pathway_definitions, scope_pathways_to_universe
        hierarchy, term_names = _make_hierarchy()
        annotations = _make_annotations()

        pw_defs = build_pathway_definitions(annotations, hierarchy, level=0, min_genes=1, term_names=term_names)
        # R1 has genes: G01, G02, G03, G06, G08
        universe = {"G01", "G02", "G03"}  # 3 of 5 in R1's full set

        scoped = scope_pathways_to_universe(pw_defs, universe)
        r1 = scoped[scoped["pathway_id"] == "R1"].iloc[0]
        assert r1["n_in_universe"] == 3
        assert r1["coverage"] == pytest.approx(3 / 5)
```

- [ ] **Step 2: Run tests — verify they fail**

```bash
uv run pytest enrich_utils/tests/test_survey.py -v -k "not TestRollUp"
```

- [ ] **Step 3: Implement survey.py**

```python
# enrich_utils/survey.py
"""Annotation landscape survey: coverage, term stats, ranking, pathway definitions.

Pure DataFrame operations. Uses hierarchy.roll_up_to_level for per-level analysis.
"""
import pandas as pd
import numpy as np
from enrich_utils.hierarchy import roll_up_to_level


def survey_ontology(
    annotations_df: pd.DataFrame,
    hierarchy_df: pd.DataFrame,
    gene_universe: set[str],
) -> dict:
    """Profile an ontology's suitability for enrichment testing.

    Returns dict with:
    - coverage: fraction of gene_universe with ≥1 annotation
    - n_annotated: count of annotated genes
    - n_unannotated: count of unannotated genes
    - per_level: list of dicts, one per hierarchy level, with:
      - level, n_terms_with_genes, min_genes, q1_genes, median_genes, q3_genes, max_genes
    """
    annotated_genes = set(annotations_df["locus_tag"]) & gene_universe
    coverage = len(annotated_genes) / len(gene_universe) if gene_universe else 0.0

    # Determine hierarchy levels
    if len(hierarchy_df) == 0:
        # Flat ontology — leaf annotations are the only level
        filtered = annotations_df[annotations_df["locus_tag"].isin(gene_universe)]
        term_sizes = filtered.groupby("term_id")["locus_tag"].nunique()
        per_level = [{
            "level": 0,
            "n_terms_with_genes": len(term_sizes),
            "min_genes": int(term_sizes.min()) if len(term_sizes) > 0 else 0,
            "q1_genes": float(term_sizes.quantile(0.25)) if len(term_sizes) > 0 else 0,
            "median_genes": float(term_sizes.median()) if len(term_sizes) > 0 else 0,
            "q3_genes": float(term_sizes.quantile(0.75)) if len(term_sizes) > 0 else 0,
            "max_genes": int(term_sizes.max()) if len(term_sizes) > 0 else 0,
        }]
    else:
        all_levels = sorted(set(hierarchy_df["parent_level"]) | set(hierarchy_df["child_level"]))
        per_level = []
        for level in all_levels:
            rolled = roll_up_to_level(annotations_df, hierarchy_df, target_level=level)
            filtered = rolled[rolled["locus_tag"].isin(gene_universe)]
            if len(filtered) == 0:
                continue
            term_sizes = filtered.groupby("term_id")["locus_tag"].nunique()
            per_level.append({
                "level": level,
                "n_terms_with_genes": len(term_sizes),
                "min_genes": int(term_sizes.min()),
                "q1_genes": float(term_sizes.quantile(0.25)),
                "median_genes": float(term_sizes.median()),
                "q3_genes": float(term_sizes.quantile(0.75)),
                "max_genes": int(term_sizes.max()),
            })

    return {
        "coverage": coverage,
        "n_annotated": len(annotated_genes),
        "n_unannotated": len(gene_universe) - len(annotated_genes),
        "per_level": per_level,
    }


def rank_ontologies(profiles: dict[str, dict]) -> pd.DataFrame:
    """Rank ontology profiles by suitability for enrichment.

    profiles: {ontology_name: survey_ontology result, ...}

    Ranking criteria (simple weighted score):
    - Coverage (higher is better)
    - Has a level with median term size 5-50 (sweet spot for power + specificity)
    """
    rows = []
    for name, profile in profiles.items():
        best_level = None
        best_median = 0
        for lp in profile.get("per_level", []):
            median = lp.get("median_genes", 0)
            if 5 <= median <= 50 and median > best_median:
                best_median = median
                best_level = lp["level"]

        rows.append({
            "ontology": name,
            "coverage": profile["coverage"],
            "best_level": best_level,
            "best_level_median_genes": best_median,
            "has_sweet_spot": best_level is not None,
        })

    df = pd.DataFrame(rows)
    # Sort: sweet spot first, then by coverage descending
    df = df.sort_values(
        ["has_sweet_spot", "coverage"], ascending=[False, False]
    ).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)
    return df


def build_pathway_definitions(
    annotations_df: pd.DataFrame,
    hierarchy_df: pd.DataFrame,
    level: int,
    min_genes: int,
    term_names: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Build pathway gene sets at a chosen hierarchy level.

    Returns DataFrame with: pathway_id, pathway_name, locus_tags (set), gene_count
    """
    rolled = roll_up_to_level(annotations_df, hierarchy_df, target_level=level)

    if len(rolled) == 0:
        return pd.DataFrame(columns=["pathway_id", "pathway_name", "locus_tags", "gene_count"])

    grouped = rolled.groupby("term_id")["locus_tag"].apply(set).reset_index()
    grouped.columns = ["pathway_id", "locus_tags"]
    grouped["gene_count"] = grouped["locus_tags"].apply(len)
    grouped = grouped[grouped["gene_count"] >= min_genes].reset_index(drop=True)

    if term_names:
        grouped["pathway_name"] = grouped["pathway_id"].map(term_names).fillna(grouped["pathway_id"])
    else:
        grouped["pathway_name"] = grouped["pathway_id"]

    return grouped[["pathway_id", "pathway_name", "locus_tags", "gene_count"]]


def scope_pathways_to_universe(
    pathway_defs: pd.DataFrame,
    gene_universe: set[str],
) -> pd.DataFrame:
    """Scope pathway definitions to a specific experiment's gene universe.

    Returns copy with added columns: n_in_universe, coverage, locus_tags_in_universe
    """
    result = pathway_defs.copy()
    result["locus_tags_in_universe"] = result["locus_tags"].apply(lambda s: s & gene_universe)
    result["n_in_universe"] = result["locus_tags_in_universe"].apply(len)
    result["coverage"] = result["n_in_universe"] / result["gene_count"]
    return result
```

- [ ] **Step 4: Run all tests**

```bash
uv run pytest enrich_utils/tests/test_survey.py -v
```
Expected: All pass (both TestRollUp and new tests).

- [ ] **Step 5: Commit**

```bash
git add enrich_utils/survey.py enrich_utils/tests/test_survey.py
git commit -m "feat: add survey.py with profiling and pathway definitions [Step 0b]"
```

---

### Task 5: Implement extraction.py — KG data extraction [Spec: Step 0c]

This is the KG-touching module. Calls Python API and `run_cypher`.

**Files:**
- Create: `enrich_utils/extraction.py`
- Create: `enrich_utils/tests/test_extraction.py`

- [ ] **Step 1: Implement extraction.py**

```python
# enrich_utils/extraction.py
"""KG data extraction: annotations and hierarchies.

Only module that calls the KG. All others are pure DataFrame operations.

Uses the multiomics_explorer Python API for annotation extraction and
run_cypher for hierarchy extraction. Edge names from ONTOLOGY_CONFIG in
multiomics_explorer/kg/queries_lib.py.
"""
import pandas as pd


# Edge names per ontology (from ONTOLOGY_CONFIG in queries_lib.py)
ONTOLOGY_EDGES = {
    "go_bp": {
        "gene_rel": "Gene_involved_in_biological_process",
        "hierarchy_rels": ["Biological_process_is_a_biological_process",
                           "Biological_process_part_of_biological_process"],
        "node_label": "BiologicalProcess",
    },
    "go_mf": {
        "gene_rel": "Gene_enables_molecular_function",
        "hierarchy_rels": ["Molecular_function_is_a_molecular_function",
                           "Molecular_function_part_of_molecular_function"],
        "node_label": "MolecularFunction",
    },
    "go_cc": {
        "gene_rel": "Gene_located_in_cellular_component",
        "hierarchy_rels": ["Cellular_component_is_a_cellular_component",
                           "Cellular_component_part_of_cellular_component"],
        "node_label": "CellularComponent",
    },
    "kegg": {
        "gene_rel": "Gene_has_kegg_ko",
        "hierarchy_rels": ["Kegg_term_is_a_kegg_term"],
        "node_label": "KeggTerm",
        "gene_connects_to_level": "ko",
    },
    "cyanorak_role": {
        "gene_rel": "Gene_has_cyanorak_role",
        "hierarchy_rels": ["Cyanorak_role_is_a_cyanorak_role"],
        "node_label": "CyanorakRole",
    },
    "tigr_role": {
        "gene_rel": "Gene_has_tigr_role",
        "hierarchy_rels": [],
        "node_label": "TigrRole",
    },
    "cog_category": {
        "gene_rel": "Gene_in_cog_category",
        "hierarchy_rels": [],
        "node_label": "CogFunctionalCategory",
    },
    "ec": {
        "gene_rel": "Gene_catalyzes_ec_number",
        "hierarchy_rels": ["Ec_number_is_a_ec_number"],
        "node_label": "EcNumber",
    },
    "pfam": {
        "gene_rel": "Gene_has_pfam",
        "hierarchy_rels": ["Pfam_in_pfam_clan"],
        "node_label": "Pfam",
    },
}


def extract_annotations(
    locus_tags: list[str],
    ontology: str,
    conn=None,
) -> pd.DataFrame:
    """Extract gene × term leaf annotations from the KG.

    Uses the Python API gene_ontology_terms(). Returns DataFrame with
    columns: locus_tag, term_id, term_name.
    """
    from multiomics_explorer import gene_ontology_terms

    result = gene_ontology_terms(
        locus_tags=locus_tags,
        ontology=ontology,
        limit=None,
        conn=conn,
    )

    if not result.get("results"):
        return pd.DataFrame(columns=["locus_tag", "term_id", "term_name"])

    df = pd.DataFrame(result["results"])
    return df[["locus_tag", "term_id", "term_name"]]


def extract_hierarchy(
    ontology: str,
    max_depth: int = 4,
    conn=None,
) -> pd.DataFrame:
    """Extract hierarchy edges from the KG via run_cypher.

    Returns DataFrame with columns: child_id, parent_id, child_level, parent_level.
    Capped at max_depth levels from root.

    For flat ontologies (no hierarchy_rels), returns empty DataFrame.
    """
    from multiomics_explorer import run_cypher

    config = ONTOLOGY_EDGES.get(ontology)
    if not config:
        raise ValueError(f"Unknown ontology: {ontology}. Valid: {sorted(ONTOLOGY_EDGES)}")

    hierarchy_rels = config["hierarchy_rels"]
    if not hierarchy_rels:
        return pd.DataFrame(columns=["child_id", "parent_id", "child_level", "parent_level"])

    node_label = config["node_label"]
    rel_pattern = "|".join(hierarchy_rels)

    # Extract all parent-child edges within max_depth from roots
    # Step 1: find roots (nodes with no outgoing hierarchy edge)
    # Step 2: traverse down to max_depth, recording edges with levels
    query = f"""
    MATCH (root:{node_label})
    WHERE NOT (root)-[:{rel_pattern}]->()
    WITH root, 0 AS root_level
    MATCH path = (child:{node_label})-[:{rel_pattern}*1..{max_depth}]->(root)
    WITH child, root,
         length(path) AS depth,
         [n IN nodes(path) | n.id] AS path_nodes
    UNWIND range(0, size(path_nodes)-2) AS i
    WITH path_nodes[i] AS child_id, path_nodes[i+1] AS parent_id,
         size(path_nodes)-1 - i AS child_depth, size(path_nodes)-1 - (i+1) AS parent_depth
    RETURN DISTINCT child_id, parent_id, child_depth AS child_level, parent_depth AS parent_level
    ORDER BY parent_level, child_level, parent_id, child_id
    """

    result = run_cypher(query=query, limit=10000, conn=conn)

    if not result.get("results"):
        return pd.DataFrame(columns=["child_id", "parent_id", "child_level", "parent_level"])

    return pd.DataFrame(result["results"])
```

**Important:** The Cypher query above is a starting point. Per the spec's Cypher verification rule, it **must** be tested interactively before trusting. Step 0c's interactive validation (Task 6) serves this purpose.

- [ ] **Step 2: Write test_extraction.py (marker gene integration tests)**

```python
# enrich_utils/tests/test_extraction.py
"""Integration tests against real KG data with marker genes.

These tests require a running Neo4j + multiomics_explorer.
Run with: uv run pytest enrich_utils/tests/test_extraction.py -v

Marker genes from v2:
- PMM0920 (glnA) — canonical N-limitation, multi-annotated
- PMM0370 (cynA) — N-scavenging transporter
- PMM0550 (rbcL) — photosynthesis
- PMM0030 — unnamed, may have sparse annotations
"""
import pytest
import pandas as pd

MARKER_GENES = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]


@pytest.fixture(scope="module")
def kg_conn():
    """Get a KG connection, skip if unavailable."""
    try:
        from multiomics_explorer.api.connection import get_connection
        conn = get_connection()
        return conn
    except Exception:
        pytest.skip("KG connection not available")


class TestExtractAnnotations:

    def test_cyanorak_annotations(self, kg_conn):
        """Marker genes should have CyanoRak annotations."""
        from enrich_utils.extraction import extract_annotations

        df = extract_annotations(MARKER_GENES, "cyanorak_role", conn=kg_conn)
        assert len(df) > 0
        assert "locus_tag" in df.columns
        assert "term_id" in df.columns
        # glnA should have at least one CyanoRak annotation
        assert "PMM0920" in df["locus_tag"].values

    def test_kegg_annotations(self, kg_conn):
        """Marker genes should have KEGG annotations."""
        from enrich_utils.extraction import extract_annotations

        df = extract_annotations(MARKER_GENES, "kegg", conn=kg_conn)
        assert len(df) > 0

    def test_no_annotation_gene(self, kg_conn):
        """PMM0030 may have no annotation in some ontologies — should not error."""
        from enrich_utils.extraction import extract_annotations

        df = extract_annotations(["PMM0030"], "cyanorak_role", conn=kg_conn)
        # May be empty — that's fine
        assert isinstance(df, pd.DataFrame)


class TestExtractHierarchy:

    def test_cyanorak_hierarchy(self, kg_conn):
        """CyanoRak should have a 3-level hierarchy."""
        from enrich_utils.extraction import extract_hierarchy

        df = extract_hierarchy("cyanorak_role", max_depth=4, conn=kg_conn)
        assert len(df) > 0
        assert set(df.columns) >= {"child_id", "parent_id", "child_level", "parent_level"}
        # Should have at least 2 levels
        assert df["child_level"].nunique() >= 2

    def test_flat_ontology_returns_empty(self, kg_conn):
        """TIGR role is flat — should return empty hierarchy."""
        from enrich_utils.extraction import extract_hierarchy

        df = extract_hierarchy("tigr_role", conn=kg_conn)
        assert len(df) == 0


class TestEndToEndRollUp:

    def test_cyanorak_rollup_matches_genes_by_ontology(self, kg_conn):
        """Roll-up result for a CyanoRak term should match genes_by_ontology.

        Pick one level-1 term, roll up our annotations, and compare
        gene sets with the MCP tool's answer.
        """
        from enrich_utils.extraction import extract_annotations, extract_hierarchy
        from enrich_utils.hierarchy import roll_up_to_level
        from multiomics_explorer import genes_by_ontology

        # Extract full MED4 annotations and hierarchy
        from multiomics_explorer import genes_by_function
        all_med4 = genes_by_function("*", organism="MED4", limit=None)
        all_locus_tags = [r["locus_tag"] for r in all_med4["results"]]

        annotations = extract_annotations(all_locus_tags, "cyanorak_role", conn=kg_conn)
        hierarchy = extract_hierarchy("cyanorak_role", conn=kg_conn)

        if len(hierarchy) == 0:
            pytest.skip("No hierarchy extracted")

        # Roll up to level 1
        rolled = roll_up_to_level(annotations, hierarchy, target_level=1)

        # Pick the first level-1 term
        if len(rolled) == 0:
            pytest.skip("No genes at level 1")
        test_term = rolled["term_id"].value_counts().index[0]

        our_genes = set(rolled[rolled["term_id"] == test_term]["locus_tag"])

        # Compare with genes_by_ontology
        mcp_result = genes_by_ontology(
            term_ids=[test_term], ontology="cyanorak_role",
            organism="MED4", limit=None, conn=kg_conn,
        )
        mcp_genes = {r["locus_tag"] for r in mcp_result["results"]}

        # Should be equal (or very close — differences = bugs)
        assert our_genes == mcp_genes, (
            f"Mismatch for {test_term}: "
            f"ours={len(our_genes)}, MCP={len(mcp_genes)}, "
            f"only_ours={our_genes - mcp_genes}, "
            f"only_mcp={mcp_genes - our_genes}"
        )
```

- [ ] **Step 3: Commit extraction module and tests**

```bash
git add enrich_utils/extraction.py enrich_utils/tests/test_extraction.py
git commit -m "feat: add extraction.py with KG annotation/hierarchy extraction [Step 0c]"
```

- [ ] **Step 4: Run integration tests interactively (requires Neo4j)**

```bash
uv run pytest enrich_utils/tests/test_extraction.py -v
```

This is the **MCP pause point**. If extraction or roll-up fails:
1. Document what went wrong in `mcp_tool_requirements.md`
2. Decide: fix Cypher, or file MCP enhancement and wait

- [ ] **Step 5: Commit validation results to notebook**

Add notebook entry to `exploration/YYYY-MM-DD-notebook.md`:
- Phase 0c results: which ontologies worked, which failed
- Marker gene annotation traces
- MCP friction points
- Decision: proceed to phase 1 / pause for MCP enhancement

```bash
git add exploration/ mcp_tool_requirements.md
git commit -m "notebook: phase 0c validation results and MCP assessment [Step 0c]"
```

---

### Task 6: Create io.py helper [Spec: enrich_utils modules]

**Files:**
- Create: `enrich_utils/io.py`

- [ ] **Step 1: Implement io.py**

```python
# enrich_utils/io.py
"""Load/save helpers for enrichment analysis data files."""
import pandas as pd
from pathlib import Path


def load_annotations(path: str | Path) -> pd.DataFrame:
    """Load annotations CSV (locus_tag, term_id, term_name)."""
    return pd.read_csv(path)


def load_hierarchy(path: str | Path) -> pd.DataFrame:
    """Load hierarchy CSV (child_id, parent_id, child_level, parent_level)."""
    return pd.read_csv(path)


def load_de(path: str | Path) -> pd.DataFrame:
    """Load DE CSV from v2 analysis."""
    return pd.read_csv(path)


def save_enrichment_results(df: pd.DataFrame, path: str | Path) -> None:
    """Save enrichment results, converting set columns to pipe-delimited strings."""
    out = df.copy()
    for col in out.columns:
        if out[col].apply(lambda x: isinstance(x, set)).any():
            out[col] = out[col].apply(lambda x: "|".join(sorted(x)) if isinstance(x, set) else x)
    out.to_csv(path, index=False)


def load_gene_universes(de_dir: str | Path) -> dict[str, set[str]]:
    """Load gene universes from all de_*.csv files in a directory.

    Returns {filename_stem: set of locus_tags}.
    """
    de_dir = Path(de_dir)
    universes = {}
    for f in sorted(de_dir.glob("de_*.csv")):
        df = pd.read_csv(f)
        universes[f.stem] = set(df["locus_tag"])
    return universes
```

- [ ] **Step 2: Commit**

```bash
git add enrich_utils/io.py
git commit -m "feat: add io.py load/save helpers"
```

---

## Phase 1: Analysis Scripts (scaffolds for interactive analysis)

Phase 1 steps are interactive — the scripts run, then the researcher explores and decides. These tasks create the script scaffolds. The actual do→show→explore→decide cycle happens during execution, not in the plan.

**Prerequisite:** Phase 0 complete and validated. MCP tools sufficient.

### Task 7: Script 01_extract_annotations.py [Spec: Phase 1 Step 1]

**Files:**
- Create: `scripts/01_extract_annotations.py`

- [ ] **Step 1: Write the extraction script**

```python
#!/usr/bin/env python3
"""Step 1: Extract MED4 annotations from KG for all ontologies.

Extracts for the full MED4 genome (not just experiment universes)
to enable genome-level coverage calculations.

Usage:
    uv run scripts/01_extract_annotations.py [--explore]
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

# Add parent dir to path for enrich_utils
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from enrich_utils.extraction import extract_annotations, extract_hierarchy, ONTOLOGY_EDGES
from enrich_utils.io import load_gene_universes

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
LOG_DIR = Path(__file__).resolve().parent.parent / "logs"
V2_DATA = Path(__file__).resolve().parent.parent.parent / "2026-04-08-1038-n_limitation_signature_v2" / "data"

ONTOLOGIES = list(ONTOLOGY_EDGES.keys())
MARKER_GENES = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]


def main():
    parser = argparse.ArgumentParser(description="Extract MED4 ontology annotations")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces")
    args = parser.parse_args()

    log_lines = []

    # 1. Collect gene universes from v2
    print("Loading gene universes from v2...")
    universes = load_gene_universes(V2_DATA)
    universe_df = pd.DataFrame([
        {"experiment": name, "n_genes": len(genes)}
        for name, genes in universes.items()
    ])
    universe_df.to_csv(DATA_DIR / "gene_universes.csv", index=False)
    log_lines.append(f"Gene universes: {len(universes)} experiments")
    for name, genes in universes.items():
        log_lines.append(f"  {name}: {len(genes)} genes")

    # 2. Get full MED4 genome
    print("Extracting full MED4 gene list...")
    from multiomics_explorer import genes_by_function
    all_med4 = genes_by_function("*", organism="MED4", limit=None)
    all_locus_tags = [r["locus_tag"] for r in all_med4["results"]]
    log_lines.append(f"\nFull MED4 genome: {len(all_locus_tags)} genes")

    # 3. Extract annotations per ontology
    for ont in ONTOLOGIES:
        print(f"Extracting {ont}...")
        df = extract_annotations(all_locus_tags, ont)
        out_path = DATA_DIR / f"annotations_{ont}.csv"
        df.to_csv(out_path, index=False)

        n_genes = df["locus_tag"].nunique()
        n_terms = df["term_id"].nunique()
        n_no_terms = len(set(all_locus_tags) - set(df["locus_tag"]))
        log_lines.append(f"\n{ont}: {n_genes} genes, {n_terms} terms, {n_no_terms} unannotated")

        if args.explore:
            for mg in MARKER_GENES:
                mg_terms = df[df["locus_tag"] == mg]
                if len(mg_terms) > 0:
                    terms_str = ", ".join(f"{r['term_id']}:{r['term_name']}" for _, r in mg_terms.iterrows())
                    print(f"  {mg} → {terms_str}")
                else:
                    print(f"  {mg} → (no {ont} annotation)")

    # 4. Extract hierarchies
    for ont in ONTOLOGIES:
        config = ONTOLOGY_EDGES[ont]
        if not config["hierarchy_rels"]:
            continue
        print(f"Extracting hierarchy for {ont}...")
        df = extract_hierarchy(ont, max_depth=4)
        out_path = DATA_DIR / f"hierarchy_{ont}.csv"
        df.to_csv(out_path, index=False)

        n_edges = len(df)
        levels = sorted(set(df["child_level"]) | set(df["parent_level"])) if n_edges > 0 else []
        log_lines.append(f"\n{ont} hierarchy: {n_edges} edges, levels: {levels}")

    # Write log
    LOG_DIR.mkdir(exist_ok=True)
    (LOG_DIR / "01_extract_annotations.log").write_text("\n".join(log_lines))
    print(f"\nDone. Log: logs/01_extract_annotations.log")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add scripts/01_extract_annotations.py
git commit -m "feat: add script 01_extract_annotations.py [Phase 1 Step 1]"
```

---

### Task 8: Script 02_survey_landscape.py [Spec: Phase 1 Step 2]

**Files:**
- Create: `scripts/02_survey_landscape.py`

- [ ] **Step 1: Write the survey script**

```python
#!/usr/bin/env python3
"""Step 2: Characterize annotation landscape across all ontologies.

Reads annotations and hierarchies from step 1, profiles each ontology,
ranks them for enrichment suitability.

Usage:
    uv run scripts/02_survey_landscape.py [--explore]
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from enrich_utils.io import load_annotations, load_hierarchy, load_gene_universes
from enrich_utils.survey import survey_ontology, rank_ontologies
from enrich_utils.extraction import ONTOLOGY_EDGES

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
LOG_DIR = Path(__file__).resolve().parent.parent / "logs"
V2_DATA = Path(__file__).resolve().parent.parent.parent / "2026-04-08-1038-n_limitation_signature_v2" / "data"


def main():
    parser = argparse.ArgumentParser(description="Survey ontology annotation landscape")
    parser.add_argument("--explore", action="store_true")
    args = parser.parse_args()

    log_lines = []

    # Load gene universes
    universes = load_gene_universes(V2_DATA)
    # Use the largest universe (RNA-seq) as the primary assessment
    primary_universe = max(universes.values(), key=len)
    log_lines.append(f"Primary universe: {len(primary_universe)} genes")

    # Profile each ontology
    profiles = {}
    for ont in ONTOLOGY_EDGES:
        ann_path = DATA_DIR / f"annotations_{ont}.csv"
        hier_path = DATA_DIR / f"hierarchy_{ont}.csv"

        if not ann_path.exists():
            log_lines.append(f"\n{ont}: SKIPPED (no annotation file)")
            continue

        annotations = load_annotations(ann_path)
        hierarchy = load_hierarchy(hier_path) if hier_path.exists() else pd.DataFrame()

        profile = survey_ontology(annotations, hierarchy, primary_universe)
        profiles[ont] = profile

        log_lines.append(f"\n{ont}:")
        log_lines.append(f"  Coverage: {profile['coverage']:.3f} ({profile['n_annotated']}/{profile['n_annotated'] + profile['n_unannotated']})")
        for lp in profile["per_level"]:
            log_lines.append(
                f"  Level {lp['level']}: {lp['n_terms_with_genes']} terms, "
                f"genes/term: {lp['min_genes']}-{lp['max_genes']} (median {lp['median_genes']:.0f})"
            )

    # Rank ontologies
    ranking = rank_ontologies(profiles)
    ranking.to_csv(RESULTS_DIR / "ontology_ranking.csv", index=False)

    # Save profiles
    profile_rows = []
    for ont, profile in profiles.items():
        for lp in profile["per_level"]:
            profile_rows.append({"ontology": ont, "coverage": profile["coverage"], **lp})
    pd.DataFrame(profile_rows).to_csv(DATA_DIR / "ontology_profiles.csv", index=False)

    log_lines.append("\n\nRanking:")
    for _, row in ranking.iterrows():
        log_lines.append(
            f"  #{row['rank']} {row['ontology']}: coverage={row['coverage']:.3f}, "
            f"best_level={row['best_level']}, median_genes={row['best_level_median_genes']:.0f}"
        )

    if args.explore:
        print("\n".join(log_lines))

    (LOG_DIR / "02_characterize_landscape.log").write_text("\n".join(log_lines))
    print(f"Done. Ranking: results/ontology_ranking.csv")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add scripts/02_survey_landscape.py
git commit -m "feat: add script 02_survey_landscape.py [Phase 1 Step 2]"
```

---

### Task 9: Scripts 03-05 (pathway definition, enrichment, plotting) [Spec: Phase 1 Steps 3-5]

**Files:**
- Create: `scripts/03_define_pathways.py`
- Create: `scripts/04_run_enrichment.py`
- Create: `scripts/05_plot_results.py`

- [ ] **Step 1: Write 03_define_pathways.py**

```python
#!/usr/bin/env python3
"""Step 3: Build pathway definitions from the selected ontology.

Requires: ontology selected in step 2 (passed as --ontology and --level args).

Usage:
    uv run scripts/03_define_pathways.py --ontology cyanorak_role --level 1 [--min-genes 5] [--explore]
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from enrich_utils.io import load_annotations, load_hierarchy, load_gene_universes
from enrich_utils.survey import build_pathway_definitions, scope_pathways_to_universe

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
LOG_DIR = Path(__file__).resolve().parent.parent / "logs"
V2_DATA = Path(__file__).resolve().parent.parent.parent / "2026-04-08-1038-n_limitation_signature_v2" / "data"

MARKER_GENES = {"PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"}


def main():
    parser = argparse.ArgumentParser(description="Build pathway definitions")
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--level", type=int, required=True)
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--explore", action="store_true")
    args = parser.parse_args()

    log_lines = []

    annotations = load_annotations(DATA_DIR / f"annotations_{args.ontology}.csv")
    hier_path = DATA_DIR / f"hierarchy_{args.ontology}.csv"
    hierarchy = load_hierarchy(hier_path) if hier_path.exists() else pd.DataFrame()

    # Build term_names from hierarchy + annotations
    term_names = dict(zip(annotations["term_id"], annotations["term_name"]))

    pw_defs = build_pathway_definitions(
        annotations, hierarchy, level=args.level,
        min_genes=args.min_genes, term_names=term_names,
    )

    log_lines.append(f"Ontology: {args.ontology}, level: {args.level}, min_genes: {args.min_genes}")
    log_lines.append(f"Pathways defined: {len(pw_defs)}")
    log_lines.append(f"Gene count: min={pw_defs['gene_count'].min()}, "
                     f"median={pw_defs['gene_count'].median():.0f}, "
                     f"max={pw_defs['gene_count'].max()}")

    # Save pathway definitions (convert sets to pipe-delimited for CSV)
    from enrich_utils.io import save_enrichment_results
    save_enrichment_results(pw_defs, DATA_DIR / "pathway_definitions.csv")

    # Scope to each experiment universe
    universes = load_gene_universes(V2_DATA)
    coverage_rows = []
    for exp_name, universe in universes.items():
        scoped = scope_pathways_to_universe(pw_defs, universe)
        for _, row in scoped.iterrows():
            coverage_rows.append({
                "experiment": exp_name,
                "pathway_id": row["pathway_id"],
                "pathway_name": row["pathway_name"],
                "n_in_universe": row["n_in_universe"],
                "n_total": row["gene_count"],
                "coverage": row["coverage"],
            })

    pd.DataFrame(coverage_rows).to_csv(DATA_DIR / "pathway_coverage_per_experiment.csv", index=False)

    if args.explore:
        print("\n".join(log_lines))
        # Trace marker genes
        for _, pw in pw_defs.iterrows():
            markers_in = pw["locus_tags"] & MARKER_GENES
            if markers_in:
                print(f"  {pw['pathway_name']}: markers {markers_in}")

    (LOG_DIR / "03_define_pathways.log").write_text("\n".join(log_lines))
    print(f"Done. {len(pw_defs)} pathways defined.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write 04_run_enrichment.py**

```python
#!/usr/bin/env python3
"""Step 4: Run Fisher's exact enrichment across all experiments.

Usage:
    uv run scripts/04_run_enrichment.py [--explore]
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from enrich_utils.io import load_de, save_enrichment_results
from enrich_utils.enrichment import run_enrichment_all_timepoints

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
LOG_DIR = Path(__file__).resolve().parent.parent / "logs"
V2_DATA = Path(__file__).resolve().parent.parent.parent / "2026-04-08-1038-n_limitation_signature_v2" / "data"

MARKER_PATHWAYS = ["nitrogen", "photosynth", "transport", "ribosom"]


def main():
    parser = argparse.ArgumentParser(description="Run pathway enrichment")
    parser.add_argument("--explore", action="store_true")
    args = parser.parse_args()

    log_lines = []

    # Load pathway definitions (convert pipe-delimited back to sets)
    pw_defs = pd.read_csv(DATA_DIR / "pathway_definitions.csv")
    pw_defs["locus_tags"] = pw_defs["locus_tags"].apply(lambda x: set(x.split("|")) if pd.notna(x) else set())

    # Run enrichment per experiment
    all_results = []
    for de_file in sorted(V2_DATA.glob("de_*.csv")):
        de_df = load_de(de_file)
        table_scope = de_df["table_scope"].iloc[0] if "table_scope" in de_df.columns else "all_detected_genes"
        exp_name = de_file.stem

        result = run_enrichment_all_timepoints(de_df, pw_defs, table_scope)
        result["experiment"] = exp_name
        all_results.append(result)

        n_sig = (result["padj"] < 0.05).sum() if "padj" in result.columns else 0
        log_lines.append(f"{exp_name}: {len(result)} tests, {n_sig} significant (padj<0.05)")

    combined = pd.concat(all_results, ignore_index=True)
    save_enrichment_results(combined, RESULTS_DIR / "enrichment_all.csv")

    significant = combined[combined["padj"] < 0.05] if "padj" in combined.columns else combined.head(0)
    save_enrichment_results(significant, RESULTS_DIR / "enrichment_significant.csv")

    log_lines.append(f"\nTotal: {len(combined)} tests, {len(significant)} significant")

    if args.explore:
        print("\n".join(log_lines))
        # Trace marker pathways
        for marker in MARKER_PATHWAYS:
            hits = significant[significant["pathway_name"].str.lower().str.contains(marker, na=False)]
            if len(hits) > 0:
                print(f"\n  '{marker}' enriched in:")
                for _, row in hits.iterrows():
                    print(f"    {row['experiment']} {row.get('timepoint', '')} {row['direction']}: "
                          f"padj={row['padj']:.4f}, OR={row['odds_ratio']:.2f}")

    (LOG_DIR / "04_run_enrichment.log").write_text("\n".join(log_lines))
    print(f"Done. Results: results/enrichment_all.csv ({len(combined)} rows)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Write 05_plot_results.py (stub — details depend on enrichment results)**

```python
#!/usr/bin/env python3
"""Step 5: Visualize enrichment results.

Heatmaps and comparison figures.

Usage:
    uv run scripts/05_plot_results.py
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


def plot_heatmap(df, direction, out_path):
    """Plot enrichment heatmap: pathways × experiments, colored by -log10(padj)."""
    subset = df[df["direction"] == direction].copy()
    if len(subset) == 0:
        print(f"No {direction} enrichment results to plot.")
        return

    subset["neg_log10_padj"] = -np.log10(subset["padj"].clip(lower=1e-10))
    subset["label"] = subset["experiment"] + " " + subset["timepoint"].fillna("").astype(str)

    pivot = subset.pivot_table(
        index="pathway_name", columns="label",
        values="neg_log10_padj", fill_value=0,
    )

    fig, ax = plt.subplots(figsize=(max(12, len(pivot.columns) * 0.8), max(6, len(pivot) * 0.4)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title(f"Pathway enrichment ({direction}): -log10(padj)")
    plt.colorbar(im, ax=ax, label="-log10(padj)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


def main():
    enrichment = pd.read_csv(RESULTS_DIR / "enrichment_all.csv")

    plot_heatmap(enrichment, "up", RESULTS_DIR / "heatmap_enrichment_up.png")
    plot_heatmap(enrichment, "down", RESULTS_DIR / "heatmap_enrichment_down.png")

    print("Done.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Commit all scripts**

```bash
git add scripts/03_define_pathways.py scripts/04_run_enrichment.py scripts/05_plot_results.py
git commit -m "feat: add scripts 03-05 for pathway definition, enrichment, and plotting [Phase 1 Steps 3-5]"
```

---

## Plan Self-Review

**Spec coverage check:**
- [x] Phase 0a (enrichment toy tests): Task 2
- [x] Phase 0b (hierarchy roll-up, survey, pathway defs): Tasks 3-4
- [x] Phase 0c (KG extraction, marker validation, MCP pause point): Task 5
- [x] Phase 1 Step 1 (extract annotations): Task 7
- [x] Phase 1 Step 2 (survey landscape): Task 8
- [x] Phase 1 Steps 3-5 (pathways, enrichment, plots): Task 9
- [x] Phase 1 Step 6 (interpret/document): Not a code task — done during interactive analysis
- [x] enrich_utils modules: enrichment.py (Task 2), hierarchy.py (Task 3), survey.py (Task 4), extraction.py (Task 5), io.py (Task 6)
- [x] Directory scaffold: Task 1
- [x] Marker genes traced: in test_extraction.py and script `--explore` flags
- [x] test_type mapping: in test_enrichment.py Case 9
- [x] DAG convergence, multi-annotation, non-leaf annotation: in test_survey.py TestRollUp
- [x] Notebook-commit gates: noted in Task 5 Step 5 and task descriptions
- [x] MCP pause point: Task 5 Step 4

**Placeholder scan:** No TBDs, no "implement later", no "similar to Task N." All code blocks are complete.

**Type consistency:** `run_enrichment` signature consistent between test and implementation. `roll_up_to_level` parameters match. `survey_ontology` return shape matches test expectations. `pathway_defs` DataFrame schema consistent across survey.py, enrichment.py, and scripts.
