"""Tests for hierarchy.py (roll-up) and survey.py (profiling + pathway defs).

Shared synthetic DAG:

    Level 0 (roots): R1, R2
    Level 1 (mid):   M1->R1, M2->R1, M3->R2, M4->R2
    Level 2 (leaves): L1->M1, L2->M1, L3->M2, L4->M2, L5->M3,
                       L6->M3, L6->M4 (DAG: L6 has two parents)

Gene annotations:
    G01 -> L1            (leaf)
    G02 -> L1, L3        (two leaves, different branches under R1)
    G03 -> M2            (mid-level directly)
    G04 -> L5            (leaf under R2)
    G05 -> L6            (leaf with two parents: M3 and M4, both under R2)
    G06 -> R1            (root directly)
    G07 -> (no annotation)
    G08 -> L1, M1        (both leaf and its parent)
"""

import pandas as pd
import pytest

from enrich_utils.hierarchy import roll_up_to_level
from enrich_utils.survey import (
    build_pathway_definitions,
    rank_ontologies,
    scope_pathways_to_universe,
    survey_ontology,
)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

GENE_UNIVERSE = {"G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08"}


@pytest.fixture
def hierarchy_df():
    """Build the synthetic DAG hierarchy as a DataFrame."""
    rows = [
        # Level 1 -> Level 0
        {"child_id": "M1", "parent_id": "R1", "child_level": 1, "parent_level": 0},
        {"child_id": "M2", "parent_id": "R1", "child_level": 1, "parent_level": 0},
        {"child_id": "M3", "parent_id": "R2", "child_level": 1, "parent_level": 0},
        {"child_id": "M4", "parent_id": "R2", "child_level": 1, "parent_level": 0},
        # Level 2 -> Level 1
        {"child_id": "L1", "parent_id": "M1", "child_level": 2, "parent_level": 1},
        {"child_id": "L2", "parent_id": "M1", "child_level": 2, "parent_level": 1},
        {"child_id": "L3", "parent_id": "M2", "child_level": 2, "parent_level": 1},
        {"child_id": "L4", "parent_id": "M2", "child_level": 2, "parent_level": 1},
        {"child_id": "L5", "parent_id": "M3", "child_level": 2, "parent_level": 1},
        {"child_id": "L6", "parent_id": "M3", "child_level": 2, "parent_level": 1},
        {"child_id": "L6", "parent_id": "M4", "child_level": 2, "parent_level": 1},  # DAG edge
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def annotations_df():
    """Gene annotations to leaf/mid/root terms."""
    rows = [
        {"locus_tag": "G01", "term_id": "L1"},
        {"locus_tag": "G02", "term_id": "L1"},
        {"locus_tag": "G02", "term_id": "L3"},
        {"locus_tag": "G03", "term_id": "M2"},
        {"locus_tag": "G04", "term_id": "L5"},
        {"locus_tag": "G05", "term_id": "L6"},
        {"locus_tag": "G06", "term_id": "R1"},
        # G07: no annotation
        {"locus_tag": "G08", "term_id": "L1"},
        {"locus_tag": "G08", "term_id": "M1"},
    ]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Roll-up tests (Task 3)
# ---------------------------------------------------------------------------

class TestRollUpToLevel:
    """7 edge-case tests for roll_up_to_level."""

    def test_leaf_propagates_up(self, annotations_df, hierarchy_df):
        """Case 1: G01 -> L1. At level 1, G01 should appear under M1.
        At level 0, G01 should appear under R1."""
        result_l1 = roll_up_to_level(annotations_df, hierarchy_df, target_level=1)
        genes_m1 = set(result_l1[result_l1["term_id"] == "M1"]["locus_tag"])
        assert "G01" in genes_m1, "G01 (annotated to L1) must propagate to M1 at level 1"

        result_l0 = roll_up_to_level(annotations_df, hierarchy_df, target_level=0)
        genes_r1 = set(result_l0[result_l0["term_id"] == "R1"]["locus_tag"])
        assert "G01" in genes_r1, "G01 must propagate to R1 at level 0"

    def test_mid_level_does_not_propagate_down(self, annotations_df, hierarchy_df):
        """Case 2: G03 -> M2. At level 2, G03 must NOT appear in L3 or L4 (no downward
        propagation). At level 1, G03 appears in M2. At level 0, G03 appears in R1."""
        result_l2 = roll_up_to_level(annotations_df, hierarchy_df, target_level=2)
        for leaf in ("L3", "L4"):
            genes = set(result_l2[result_l2["term_id"] == leaf]["locus_tag"])
            assert "G03" not in genes, f"G03 should NOT appear in {leaf} (no downward propagation)"

        result_l1 = roll_up_to_level(annotations_df, hierarchy_df, target_level=1)
        genes_m2 = set(result_l1[result_l1["term_id"] == "M2"]["locus_tag"])
        assert "G03" in genes_m2, "G03 (annotated to M2) must appear at M2 level 1"

        result_l0 = roll_up_to_level(annotations_df, hierarchy_df, target_level=0)
        genes_r1 = set(result_l0[result_l0["term_id"] == "R1"]["locus_tag"])
        assert "G03" in genes_r1, "G03 must appear at R1 level 0"

    def test_two_branch_gene(self, annotations_df, hierarchy_df):
        """Case 3: G02 -> L1 (under M1) and L3 (under M2), both under R1.
        At level 1, G02 must appear in both M1 and M2."""
        result_l1 = roll_up_to_level(annotations_df, hierarchy_df, target_level=1)
        genes_m1 = set(result_l1[result_l1["term_id"] == "M1"]["locus_tag"])
        genes_m2 = set(result_l1[result_l1["term_id"] == "M2"]["locus_tag"])
        assert "G02" in genes_m1, "G02 must appear in M1 at level 1"
        assert "G02" in genes_m2, "G02 must appear in M2 at level 1"

    def test_leaf_and_parent_counted_once(self, annotations_df, hierarchy_df):
        """Case 4: G08 -> L1, M1. L1's parent is M1. At level 1, G08 must appear
        exactly once under M1."""
        result_l1 = roll_up_to_level(annotations_df, hierarchy_df, target_level=1)
        m1_rows = result_l1[(result_l1["term_id"] == "M1") & (result_l1["locus_tag"] == "G08")]
        assert len(m1_rows) == 1, f"G08 should appear exactly once in M1; got {len(m1_rows)}"

    def test_dag_convergence_counted_once(self, annotations_df, hierarchy_df):
        """Case 5: G05 -> L6. L6 has two parents: M3 and M4, both under R2.
        At level 1, G05 must appear in both M3 and M4 (one entry each).
        At level 0, G05 must appear exactly once under R2."""
        result_l1 = roll_up_to_level(annotations_df, hierarchy_df, target_level=1)
        genes_m3 = set(result_l1[result_l1["term_id"] == "M3"]["locus_tag"])
        genes_m4 = set(result_l1[result_l1["term_id"] == "M4"]["locus_tag"])
        assert "G05" in genes_m3, "G05 must appear in M3 at level 1 (first parent of L6)"
        assert "G05" in genes_m4, "G05 must appear in M4 at level 1 (second parent of L6)"

        result_l0 = roll_up_to_level(annotations_df, hierarchy_df, target_level=0)
        r2_rows = result_l0[(result_l0["term_id"] == "R2") & (result_l0["locus_tag"] == "G05")]
        assert len(r2_rows) == 1, f"G05 should appear exactly once in R2; got {len(r2_rows)}"

    def test_root_annotation_present_at_level0_absent_below(self, annotations_df, hierarchy_df):
        """Case 6: G06 -> R1. At level 0, G06 appears under R1.
        At level 1 and level 2, G06 must not appear (no downward propagation)."""
        result_l0 = roll_up_to_level(annotations_df, hierarchy_df, target_level=0)
        genes_r1 = set(result_l0[result_l0["term_id"] == "R1"]["locus_tag"])
        assert "G06" in genes_r1, "G06 (annotated to R1) must appear at level 0"

        result_l1 = roll_up_to_level(annotations_df, hierarchy_df, target_level=1)
        all_l1_genes = set(result_l1["locus_tag"])
        assert "G06" not in all_l1_genes, "G06 must not appear at level 1 (no downward propagation)"

        result_l2 = roll_up_to_level(annotations_df, hierarchy_df, target_level=2)
        all_l2_genes = set(result_l2["locus_tag"])
        assert "G06" not in all_l2_genes, "G06 must not appear at level 2"

    def test_unannotated_gene_absent(self, annotations_df, hierarchy_df):
        """Case 7: G07 has no annotation. Must be absent from all roll-up levels."""
        for level in (0, 1, 2):
            result = roll_up_to_level(annotations_df, hierarchy_df, target_level=level)
            assert "G07" not in set(result["locus_tag"]), \
                f"Unannotated G07 must not appear at level {level}"


# ---------------------------------------------------------------------------
# Survey tests (Task 4)
# ---------------------------------------------------------------------------

class TestSurveyOntology:
    def test_coverage(self, annotations_df, hierarchy_df):
        """7/8 genes annotated -> coverage = 0.875."""
        profile = survey_ontology(annotations_df, hierarchy_df, GENE_UNIVERSE)
        assert abs(profile["coverage"] - 7 / 8) < 1e-9
        assert profile["n_annotated"] == 7
        assert profile["n_unannotated"] == 1

    def test_level0_term_count(self, annotations_df, hierarchy_df):
        """At level 0, exactly 2 terms have genes (R1 and R2)."""
        profile = survey_ontology(annotations_df, hierarchy_df, GENE_UNIVERSE)
        l0 = next(p for p in profile["per_level"] if p["level"] == 0)
        assert l0["n_terms_with_genes"] == 2

    def test_level0_term_sizes(self, annotations_df, hierarchy_df):
        """R1 has genes {G01,G02,G03,G06,G08} = 5; R2 has {G04,G05} = 2."""
        profile = survey_ontology(annotations_df, hierarchy_df, GENE_UNIVERSE)
        l0 = next(p for p in profile["per_level"] if p["level"] == 0)
        assert l0["min_genes"] == 2
        assert l0["max_genes"] == 5


class TestSurveyGeneCoverage:
    def test_gene_coverage_at_level(self, annotations_df, hierarchy_df):
        """gene_coverage per level: level 0 should have all 7 annotated genes."""
        profile = survey_ontology(annotations_df, hierarchy_df, GENE_UNIVERSE)
        l0 = next(p for p in profile["per_level"] if p["level"] == 0)
        assert l0["n_genes_at_level"] == 7  # all annotated genes roll up to roots
        assert abs(l0["gene_coverage"] - 1.0) < 1e-9

    def test_gene_coverage_decreases_at_deeper_levels(self, annotations_df, hierarchy_df):
        """Deeper levels may cover fewer genes (e.g., G06 at root only)."""
        profile = survey_ontology(annotations_df, hierarchy_df, GENE_UNIVERSE)
        l0 = next(p for p in profile["per_level"] if p["level"] == 0)
        l1 = next(p for p in profile["per_level"] if p["level"] == 1)
        # G06 is at root only — absent from level 1
        assert l1["n_genes_at_level"] < l0["n_genes_at_level"]
        assert l1["gene_coverage"] < l0["gene_coverage"]


class TestRankOntologies:
    def test_ranking_returns_dataframe(self):
        """rank_ontologies takes a dict of profiles and returns a DataFrame."""
        profiles = {
            "GO": {"coverage": 0.9, "per_level": [{"level": 1, "n_genes_at_level": 900,
                                                     "gene_coverage": 0.9,
                                                     "n_terms_with_genes": 10,
                                                     "min_genes": 1, "q1_genes": 5,
                                                     "median_genes": 20, "q3_genes": 50,
                                                     "max_genes": 200}]},
            "KEGG": {"coverage": 0.5, "per_level": [{"level": 1, "n_genes_at_level": 450,
                                                       "gene_coverage": 0.9,
                                                       "n_terms_with_genes": 5,
                                                       "min_genes": 3, "q1_genes": 8,
                                                       "median_genes": 15, "q3_genes": 30,
                                                       "max_genes": 100}]},
        }
        result = rank_ontologies(profiles)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert "ontology" in result.columns

    def test_higher_coverage_ranks_first(self):
        """Ontology with higher coverage should rank first when gene_coverage is equal."""
        profiles = {
            "HighCov": {"coverage": 0.9, "per_level": [{"level": 1, "n_genes_at_level": 800,
                                                          "gene_coverage": 0.9,
                                                          "n_terms_with_genes": 5,
                                                          "min_genes": 2, "q1_genes": 5,
                                                          "median_genes": 15, "q3_genes": 30,
                                                          "max_genes": 100}]},
            "LowCov": {"coverage": 0.3, "per_level": [{"level": 1, "n_genes_at_level": 250,
                                                         "gene_coverage": 0.9,
                                                         "n_terms_with_genes": 5,
                                                         "min_genes": 2, "q1_genes": 5,
                                                         "median_genes": 15, "q3_genes": 30,
                                                         "max_genes": 100}]},
        }
        result = rank_ontologies(profiles)
        assert result.iloc[0]["ontology"] == "HighCov"

    def test_low_gene_coverage_level_rejected(self):
        """A level with good median but low gene coverage should not be selected."""
        profiles = {
            "DeepOnly": {"coverage": 0.8, "per_level": [
                {"level": 0, "n_genes_at_level": 800, "gene_coverage": 1.0,
                 "n_terms_with_genes": 3, "min_genes": 200, "q1_genes": 250,
                 "median_genes": 300, "q3_genes": 350, "max_genes": 400},  # too broad
                {"level": 2, "n_genes_at_level": 160, "gene_coverage": 0.2,
                 "n_terms_with_genes": 20, "min_genes": 5, "q1_genes": 8,
                 "median_genes": 12, "q3_genes": 20, "max_genes": 40},  # good median, bad coverage
            ]},
        }
        result = rank_ontologies(profiles)
        # Level 2 has gene_coverage 0.2 (< 0.5 threshold) — should not qualify
        assert result.iloc[0]["best_level"] is None or pd.isna(result.iloc[0]["best_level"])


class TestBuildPathwayDefinitions:
    def test_level1_m1_contains_correct_genes(self, annotations_df, hierarchy_df):
        """At level 1 with min_genes=1, M1 must contain G01, G02, G08."""
        defs = build_pathway_definitions(annotations_df, hierarchy_df, level=1, min_genes=1)
        m1_row = defs[defs["pathway_id"] == "M1"]
        assert len(m1_row) == 1
        m1_genes = m1_row.iloc[0]["locus_tags"]
        assert {"G01", "G02", "G08"}.issubset(m1_genes), \
            f"M1 should contain G01, G02, G08; got {m1_genes}"

    def test_min_genes_filter(self, annotations_df, hierarchy_df):
        """With min_genes=3, only pathways with >=3 genes survive."""
        defs = build_pathway_definitions(annotations_df, hierarchy_df, level=1, min_genes=3)
        assert (defs["gene_count"] >= 3).all()
        # M2 has G02, G03 = 2 genes; should be filtered out
        assert "M2" not in set(defs["pathway_id"]) or \
               defs[defs["pathway_id"] == "M2"].iloc[0]["gene_count"] >= 3

    def test_returns_required_columns(self, annotations_df, hierarchy_df):
        """build_pathway_definitions must return pathway_id, pathway_name, locus_tags, gene_count."""
        defs = build_pathway_definitions(annotations_df, hierarchy_df, level=1, min_genes=1)
        for col in ("pathway_id", "pathway_name", "locus_tags", "gene_count"):
            assert col in defs.columns, f"Missing column: {col}"

    def test_term_names_used_when_provided(self, annotations_df, hierarchy_df):
        """When term_names dict is provided, pathway_name should come from it."""
        term_names = {"M1": "Metabolism group 1", "M2": "Metabolism group 2"}
        defs = build_pathway_definitions(
            annotations_df, hierarchy_df, level=1, min_genes=1, term_names=term_names
        )
        m1_row = defs[defs["pathway_id"] == "M1"]
        if len(m1_row) > 0:
            assert m1_row.iloc[0]["pathway_name"] == "Metabolism group 1"

    def test_gene_count_matches_locus_tags_length(self, annotations_df, hierarchy_df):
        """gene_count must equal len(locus_tags) for every row."""
        defs = build_pathway_definitions(annotations_df, hierarchy_df, level=1, min_genes=1)
        for _, row in defs.iterrows():
            assert row["gene_count"] == len(row["locus_tags"]), \
                f"gene_count mismatch for {row['pathway_id']}"


class TestScopePathwaysToUniverse:
    def test_r1_with_small_universe(self, annotations_df, hierarchy_df):
        """R1 contains {G01,G02,G03,G06,G08} = 5 genes.
        Universe {G01,G02,G03} -> n_in_universe=3, coverage=3/5."""
        defs_l0 = build_pathway_definitions(annotations_df, hierarchy_df, level=0, min_genes=1)
        small_universe = {"G01", "G02", "G03"}
        scoped = scope_pathways_to_universe(defs_l0, small_universe)

        r1_row = scoped[scoped["pathway_id"] == "R1"]
        assert len(r1_row) == 1
        row = r1_row.iloc[0]
        assert row["n_in_universe"] == 3
        assert abs(row["coverage"] - 3 / 5) < 1e-9
        assert row["locus_tags_in_universe"] == {"G01", "G02", "G03"}

    def test_scope_adds_required_columns(self, annotations_df, hierarchy_df):
        """scope_pathways_to_universe must add n_in_universe, coverage, locus_tags_in_universe."""
        defs = build_pathway_definitions(annotations_df, hierarchy_df, level=0, min_genes=1)
        scoped = scope_pathways_to_universe(defs, GENE_UNIVERSE)
        for col in ("n_in_universe", "coverage", "locus_tags_in_universe"):
            assert col in scoped.columns, f"Missing column: {col}"

    def test_scope_does_not_exceed_universe(self, annotations_df, hierarchy_df):
        """locus_tags_in_universe must be a subset of gene_universe."""
        defs = build_pathway_definitions(annotations_df, hierarchy_df, level=1, min_genes=1)
        small_universe = {"G01", "G02"}
        scoped = scope_pathways_to_universe(defs, small_universe)
        for _, row in scoped.iterrows():
            assert row["locus_tags_in_universe"].issubset(small_universe)
