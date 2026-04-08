"""Toy data tests for signature construction.

Hand-calculated expected values for every operation.
Seed for a real test suite during productization.
"""

import pandas as pd
import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from sig_utils.signature import summarize_per_gene, intersect_references


# === Toy DE data ===

def make_toy_study_a():
    """Study A: 5 genes, 2 timepoints.

    G1: up at both tp1 (rank_up=2) and tp2 (rank_up=1) → up, best_rank=1
    G2: down at tp1 (rank_down=3), not sig at tp2 → down, best_rank=3
    G3: up at tp1 (rank_up=4), not sig at tp2 → up, best_rank=4
    G4: up at tp1 (rank_up=3), down at tp2 (rank_down=2) → tie, break by rank: down wins (2<3)
    G5: not significant at either → excluded
    """
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "geneA", "timepoint": "tp1",
         "log2fc": 2.0, "expression_status": "significant_up",
         "rank": 3, "rank_up": 2, "rank_down": pd.NA},
        {"locus_tag": "G1", "gene_name": "geneA", "timepoint": "tp2",
         "log2fc": 3.0, "expression_status": "significant_up",
         "rank": 1, "rank_up": 1, "rank_down": pd.NA},
        {"locus_tag": "G2", "gene_name": "geneB", "timepoint": "tp1",
         "log2fc": -1.5, "expression_status": "significant_down",
         "rank": 5, "rank_up": pd.NA, "rank_down": 3},
        {"locus_tag": "G2", "gene_name": "geneB", "timepoint": "tp2",
         "log2fc": -0.3, "expression_status": "not_significant",
         "rank": 10, "rank_up": pd.NA, "rank_down": pd.NA},
        {"locus_tag": "G3", "gene_name": "geneC", "timepoint": "tp1",
         "log2fc": 1.0, "expression_status": "significant_up",
         "rank": 7, "rank_up": 4, "rank_down": pd.NA},
        {"locus_tag": "G3", "gene_name": "geneC", "timepoint": "tp2",
         "log2fc": 0.1, "expression_status": "not_significant",
         "rank": 20, "rank_up": pd.NA, "rank_down": pd.NA},
        {"locus_tag": "G4", "gene_name": "geneD", "timepoint": "tp1",
         "log2fc": 1.2, "expression_status": "significant_up",
         "rank": 6, "rank_up": 3, "rank_down": pd.NA},
        {"locus_tag": "G4", "gene_name": "geneD", "timepoint": "tp2",
         "log2fc": -1.8, "expression_status": "significant_down",
         "rank": 4, "rank_up": pd.NA, "rank_down": 2},
        {"locus_tag": "G5", "gene_name": "geneE", "timepoint": "tp1",
         "log2fc": 0.05, "expression_status": "not_significant",
         "rank": 50, "rank_up": pd.NA, "rank_down": pd.NA},
        {"locus_tag": "G5", "gene_name": "geneE", "timepoint": "tp2",
         "log2fc": -0.1, "expression_status": "not_significant",
         "rank": 45, "rank_up": pd.NA, "rank_down": pd.NA},
    ])


def make_toy_study_b():
    """Study B: 4 genes, 1 timepoint. G3 absent from dataset entirely.

    G1: up, rank_up=1 → concordant with study A
    G2: down, rank_down=2 → concordant with study A
    G4: up, rank_up=2 → may be discordant depending on A's tie-break
    G5: not significant
    G6: down, rank_down=1 → only in study B
    """
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "geneA", "timepoint": "tp1",
         "log2fc": 4.0, "expression_status": "significant_up",
         "rank": 1, "rank_up": 1, "rank_down": pd.NA},
        {"locus_tag": "G2", "gene_name": "geneB", "timepoint": "tp1",
         "log2fc": -2.0, "expression_status": "significant_down",
         "rank": 3, "rank_up": pd.NA, "rank_down": 2},
        {"locus_tag": "G4", "gene_name": "geneD", "timepoint": "tp1",
         "log2fc": 1.5, "expression_status": "significant_up",
         "rank": 4, "rank_up": 2, "rank_down": pd.NA},
        {"locus_tag": "G5", "gene_name": "geneE", "timepoint": "tp1",
         "log2fc": 0.1, "expression_status": "not_significant",
         "rank": 20, "rank_up": pd.NA, "rank_down": pd.NA},
        {"locus_tag": "G6", "gene_name": "geneF", "timepoint": "tp1",
         "log2fc": -3.0, "expression_status": "significant_down",
         "rank": 2, "rank_up": pd.NA, "rank_down": 1},
    ])


class TestSummarizePerGene:
    def test_clear_majority_direction(self):
        """G1 is up at both timepoints → direction = up."""
        result = summarize_per_gene(make_toy_study_a())
        g1 = result[result["locus_tag"] == "G1"].iloc[0]
        assert g1["direction"] == "up"

    def test_best_directional_rank(self):
        """G1 has rank_up=2 at tp1, rank_up=1 at tp2 → best_dir_rank = 1."""
        result = summarize_per_gene(make_toy_study_a())
        g1 = result[result["locus_tag"] == "G1"].iloc[0]
        assert g1["best_dir_rank"] == 1

    def test_single_timepoint_significant(self):
        """G2 is significant_down at tp1 only → direction = down."""
        result = summarize_per_gene(make_toy_study_a())
        g2 = result[result["locus_tag"] == "G2"].iloc[0]
        assert g2["direction"] == "down"

    def test_not_significant_excluded(self):
        """G5 is not significant at any timepoint → excluded from summary."""
        result = summarize_per_gene(make_toy_study_a())
        assert "G5" not in result["locus_tag"].values

    def test_gene_count(self):
        """4 genes significant in at least one timepoint (G1, G2, G3, G4)."""
        result = summarize_per_gene(make_toy_study_a())
        assert len(result) == 4

    def test_tie_broken_by_best_rank(self):
        """G4: up at tp1 (rank_up=3), down at tp2 (rank_down=2).
        Tie 1:1 → break by best rank: down=2 < up=3 → direction = down.
        """
        result = summarize_per_gene(make_toy_study_a())
        g4 = result[result["locus_tag"] == "G4"].iloc[0]
        assert g4["direction"] == "down"
        assert g4["best_dir_rank"] == 2

    def test_peak_timepoint(self):
        """G1 best rank at tp2 → peak_timepoint = tp2."""
        result = summarize_per_gene(make_toy_study_a())
        g1 = result[result["locus_tag"] == "G1"].iloc[0]
        assert g1["peak_timepoint"] == "tp2"


class TestIntersectReferences:
    def test_concordant_core(self):
        """G1 (up in both) and G2 (down in both) → core."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        core, _ = intersect_references(a_summary, b_summary)
        core_tags = set(core["locus_tag"])
        assert "G1" in core_tags
        assert "G2" in core_tags

    def test_core_direction(self):
        """Core genes have correct direction."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        core, _ = intersect_references(a_summary, b_summary)
        g1 = core[core["locus_tag"] == "G1"].iloc[0]
        assert g1["direction"] == "up"
        g2 = core[core["locus_tag"] == "G2"].iloc[0]
        assert g2["direction"] == "down"

    def test_cross_study_best_rank(self):
        """G1: study_a rank=1, study_b rank=1 → cross_study = 1.
        G2: study_a rank=3, study_b rank=2 → cross_study = 2.
        """
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        core, _ = intersect_references(a_summary, b_summary)
        g1 = core[core["locus_tag"] == "G1"].iloc[0]
        assert g1["cross_study_best_dir_rank"] == 1
        g2 = core[core["locus_tag"] == "G2"].iloc[0]
        assert g2["cross_study_best_dir_rank"] == 2

    def test_discordant_excluded(self):
        """G4: down in A (tie-break), up in B → discordant → not in core or extended."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        core, extended = intersect_references(a_summary, b_summary)
        all_tags = set(core["locus_tag"]) | set(extended["locus_tag"])
        assert "G4" not in all_tags

    def test_b_only_extended(self):
        """G6 is only in study B → extended, study_b_only."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        _, extended = intersect_references(a_summary, b_summary)
        g6 = extended[extended["locus_tag"] == "G6"]
        assert len(g6) == 1
        assert "study_b_only" in g6.iloc[0]["signature_type"]

    def test_a_only_b_absent(self):
        """G3 is in study A but absent from study B dataset → a_only_b_absent."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        study_b_all_tags = {"G1", "G2", "G4", "G5", "G6"}  # G3 absent
        _, extended = intersect_references(
            a_summary, b_summary,
            study_b_all_locus_tags=study_b_all_tags,
        )
        g3 = extended[extended["locus_tag"] == "G3"]
        assert len(g3) == 1
        assert "absent" in g3.iloc[0]["signature_type"]

    def test_a_only_b_ns(self):
        """If G3 were in study B's dataset (but not significant), classify as ns."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        study_b_all_tags = {"G1", "G2", "G3", "G4", "G5", "G6"}  # G3 present
        _, extended = intersect_references(
            a_summary, b_summary,
            study_b_all_locus_tags=study_b_all_tags,
        )
        g3 = extended[extended["locus_tag"] == "G3"]
        assert len(g3) == 1
        assert "_ns" in g3.iloc[0]["signature_type"]
