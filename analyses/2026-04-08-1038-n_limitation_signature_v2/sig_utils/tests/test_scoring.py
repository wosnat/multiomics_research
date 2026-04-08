"""Toy data tests for scoring and permutation.

Hand-calculated expected values for rank score and permutation test.
"""

import pandas as pd
import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from sig_utils.scoring import apply_signature, rank_score, permutation_test, score_with_significance


def make_toy_signature():
    """4-gene signature: 2 up, 2 down."""
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "up1", "direction": "up"},
        {"locus_tag": "G2", "gene_name": "down1", "direction": "down"},
        {"locus_tag": "G3", "gene_name": "up2", "direction": "up"},
        {"locus_tag": "G4", "gene_name": "down2", "direction": "down"},
    ])


def _bg_genes(n=96):
    """Background genes to pad a DE dataset."""
    return [
        {"locus_tag": f"BG{i}", "gene_name": f"bg{i}", "log2fc": 0.1,
         "expression_status": "not_significant", "rank": 20 + i,
         "rank_up": pd.NA, "rank_down": pd.NA, "timepoint": "tp1"}
        for i in range(n)
    ]


def make_toy_de_all_concordant():
    """All 4 signature genes significant in expected direction.

    G1: up, rank_up=1. G2: down, rank_down=1.
    G3: up, rank_up=2. G4: down, rank_down=3.
    Total genes = 100.
    """
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "up1", "log2fc": 3.0,
         "expression_status": "significant_up", "rank": 2,
         "rank_up": 1, "rank_down": pd.NA, "timepoint": "tp1"},
        {"locus_tag": "G2", "gene_name": "down1", "log2fc": -2.5,
         "expression_status": "significant_down", "rank": 3,
         "rank_up": pd.NA, "rank_down": 1, "timepoint": "tp1"},
        {"locus_tag": "G3", "gene_name": "up2", "log2fc": 1.5,
         "expression_status": "significant_up", "rank": 5,
         "rank_up": 2, "rank_down": pd.NA, "timepoint": "tp1"},
        {"locus_tag": "G4", "gene_name": "down2", "log2fc": -1.0,
         "expression_status": "significant_down", "rank": 8,
         "rank_up": pd.NA, "rank_down": 3, "timepoint": "tp1"},
        *_bg_genes(96),
    ])


def make_toy_de_mixed():
    """G1 concordant, G2 concordant, G3 reversed, G4 not significant.

    G1: sig_up rank_up=2. G2: sig_down rank_down=3.
    G3: sig_down rank_down=5 (REVERSED — signature says up).
    G4: not_significant.
    Total genes = 100.
    """
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "up1", "log2fc": 2.0,
         "expression_status": "significant_up", "rank": 3,
         "rank_up": 2, "rank_down": pd.NA, "timepoint": "tp1"},
        {"locus_tag": "G2", "gene_name": "down1", "log2fc": -1.5,
         "expression_status": "significant_down", "rank": 5,
         "rank_up": pd.NA, "rank_down": 3, "timepoint": "tp1"},
        {"locus_tag": "G3", "gene_name": "up2", "log2fc": -1.0,
         "expression_status": "significant_down", "rank": 8,
         "rank_up": pd.NA, "rank_down": 5, "timepoint": "tp1"},
        {"locus_tag": "G4", "gene_name": "down2", "log2fc": 0.2,
         "expression_status": "not_significant", "rank": 50,
         "rank_up": pd.NA, "rank_down": pd.NA, "timepoint": "tp1"},
        *_bg_genes(96),
    ])


def make_toy_de_gene_absent():
    """G1 and G3 present, G2 absent from dataset entirely, G4 not sig.

    Total genes = 99 (G2 missing).
    """
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "up1", "log2fc": 2.0,
         "expression_status": "significant_up", "rank": 3,
         "rank_up": 2, "rank_down": pd.NA, "timepoint": "tp1"},
        {"locus_tag": "G3", "gene_name": "up2", "log2fc": 1.0,
         "expression_status": "significant_up", "rank": 5,
         "rank_up": 3, "rank_down": pd.NA, "timepoint": "tp1"},
        {"locus_tag": "G4", "gene_name": "down2", "log2fc": -0.5,
         "expression_status": "not_significant", "rank": 40,
         "rank_up": pd.NA, "rank_down": pd.NA, "timepoint": "tp1"},
        *_bg_genes(96),
    ])


class TestApplySignature:
    def test_all_present(self):
        """All 4 signature genes found in DE data."""
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        applied = apply_signature(sig, de)
        assert len(applied) == 4
        assert set(applied["locus_tag"]) == {"G1", "G2", "G3", "G4"}

    def test_concordance_correct(self):
        """Mixed DE: G1 concordant (+1), G3 reversed (-1), G4 not sig (0)."""
        sig = make_toy_signature()
        de = make_toy_de_mixed()
        applied = apply_signature(sig, de)
        g1 = applied[applied["locus_tag"] == "G1"].iloc[0]
        g3 = applied[applied["locus_tag"] == "G3"].iloc[0]
        g4 = applied[applied["locus_tag"] == "G4"].iloc[0]
        assert g1["concordance"] == 1
        assert g3["concordance"] == -1
        assert g4["concordance"] == 0

    def test_absent_gene(self):
        """G2 absent from dataset → row present with NaN DE values."""
        sig = make_toy_signature()
        de = make_toy_de_gene_absent()
        applied = apply_signature(sig, de)
        g2 = applied[applied["locus_tag"] == "G2"]
        assert len(g2) == 1
        assert pd.isna(g2.iloc[0]["log2fc"])


class TestRankScore:
    def test_all_concordant_hand_calculated(self):
        """All concordant, 100 total genes.
        G1: +1 * (1 - 1/100) = 0.990
        G2: +1 * (1 - 1/100) = 0.990
        G3: +1 * (1 - 2/100) = 0.980
        G4: +1 * (1 - 3/100) = 0.970
        mean = (0.990 + 0.990 + 0.980 + 0.970) / 4 = 0.9825
        """
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        applied = apply_signature(sig, de)
        score = rank_score(applied, total_genes=100)
        assert abs(score - 0.9825) < 0.001

    def test_mixed_hand_calculated(self):
        """G1 conc, G2 conc, G3 reversed, G4 not sig. 100 total genes.
        G1: +1 * (1 - 2/100) = +0.980
        G2: +1 * (1 - 3/100) = +0.970
        G3: -1 * (1 - 5/100) = -0.950
        G4:  0 * 0            =  0.000
        mean = (0.980 + 0.970 - 0.950 + 0) / 4 = 0.250
        """
        sig = make_toy_signature()
        de = make_toy_de_mixed()
        applied = apply_signature(sig, de)
        score = rank_score(applied, total_genes=100)
        assert abs(score - 0.250) < 0.001

    def test_absent_gene_excluded(self):
        """G2 absent → score computed over 3 genes, not 4.
        G1: +1 * (1 - 2/99) = +0.9798
        G3: +1 * (1 - 3/99) = +0.9697
        G4:  0 * 0           =  0.000
        mean over 3 = (0.9798 + 0.9697 + 0) / 3 = 0.6498
        """
        sig = make_toy_signature()
        de = make_toy_de_gene_absent()
        applied = apply_signature(sig, de)
        score = rank_score(applied, total_genes=99)
        expected = (0.9798 + 0.9697 + 0) / 3
        assert abs(score - expected) < 0.01

    def test_mixed_lower_than_concordant(self):
        """Mixed concordance → lower score than all concordant."""
        sig = make_toy_signature()
        applied_conc = apply_signature(sig, make_toy_de_all_concordant())
        applied_mixed = apply_signature(sig, make_toy_de_mixed())
        score_conc = rank_score(applied_conc, total_genes=100)
        score_mixed = rank_score(applied_mixed, total_genes=100)
        assert score_conc > score_mixed


class TestPermutationTest:
    def test_too_few_genes_returns_nan(self):
        """4 matched genes < 30 minimum → p-value is NaN, 0 permutations."""
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        result = permutation_test(sig, de, n_perms=200, seed=42)
        assert np.isnan(result["empirical_p"])
        assert result["n_permutations"] == 0
        assert not np.isnan(result["observed"])  # score still computed

    def test_large_signature_runs(self):
        """With enough genes (>30), permutation test runs and returns p-value."""
        # Build a 40-gene signature
        sig = pd.DataFrame([
            {"locus_tag": f"SIG{i}", "gene_name": f"sig{i}",
             "direction": "up" if i % 2 == 0 else "down"}
            for i in range(40)
        ])
        # Build DE with all 40 concordant + 60 background
        de_rows = []
        for i in range(40):
            d = "up" if i % 2 == 0 else "down"
            de_rows.append({
                "locus_tag": f"SIG{i}", "gene_name": f"sig{i}",
                "log2fc": 2.0 if d == "up" else -2.0,
                "expression_status": f"significant_{d}",
                "rank": i + 1,
                "rank_up": (i // 2 + 1) if d == "up" else pd.NA,
                "rank_down": (i // 2 + 1) if d == "down" else pd.NA,
                "timepoint": "tp1",
            })
        de_rows.extend(_bg_genes(60))
        de = pd.DataFrame(de_rows)
        result = permutation_test(sig, de, n_perms=100, seed=42)
        assert result["n_permutations"] == 100
        assert not np.isnan(result["empirical_p"])
        assert result["empirical_p"] < 0.1  # strong signal should be significant

    def test_returns_expected_keys(self):
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        result = permutation_test(sig, de, n_perms=50, seed=42)
        assert "observed" in result
        assert "empirical_p" in result
        assert "n_permutations" in result


class TestScoreWithSignificance:
    def test_returns_all_fields(self):
        """Wrapper returns score, p-value, and breakdown."""
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        result = score_with_significance(sig, de, n_perms=50)
        assert "rank_score" in result
        assert "empirical_p" in result
        assert "n_matched" in result
        assert "n_absent" in result
        assert "n_concordant" in result
        assert "n_reversed" in result
        assert "n_not_significant" in result
        assert "hit_rate" in result

    def test_all_concordant_counts(self):
        """4 matched, 0 absent, 4 concordant, 0 reversed, 0 not sig."""
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        result = score_with_significance(sig, de, n_perms=50)
        assert result["n_matched"] == 4
        assert result["n_absent"] == 0
        assert result["n_concordant"] == 4
        assert result["n_reversed"] == 0
        assert result["n_not_significant"] == 0
        assert result["hit_rate"] == 1.0

    def test_mixed_counts(self):
        """4 matched, 0 absent, 2 concordant, 1 reversed, 1 not sig."""
        sig = make_toy_signature()
        de = make_toy_de_mixed()
        result = score_with_significance(sig, de, n_perms=50)
        assert result["n_matched"] == 4
        assert result["n_concordant"] == 2
        assert result["n_reversed"] == 1
        assert result["n_not_significant"] == 1
        assert abs(result["hit_rate"] - 0.5) < 0.01

    def test_absent_gene_counts(self):
        """3 matched, 1 absent."""
        sig = make_toy_signature()
        de = make_toy_de_gene_absent()
        result = score_with_significance(sig, de, n_perms=50)
        assert result["n_matched"] == 3
        assert result["n_absent"] == 1
