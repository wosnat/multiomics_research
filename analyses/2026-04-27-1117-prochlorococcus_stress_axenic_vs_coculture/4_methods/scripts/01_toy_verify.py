"""Toy-data verification for axis_stress_score.

Hand-computed expected values:
  axis genes A, B, C all direction +1; log2fc 2.0, 1.5, 1.0 → signed_mean = 1.5
  background D, E, F, G with log2fc 0.1, 0.3, -0.2, 0.0
    background_mean = 0.05
    background_var (population, ddof=0) = ((0.1-0.05)^2 + (0.3-0.05)^2 + (-0.2-0.05)^2 + (0-0.05)^2) / 4
                                       = (0.0025 + 0.0625 + 0.0625 + 0.0025) / 4
                                       = 0.13 / 4 = 0.0325
    background_sd = sqrt(0.0325) ≈ 0.180278
  axis_score = (1.5 - 0.05) / 0.180278 ≈ 8.0432

(Earlier docstring said sqrt(0.13/3) = 0.20817, but the implementation uses
ddof=0 — population variance — so dividing by N (=4), not N-1. The number
in the module docstring matches ddof=0 which is what the implementation
does. This script is the source of truth for the expected value.)

Also tests:
  - direction = -1 case: axis genes "down under stress" should give the same
    score if their log2fc is multiplied by -1 (the function's job)
  - missing direction key raises KeyError
  - NaN log2fc is excluded
  - axis with 1 gene works (n_axis>=1 minimum)

Run:
  uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/4_methods/scripts/01_toy_verify.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Import the module under test from its analysis-local location
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from stress_score import axis_stress_score  # noqa: E402


def case_basic_positive_direction() -> None:
    de = pd.DataFrame({
        "locus_tag": ["A", "B", "C", "D", "E", "F", "G"],
        "log2fc":    [2.0, 1.5, 1.0, 0.1, 0.3, -0.2, 0.0],
    })
    direction = {"A": 1, "B": 1, "C": 1}
    res = axis_stress_score(de, axis_genes=["A", "B", "C"], direction=direction)

    expected_axis_mean = 1.5
    expected_bg_mean = 0.05
    expected_bg_sd = math.sqrt(0.0325)
    expected_score = (expected_axis_mean - expected_bg_mean) / expected_bg_sd

    print("[case 1] basic, all positive direction")
    print(f"  axis_mean        = {res['axis_mean']:.6f}  (expected {expected_axis_mean:.6f})")
    print(f"  background_mean  = {res['background_mean']:.6f}  (expected {expected_bg_mean:.6f})")
    print(f"  background_sd    = {res['background_sd']:.6f}  (expected {expected_bg_sd:.6f})")
    print(f"  axis_score       = {res['axis_score']:.6f}  (expected {expected_score:.6f})")
    assert math.isclose(res["axis_mean"], expected_axis_mean, rel_tol=1e-12)
    assert math.isclose(res["background_mean"], expected_bg_mean, rel_tol=1e-12)
    assert math.isclose(res["background_sd"], expected_bg_sd, rel_tol=1e-12)
    assert math.isclose(res["axis_score"], expected_score, rel_tol=1e-12)
    assert res["n_axis"] == 3 and res["n_background"] == 4
    print("  OK\n")


def case_negative_direction_symmetry() -> None:
    """If genes A, B, C go DOWN under stress and we calibrate their direction = -1,
    then negating their log2fc gives the same axis_score as the all-+1 case."""
    de = pd.DataFrame({
        "locus_tag": ["A", "B", "C", "D", "E", "F", "G"],
        "log2fc":    [-2.0, -1.5, -1.0, 0.1, 0.3, -0.2, 0.0],
    })
    direction = {"A": -1, "B": -1, "C": -1}
    res = axis_stress_score(de, axis_genes=["A", "B", "C"], direction=direction)
    print("[case 2] all negative direction with negative log2fc → same axis_score")
    print(f"  axis_score = {res['axis_score']:.6f}")
    expected = (1.5 - 0.05) / math.sqrt(0.0325)
    assert math.isclose(res["axis_score"], expected, rel_tol=1e-12), (
        f"expected {expected}, got {res['axis_score']}"
    )
    print("  OK\n")


def case_mixed_direction_bidirectional_axis() -> None:
    """For bidirectional axes (e.g., photo: HLI up, psbA down), each gene gets
    its own calibrated direction."""
    de = pd.DataFrame({
        "locus_tag": ["HLI1", "psbA",  "psbD",  "BG1", "BG2", "BG3", "BG4"],
        "log2fc":    [4.0,    -3.0,    -2.5,    0.1,   -0.05, 0.2,   -0.15],
    })
    # All three "stress engaged" — direction calibration converts to all-positive contribution
    direction = {"HLI1": 1, "psbA": -1, "psbD": -1}
    res = axis_stress_score(de, axis_genes=["HLI1", "psbA", "psbD"], direction=direction)

    # signed_lfc = +4.0 (HLI1), +3.0 (psbA), +2.5 (psbD) → mean = 3.1666...
    expected_axis_mean = (4.0 + 3.0 + 2.5) / 3
    bg = np.array([0.1, -0.05, 0.2, -0.15])
    expected_bg_mean = float(bg.mean())
    expected_bg_sd = float(np.std(bg, ddof=0))
    expected_score = (expected_axis_mean - expected_bg_mean) / expected_bg_sd

    print("[case 3] bidirectional axis (signs differ per gene)")
    print(f"  axis_mean   = {res['axis_mean']:.6f}  (expected {expected_axis_mean:.6f})")
    print(f"  axis_score  = {res['axis_score']:.6f}  (expected {expected_score:.6f})")
    assert math.isclose(res["axis_mean"], expected_axis_mean, rel_tol=1e-10)
    assert math.isclose(res["axis_score"], expected_score, rel_tol=1e-10)
    print("  OK\n")


def case_missing_direction_raises() -> None:
    de = pd.DataFrame({"locus_tag": ["A", "B"], "log2fc": [1.0, 0.5]})
    print("[case 4] missing direction entry must raise KeyError")
    try:
        axis_stress_score(de, axis_genes=["A", "B"], direction={"A": 1})
    except KeyError as e:
        print(f"  raised KeyError: {e}")
        print("  OK\n")
        return
    raise AssertionError("expected KeyError, did not raise")


def case_nan_excluded() -> None:
    """NaN log2fc rows should be excluded from both axis and background."""
    de = pd.DataFrame({
        "locus_tag": ["A", "B", "C", "D", "E"],
        "log2fc":    [2.0, np.nan, 1.0, 0.1, np.nan],
    })
    direction = {"A": 1, "B": 1, "C": 1}
    res = axis_stress_score(de, axis_genes=["A", "B", "C"], direction=direction)
    # Effective: axis A, C with log2fc 2.0, 1.0 → mean 1.5
    #            background D with log2fc 0.1 — only one row, n_background=1 < 2 → score=NaN
    print("[case 5] NaN log2fc excluded; n_background=1 returns NaN score")
    assert res["n_axis"] == 2
    assert res["n_background"] == 1
    assert math.isnan(res["axis_score"])
    print(f"  n_axis={res['n_axis']}  n_background={res['n_background']}  axis_score={res['axis_score']}")
    print("  OK\n")


def case_single_gene_axis() -> None:
    """An axis can be 1 gene — the function should not require >= 2."""
    de = pd.DataFrame({
        "locus_tag": ["A", "B", "C", "D", "E"],
        "log2fc":    [3.0, 0.0, 0.1, -0.1, 0.05],
    })
    direction = {"A": 1}
    res = axis_stress_score(de, axis_genes=["A"], direction=direction)
    print("[case 6] single-gene axis")
    print(f"  axis_score={res['axis_score']:.4f}  n_axis={res['n_axis']}  n_background={res['n_background']}")
    assert res["n_axis"] == 1
    assert res["n_background"] == 4
    assert not math.isnan(res["axis_score"])
    print("  OK\n")


def main() -> None:
    case_basic_positive_direction()
    case_negative_direction_symmetry()
    case_mixed_direction_bidirectional_axis()
    case_missing_direction_raises()
    case_nan_excluded()
    case_single_gene_axis()
    print("All toy cases passed.")


if __name__ == "__main__":
    main()
