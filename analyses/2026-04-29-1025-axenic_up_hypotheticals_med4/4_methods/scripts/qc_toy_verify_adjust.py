"""
Toy-data verification of the response-profile adjust+flag math (`_adjust_nitrogen_row`).

The adjustment subtracts the locked experiment's contribution from the nitrogen
row's clean-subtractable counts (experiments_*, timepoints_*) and keeps the
extremes (up_best_rank, up_max_log2fc, down_*) flagged in place.

Cases tested (hand-calculated):
  T1: Locked is sig_up; full has 8 N-experiments, 4 up across 14 timepoints
      Expected: other-than-locked = 7 experiments, 3 ups across 13 of 24 TPs
  T2: Locked is sig_up but only sub-sample of N experiments tested in cross-study
      Mirrors PMM1898's actual situation — 2 of 8 N tested, 1 up, 1 down
  T3: Locked is sig_up; cross-study has 0 other N experiments tested
      Edge case — adjustment produces all-zero "other" row
  T4: full_summary missing — function returns None
  T5: locked_summary missing — function returns full unchanged with no extremes flag
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dossier import _adjust_nitrogen_row


def case(name: str, full: dict | None, locked: dict | None, expected: dict | None) -> None:
    actual = _adjust_nitrogen_row(full, locked)
    if expected is None:
        assert actual is None, f"{name}: expected None, got {actual}"
        print(f"OK  {name}: returned None as expected")
        return
    if actual is None:
        raise AssertionError(f"{name}: expected dict, got None")
    other_actual = actual.get("other_experiments")
    other_expected = expected.get("other_experiments")
    if other_expected is not None:
        for k, v in other_expected.items():
            assert other_actual.get(k) == v, f"{name}: other_experiments.{k} expected {v}, got {other_actual.get(k)}"
    if "extremes_may_include_locked" in expected:
        assert actual.get("extremes_may_include_locked") == expected["extremes_may_include_locked"], (
            f"{name}: extremes_may_include_locked mismatch"
        )
    print(f"OK  {name}: other_experiments={other_actual} extremes_may_include_locked={actual.get('extremes_may_include_locked')}")


def make_summary(experiments_total=0, experiments_tested=0, experiments_up=0, experiments_down=0,
                 timepoints_total=0, timepoints_tested=0, timepoints_up=0, timepoints_down=0,
                 up_best_rank=None, up_max_log2fc=None, down_best_rank=None, down_max_log2fc=None) -> dict:
    return {
        "experiments_total": experiments_total,
        "experiments_tested": experiments_tested,
        "experiments_up": experiments_up,
        "experiments_down": experiments_down,
        "timepoints_total": timepoints_total,
        "timepoints_tested": timepoints_tested,
        "timepoints_up": timepoints_up,
        "timepoints_down": timepoints_down,
        "up_best_rank": up_best_rank,
        "up_max_log2fc": up_max_log2fc,
        "down_best_rank": down_best_rank,
        "down_max_log2fc": down_max_log2fc,
    }


def main() -> None:
    # T1: PMM0958-like — locked is 1 sig_up, cross-study has 8 N experiments, 4 up across 14 of 25 TPs
    full_T1 = make_summary(
        experiments_total=8, experiments_tested=8, experiments_up=4, experiments_down=0,
        timepoints_total=25, timepoints_tested=25, timepoints_up=14, timepoints_down=0,
        up_best_rank=1, up_max_log2fc=6.018,
    )
    locked_T1 = make_summary(
        experiments_total=1, experiments_tested=1, experiments_up=1, experiments_down=0,
        timepoints_total=1, timepoints_tested=1, timepoints_up=1, timepoints_down=0,
        up_best_rank=7, up_max_log2fc=5.740,
    )
    case("T1 PMM0958-like (8 N exps, 4 up, 14/25 TPs up)", full_T1, locked_T1, {
        "other_experiments": {
            "experiments_total": 7, "experiments_tested": 7,
            "experiments_up": 3, "experiments_down": 0,
            "timepoints_total": 24, "timepoints_tested": 24,
            "timepoints_up": 13, "timepoints_down": 0,
        },
        "extremes_may_include_locked": True,
    })

    # T2: PMM1898-like — only 2 of 8 N tested, 1 up, 1 down
    full_T2 = make_summary(
        experiments_total=8, experiments_tested=2, experiments_up=1, experiments_down=1,
        timepoints_total=25, timepoints_tested=6, timepoints_up=1, timepoints_down=1,
        up_best_rank=16, up_max_log2fc=4.680,
        down_best_rank=35, down_max_log2fc=-2.098,
    )
    locked_T2 = make_summary(
        experiments_total=1, experiments_tested=1, experiments_up=1, experiments_down=0,
        timepoints_total=1, timepoints_tested=1, timepoints_up=1, timepoints_down=0,
    )
    case("T2 PMM1898-like (2 of 8 tested, 1 up + 1 down)", full_T2, locked_T2, {
        "other_experiments": {
            "experiments_total": 7, "experiments_tested": 1,
            "experiments_up": 0, "experiments_down": 1,
            "timepoints_total": 24, "timepoints_tested": 5,
            "timepoints_up": 0, "timepoints_down": 1,
        },
        "extremes_may_include_locked": True,
    })

    # T3: only locked tested (no other N experiments tested)
    full_T3 = make_summary(
        experiments_total=8, experiments_tested=1, experiments_up=1, experiments_down=0,
        timepoints_total=25, timepoints_tested=1, timepoints_up=1, timepoints_down=0,
        up_best_rank=42, up_max_log2fc=2.5,
    )
    locked_T3 = make_summary(
        experiments_total=1, experiments_tested=1, experiments_up=1, experiments_down=0,
        timepoints_total=1, timepoints_tested=1, timepoints_up=1, timepoints_down=0,
    )
    case("T3 only locked tested (no cross-study N evidence)", full_T3, locked_T3, {
        "other_experiments": {
            "experiments_total": 7, "experiments_tested": 0,
            "experiments_up": 0, "experiments_down": 0,
            "timepoints_total": 24, "timepoints_tested": 0,
            "timepoints_up": 0, "timepoints_down": 0,
        },
        "extremes_may_include_locked": True,
    })

    # T4: full_summary missing
    case("T4 full missing", None, locked_T1, None)

    # T5: locked_summary missing
    res = _adjust_nitrogen_row(full_T1, None)
    assert res is not None
    assert res.get("extremes_may_include_locked") is False
    assert res.get("other_experiments") == full_T1
    print(f"OK  T5 locked missing: returned full unchanged, extremes_may_include_locked=False")

    print("\nAll toy-data verifications passed.")


if __name__ == "__main__":
    main()
