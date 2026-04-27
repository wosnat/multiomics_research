"""Direction-aware per-axis stress score.

Definition
----------
For a given (experiment × timepoint) and a given stress axis with gene set G
and per-gene direction d_g ∈ {+1, -1}, the axis stress score is a signed
Z-score against the genome-wide background log2FC distribution at that TP:

    signed_axis_lfc[g] = d_g * log2fc[g]    for g in G
    axis_mean          = mean over g of signed_axis_lfc[g]

    background         = log2fc for all genes quantified at this (experiment, TP)
                         EXCLUDING the axis genes
    background_mean    = mean(background)
    background_sd      = SD(background)  (population SD, ddof=0)

    axis_score         = (axis_mean - background_mean) / background_sd

Interpretation
--------------
`axis_score` is in units of "background SD" of log2FC.

  axis_score = +z (positive) → axis is engaged in the calibrated direction,
               more strongly than the typical gene at this TP. The bigger
               z is, the more the axis stands out against the genome-wide
               response.
  axis_score ≈ 0             → axis tracks the genome-wide background; no
               axis-specific signal beyond what the global response shows.
  axis_score = -z (negative) → axis moves in the OPPOSITE direction to its
               calibrated direction. Worth investigating: miscalibration,
               condition-specific reversal, or biology (e.g., transcription
               shutting down in axenic death-phase including stress-axis
               genes).

Direction calibration (the d_g values) comes from step 3 validation, not
from a textbook default. For n_stress all 5 positives are +1 (validated UP).
For photo, HLI proteins are +1 (UP) but psbA/psbD/ftsH2 are -1 (DOWN —
PSII disassembly under N-starvation in this dataset).

Worked example (hand-computed, also asserted in 4_methods/scripts/01_toy_verify.py)
-----------------------------------------------------------------------------------
Genes A, B, C in axis (all direction +1); genes D, E, F, G in background.

    log2fc = {A: 2.0, B: 1.5, C: 1.0, D: 0.1, E: 0.3, F: -0.2, G: 0.0}

    axis_mean       = mean(2.0, 1.5, 1.0)  = 1.5
    background      = [0.1, 0.3, -0.2, 0.0]
    background_mean = mean(...)            = 0.05
    background_var (ddof=0; divide by N=4) = ((0.1-0.05)^2 + (0.3-0.05)^2
                                            + (-0.2-0.05)^2 + (0-0.05)^2) / 4
                                          = 0.13 / 4
                                          = 0.0325
    background_sd                          = sqrt(0.0325) ≈ 0.180278

    axis_score = (1.5 - 0.05) / 0.180278   ≈ 8.043153

API
---
    axis_stress_score(
        de_df: pd.DataFrame,        # required cols: locus_tag, log2fc
        axis_genes: list[str],      # locus tags in the axis
        direction: dict[str, int],  # {locus_tag: +1 | -1}
    ) -> dict

returns

    {
      "axis_score":     float,
      "axis_mean":      float,    # mean of signed log2fc over axis genes
      "background_mean": float,
      "background_sd":   float,
      "n_axis":         int,      # axis genes with non-NaN log2fc
      "n_background":   int,      # background genes with non-NaN log2fc
      "axis_genes_with_data":   list[str],
      "axis_genes_missing_data": list[str],
    }
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd


def axis_stress_score(
    de_df: pd.DataFrame,
    axis_genes: Iterable[str],
    direction: dict[str, int],
) -> dict:
    """Direction-aware signed-Z score for a stress axis at a single (experiment, TP).

    See module docstring for definition and worked example.

    Parameters
    ----------
    de_df : DataFrame with at minimum columns ``locus_tag`` and ``log2fc``.
            Should be filtered to a single (experiment, timepoint) before
            being passed in (this function does not group).
    axis_genes : iterable of locus_tag strings — the axis gene set.
    direction : mapping locus_tag -> +1 or -1. Must include every gene in
            ``axis_genes`` (raises KeyError otherwise — fail loudly rather
            than silently default to +1).

    Returns
    -------
    dict — see module docstring for fields. NaN log2fc values are excluded.
    Returns axis_score=NaN if there are <1 axis gene with data or <2
    background genes with data.

    Notes
    -----
    SD uses ddof=0 (population SD), matching the worked example. Switch to
    ddof=1 (sample SD) only if needed for inferential symmetry.
    """
    axis_genes = list(axis_genes)
    missing_dir = [g for g in axis_genes if g not in direction]
    if missing_dir:
        raise KeyError(
            f"direction dict missing entries for axis genes: {missing_dir}. "
            "Direction must be explicitly specified per gene; do not default to +1."
        )

    if "locus_tag" not in de_df.columns or "log2fc" not in de_df.columns:
        raise ValueError(
            f"de_df must have 'locus_tag' and 'log2fc' columns; got {list(de_df.columns)}"
        )

    axis_set = set(axis_genes)

    # Axis: signed log2fc, drop NaN
    axis_rows = de_df[de_df["locus_tag"].isin(axis_set)][["locus_tag", "log2fc"]].copy()
    axis_rows["log2fc"] = pd.to_numeric(axis_rows["log2fc"], errors="coerce")
    axis_rows = axis_rows.dropna(subset=["log2fc"])
    axis_rows["signed_lfc"] = axis_rows["locus_tag"].map(direction).astype(float) * axis_rows["log2fc"]

    axis_genes_with_data = sorted(axis_rows["locus_tag"].tolist())
    axis_genes_missing_data = sorted(set(axis_genes) - set(axis_genes_with_data))

    # Background: all non-axis genes, drop NaN
    bg_rows = de_df[~de_df["locus_tag"].isin(axis_set)][["locus_tag", "log2fc"]].copy()
    bg_rows["log2fc"] = pd.to_numeric(bg_rows["log2fc"], errors="coerce")
    bg_rows = bg_rows.dropna(subset=["log2fc"])

    n_axis = len(axis_rows)
    n_background = len(bg_rows)

    if n_axis < 1 or n_background < 2:
        return {
            "axis_score": float("nan"),
            "axis_mean": float("nan"),
            "background_mean": float("nan"),
            "background_sd": float("nan"),
            "n_axis": n_axis,
            "n_background": n_background,
            "axis_genes_with_data": axis_genes_with_data,
            "axis_genes_missing_data": axis_genes_missing_data,
        }

    axis_mean = float(axis_rows["signed_lfc"].mean())
    background_mean = float(bg_rows["log2fc"].mean())
    background_sd = float(np.std(bg_rows["log2fc"].values, ddof=0))

    if background_sd == 0:
        return {
            "axis_score": float("nan"),
            "axis_mean": axis_mean,
            "background_mean": background_mean,
            "background_sd": 0.0,
            "n_axis": n_axis,
            "n_background": n_background,
            "axis_genes_with_data": axis_genes_with_data,
            "axis_genes_missing_data": axis_genes_missing_data,
        }

    axis_score = (axis_mean - background_mean) / background_sd

    return {
        "axis_score": float(axis_score),
        "axis_mean": axis_mean,
        "background_mean": background_mean,
        "background_sd": background_sd,
        "n_axis": n_axis,
        "n_background": n_background,
        "axis_genes_with_data": axis_genes_with_data,
        "axis_genes_missing_data": axis_genes_missing_data,
    }
