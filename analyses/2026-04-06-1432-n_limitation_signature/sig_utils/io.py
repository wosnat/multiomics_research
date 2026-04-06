"""I/O helpers for signature analysis.

Load/save DE data and signature CSVs with consistent column types.
"""

import pandas as pd
from pathlib import Path

# Base path for data/ and results/ directories
ANALYSIS_DIR = Path(__file__).parent.parent
DATA_DIR = ANALYSIS_DIR / "data"
RESULTS_DIR = ANALYSIS_DIR / "results"


def load_de_csv(filename: str) -> pd.DataFrame:
    """Load a DE extract CSV from data/.

    Ensures consistent types for rank columns (nullable int)
    and expression_status (string).
    """
    path = DATA_DIR / filename
    df = pd.read_csv(path)

    # rank_up and rank_down are nullable int (None when not significant in that direction)
    for col in ["rank", "rank_up", "rank_down"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def load_signature_csv(filename: str) -> pd.DataFrame:
    """Load a signature definition CSV from data/.

    Expected columns: locus_tag, gene_name, direction, signature_type,
    plus rank columns.
    """
    path = DATA_DIR / filename
    df = pd.read_csv(path)

    for col in df.columns:
        if "rank" in col:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def save_scores_csv(df: pd.DataFrame, filename: str) -> Path:
    """Save a scores DataFrame to results/."""
    RESULTS_DIR.mkdir(exist_ok=True)
    path = RESULTS_DIR / filename
    df.to_csv(path, index=False)
    return path


def save_data_csv(df: pd.DataFrame, filename: str) -> Path:
    """Save a DataFrame to data/."""
    DATA_DIR.mkdir(exist_ok=True)
    path = DATA_DIR / filename
    df.to_csv(path, index=False)
    return path
