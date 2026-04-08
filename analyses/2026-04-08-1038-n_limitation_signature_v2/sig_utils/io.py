"""I/O helpers for signature analysis.

Load/save DE data, signatures, and score tables with consistent column types.
"""

import pandas as pd
from pathlib import Path

# Base path for data/ and results/ directories
ANALYSIS_DIR = Path(__file__).parent.parent
DATA_DIR = ANALYSIS_DIR / "data"
RESULTS_DIR = ANALYSIS_DIR / "results"
LOGS_DIR = ANALYSIS_DIR / "logs"


def load_de_csv(path: str | Path) -> pd.DataFrame:
    """Load a DE extract CSV.

    Ensures consistent types for rank columns (nullable int)
    and expression_status (string).

    Args:
        path: Path to CSV file (absolute or relative to DATA_DIR).
    """
    path = Path(path)
    if not path.is_absolute():
        path = DATA_DIR / path

    df = pd.read_csv(path)

    for col in ["rank", "rank_up", "rank_down"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def load_signature_csv(path: str | Path) -> pd.DataFrame:
    """Load a signature definition CSV.

    Expected columns: locus_tag, gene_name, direction, signature_type,
    plus rank columns.

    Args:
        path: Path to CSV file (absolute or relative to DATA_DIR).
    """
    path = Path(path)
    if not path.is_absolute():
        path = DATA_DIR / path

    df = pd.read_csv(path)

    for col in df.columns:
        if "rank" in col:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def save_csv(df: pd.DataFrame, path: str | Path, subdir: str = "data") -> Path:
    """Save a DataFrame to a CSV in the analysis directory.

    Args:
        df: DataFrame to save.
        path: Filename or path. If relative, placed under subdir.
        subdir: "data" or "results". Ignored if path is absolute.

    Returns:
        The absolute path where the file was saved.
    """
    path = Path(path)
    if not path.is_absolute():
        base = DATA_DIR if subdir == "data" else RESULTS_DIR
        base.mkdir(exist_ok=True)
        path = base / path

    df.to_csv(path, index=False)
    return path


def write_log(content: str, filename: str) -> Path:
    """Write diagnostic log content to logs/.

    Args:
        content: Log text to write.
        filename: Log filename (e.g., '01_discover_experiments.log').

    Returns:
        The absolute path where the log was saved.
    """
    LOGS_DIR.mkdir(exist_ok=True)
    path = LOGS_DIR / filename
    path.write_text(content)
    return path
