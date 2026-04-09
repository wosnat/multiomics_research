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
            out[col] = out[col].apply(
                lambda x: "|".join(sorted(x)) if isinstance(x, set) else x
            )
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
