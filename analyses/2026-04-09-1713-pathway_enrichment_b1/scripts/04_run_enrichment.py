"""Step 4: Run Fisher's exact enrichment across all v2 experiments.

Loads pathway definitions, runs enrichment per experiment × timepoint,
concatenates results, and saves full and significant subsets.

Usage:
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/04_run_enrichment.py
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/04_run_enrichment.py \\
        --table-scope all_detected_genes
    uv run analyses/2026-04-09-1713-pathway_enrichment_b1/scripts/04_run_enrichment.py --explore
"""

import argparse
import datetime
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from enrich_utils.enrichment import run_enrichment_all_timepoints
from enrich_utils.io import load_de, save_enrichment_results

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ANALYSIS_DIR / "data"
RESULTS_DIR = ANALYSIS_DIR / "results"
LOGS_DIR = ANALYSIS_DIR / "logs"
V2_DATA_DIR = Path(__file__).resolve().parents[3] / "2026-04-08-1038-n_limitation_signature_v2" / "data"

# Keywords for marker pathway --explore traces
MARKER_PATHWAY_KEYWORDS = ["nitrogen", "photosynthesis", "transport", "ribosomal", "ribosome"]

PADJ_THRESHOLD = 0.05

VALID_TABLE_SCOPES = [
    "all_detected_genes",
    "filtered_subset",
    "significant_any_timepoint",
    "significant_only",
]


def load_pathway_defs(path: Path):
    """Load pathway definitions CSV, converting pipe-delimited locus_tags back to sets."""
    import pandas as pd

    df = pd.read_csv(path)
    df["locus_tags"] = df["locus_tags"].apply(
        lambda x: set(x.split("|")) if isinstance(x, str) and x else set()
    )
    return df


def main(table_scope: str, explore: bool = False) -> None:
    import pandas as pd

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    log_path = LOGS_DIR / "04_run_enrichment.log"
    log_lines = []

    def log(msg: str) -> None:
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        log_lines.append(line)

    log(f"=== Step 4: Run Enrichment (table_scope={table_scope}) ===")

    # ------------------------------------------------------------------
    # 1. Load pathway definitions
    # ------------------------------------------------------------------
    defs_path = DATA_DIR / "pathway_definitions.csv"
    if not defs_path.exists():
        raise FileNotFoundError(
            f"Pathway definitions not found: {defs_path}\n"
            f"Run 03_define_pathways.py first."
        )
    pathway_defs = load_pathway_defs(defs_path)
    log(f"Loaded {len(pathway_defs)} pathway definitions")

    # ------------------------------------------------------------------
    # 2. Run enrichment for each DE file
    # ------------------------------------------------------------------
    de_files = sorted(V2_DATA_DIR.glob("de_*.csv"))
    if not de_files:
        raise FileNotFoundError(f"No de_*.csv files found in {V2_DATA_DIR}")
    log(f"Found {len(de_files)} DE files")

    all_results = []
    total_tests = 0

    for de_path in de_files:
        exp_name = de_path.stem  # e.g. de_weissberg_rnaseq_axenic
        log(f"  Processing {exp_name}...")

        de_df = load_de(de_path)

        # Check for required columns
        required_cols = {"locus_tag", "expression_status"}
        missing = required_cols - set(de_df.columns)
        if missing:
            log(f"    WARNING: missing columns {missing} — skipping")
            continue

        if "timepoint" not in de_df.columns:
            log(f"    WARNING: no 'timepoint' column — skipping (run_enrichment_all_timepoints requires it)")
            continue

        result_df = run_enrichment_all_timepoints(de_df, pathway_defs, table_scope)
        result_df["experiment"] = exp_name
        all_results.append(result_df)

        n_tests = len(result_df)
        n_sig = (result_df["padj"] < PADJ_THRESHOLD).sum() if "padj" in result_df.columns else 0
        total_tests += n_tests
        log(f"    tests={n_tests}, significant (padj<{PADJ_THRESHOLD})={n_sig}")

        # Per-timepoint breakdown
        if "timepoint" in result_df.columns:
            for tp, tp_df in result_df.groupby("timepoint"):
                n_tp_sig = (tp_df["padj"] < PADJ_THRESHOLD).sum()
                log(f"      timepoint={tp}: tests={len(tp_df)}, significant={n_tp_sig}")

    # ------------------------------------------------------------------
    # 3. Concatenate and save
    # ------------------------------------------------------------------
    if not all_results:
        log("ERROR: No results to save. Check DE files and pathway definitions.")
        return

    combined = pd.concat(all_results, ignore_index=True)
    log(f"Total tests: {total_tests}")
    n_sig_total = (combined["padj"] < PADJ_THRESHOLD).sum()
    log(f"Total significant (padj<{PADJ_THRESHOLD}): {n_sig_total}")

    all_path = RESULTS_DIR / "enrichment_all.csv"
    save_enrichment_results(combined, all_path)
    log(f"Saved enrichment_all.csv ({len(combined)} rows)")

    sig_df = combined[combined["padj"] < PADJ_THRESHOLD].copy()
    sig_path = RESULTS_DIR / "enrichment_significant.csv"
    save_enrichment_results(sig_df, sig_path)
    log(f"Saved enrichment_significant.csv ({len(sig_df)} rows)")

    # ------------------------------------------------------------------
    # 4. Explore: trace marker pathways
    # ------------------------------------------------------------------
    if explore and len(sig_df) > 0:
        log("--- Marker pathway traces ---")
        for keyword in MARKER_PATHWAY_KEYWORDS:
            mask = (
                sig_df["pathway_name"].str.lower().str.contains(keyword, na=False)
                | sig_df["pathway_id"].str.lower().str.contains(keyword, na=False)
            )
            hits = sig_df[mask]
            if len(hits) == 0:
                log(f"  '{keyword}': no significant hits")
            else:
                log(f"  '{keyword}': {len(hits)} significant rows")
                for _, row in hits.iterrows():
                    log(
                        f"    {row.get('experiment', '?')} tp={row.get('timepoint', '?')} "
                        f"dir={row['direction']} {row['pathway_id']} | {row['pathway_name']} "
                        f"padj={row['padj']:.3e} OR={row['odds_ratio']:.2f}"
                    )

    # Write log file
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Fisher's exact enrichment for all v2 experiments"
    )
    parser.add_argument(
        "--table-scope",
        default="all_detected_genes",
        choices=VALID_TABLE_SCOPES,
        dest="table_scope",
        help=(
            "Gene universe definition for enrichment tests. "
            "Default: all_detected_genes. "
            "Options: " + ", ".join(VALID_TABLE_SCOPES)
        ),
    )
    parser.add_argument(
        "--explore",
        action="store_true",
        help="Trace marker pathways (nitrogen, photosynthesis, transport, ribosomal)",
    )
    args = parser.parse_args()
    main(table_scope=args.table_scope, explore=args.explore)
