"""Derive MED4-only reference N-limitation signature.

I/O + logging runner. Pure derivation logic lives in signature.py.

Temporal filter: keep clusters with timepoint_hours > TIMEPOINT_HOURS_CUTOFF
(strict greater-than) per decisions.md D1. Excludes 0h AND 3h.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from signature import (
    CORE_SUPPORT,
    FALLBACK_SUPPORT,
    TIMEPOINT_HOURS_CUTOFF,
    derive_for_ontology,
)

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


MED4_ORGANISM_SUBSTRING = "MED4"  # organism_name filter (substring, case-insensitive)


def _append_fallback_to_decisions(ontology: str, core_size: int, analysis_dir: Path) -> None:
    """Append a fallback-rule entry to decisions.md (spec §5 Step 3)."""
    decisions = analysis_dir / "decisions.md"
    from datetime import datetime
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    section = (
        f"\n## D-auto — Signature derivation fallback ({ontology})\n\n"
        f"**Applied:** {stamp}\n\n"
        f"Ontology `{ontology}` produced {core_size} signature pathways under the "
        f"core rule (≥{CORE_SUPPORT} same-direction R clusters). Because {core_size} < 5, "
        f"the fallback threshold (≥{FALLBACK_SUPPORT}) was applied. The resulting "
        f"signature is weaker — document in `caveats.md`.\n"
    )
    with decisions.open("a") as f:
        f.write(section)


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step3.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step3")

    enrichment = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")

    # Restrict to MED4 R clusters — match MED4 by substring so we're robust to
    # display-name variants ("Prochlorococcus MED4", "MED4", case drift, etc.).
    organism_mask = classified["organism_name"].fillna("").str.contains(
        MED4_ORGANISM_SUBSTRING, case=False, regex=False
    )
    r_exp_ids = set(
        classified.loc[(classified["class"] == "R") & organism_mask, "experiment_id"]
    )
    if not r_exp_ids:
        log.error(
            f"No MED4 R experiments found (filter: class=R AND organism_name contains "
            f"{MED4_ORGANISM_SUBSTRING!r}). Check experiments_classified.csv classifications "
            f"and organism_name spelling."
        )
        return 1

    r_clusters = enrichment[enrichment["experiment_id"].isin(r_exp_ids)].copy()
    log.info(
        f"MED4 R experiments: {len(r_exp_ids)}, "
        f"clusters: {r_clusters['cluster'].nunique()}"
    )
    if r_clusters.empty:
        log.error(
            f"MED4 R experiments ({sorted(r_exp_ids)}) classified but produced zero "
            f"enrichment rows. Verify enrichment_all.csv covers these experiments."
        )
        return 1

    # Temporal filter: keep clusters with timepoint_hours > cutoff (strict >,
    # per decisions.md D1 — excludes 0h AND 3h). NaN hours are retained
    # (no temporal metadata → cannot filter; fall back to keeping).
    if "timepoint_hours" in r_clusters.columns:
        hours_col = r_clusters["timepoint_hours"]
        mask_keep = hours_col.isna() | (hours_col > TIMEPOINT_HOURS_CUTOFF)
        dropped_early = int((~mask_keep).sum())
        dropped_clusters = sorted(
            r_clusters.loc[~mask_keep, "cluster"].unique().tolist()
        )
        r_clusters = r_clusters[mask_keep].copy()
        log.info(
            f"Temporal filter (hours > {TIMEPOINT_HOURS_CUTOFF}, D1): "
            f"dropped {dropped_early} early-cluster rows "
            f"from {len(dropped_clusters)} clusters: {dropped_clusters}"
        )
    else:
        log.warning(
            "timepoint_hours column missing from enrichment_all.csv; temporal filter skipped. "
            "Verify 03_run_enrichment.py passes cluster metadata through."
        )

    if r_clusters.empty:
        log.error("After temporal filter, r_clusters is empty — cannot derive signature.")
        return 1

    # Per-ontology derivation.
    all_sig = []
    all_dropped = []
    fallback_ontologies: list[tuple[str, int]] = []  # (ontology, core_size)
    for ontology, ont_df in r_clusters.groupby("ontology"):
        sig_df, drop_df = derive_for_ontology(ont_df, CORE_SUPPORT)
        rule = "core"
        if len(sig_df) < 5:
            core_size = len(sig_df)
            log.warning(
                f"Ontology {ontology}: {core_size} < 5 signature pathways "
                f"under core rule. Applying fallback (>={FALLBACK_SUPPORT})."
            )
            sig_df, drop_df = derive_for_ontology(ont_df, FALLBACK_SUPPORT)
            rule = "fallback"
            fallback_ontologies.append((ontology, core_size))
        sig_df["ontology"] = ontology
        sig_df["rule_applied"] = rule
        drop_df["ontology"] = ontology
        all_sig.append(sig_df)
        all_dropped.append(drop_df)
        log.info(
            f"Ontology {ontology}: signature size = {len(sig_df)} ({rule}), "
            f"dropped = {len(drop_df)}"
        )

    signature_df = pd.concat(all_sig, ignore_index=True) if all_sig else pd.DataFrame()
    dropped = pd.concat(all_dropped, ignore_index=True) if all_dropped else pd.DataFrame()

    if signature_df.empty:
        log.error(
            "Derived signature is empty across all ontologies — cannot proceed. "
            "Either R set is too small, or no pathway meets the support threshold."
        )
        return 1

    signature_df.to_csv(ANALYSIS_DIR / "data" / "reference_signature.csv", index=False)
    dropped.to_csv(ANALYSIS_DIR / "data" / "signature_dropped.csv", index=False)
    log.info(f"Wrote {len(signature_df)} signature, {len(dropped)} dropped")

    # Record any fallback applications to decisions.md (spec §5 Step 3).
    for ont, core_size in fallback_ontologies:
        _append_fallback_to_decisions(ont, core_size, ANALYSIS_DIR)
        log.info(f"Appended fallback-rule entry for ontology={ont} to decisions.md")

    print("\n=== Signature sizes per ontology ===")
    print(signature_df.groupby(["ontology", "rule_applied", "direction"]).size().to_string())
    print("\n=== Dropped sizes per ontology ===")
    if not dropped.empty:
        print(dropped.groupby(["ontology", "drop_reason"]).size().to_string())
    else:
        print("(none)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
