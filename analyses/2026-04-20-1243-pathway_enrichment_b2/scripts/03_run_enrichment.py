"""Run pathway_enrichment per (organism, ontology, background_used).

Splits experiments by table_scope (all_detected_genes -> 'table_scope' background;
restricted scope -> 'organism' background) so Fisher's 2x2 is valid in each call.

Spec §5 Step 2. Two-stage pickle round-trip check:
  Stage 1 (before main loop): single-object probe — verify .explain() works
                              on a reloaded EnrichmentResult. Fail-fast.
  Stage 2 (after main loop): dict round-trip — pickle full results dict,
                             reload, verify .explain() on a non-empty value.
"""
from __future__ import annotations

import logging
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from multiomics_explorer import pathway_enrichment

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

RESTRICTED_SCOPES = {
    "significant_only",
    "significant_any_timepoint",
    "top_n",
    "filtered_subset",
}

# Locked in Step 1b (see ontology_selection.md, exploration/key_pathways.csv).
ONTOLOGIES: list[dict[str, object]] = [
    {"ontology": "cyanorak_role", "level": 1, "tree": None},
    {"ontology": "kegg", "level": 2, "tree": None},
]


def run_one(organism, ontology, level, experiment_ids, background, tree=None):
    kwargs = dict(
        organism=organism,
        experiment_ids=experiment_ids,
        ontology=ontology,
        level=level,
        background=background,
        direction="both",
        significant_only=True,
    )
    if tree is not None:
        kwargs["tree"] = tree
    return pathway_enrichment(**kwargs)


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step2.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step2")
    log.info("Starting 03_run_enrichment.py")

    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    if classified.empty:
        log.error("experiments_classified.csv is empty; cannot proceed.")
        return 1
    log.info(f"Loaded {len(classified)} rows "
             f"({classified['experiment_id'].nunique()} unique experiments)")

    if not ONTOLOGIES:
        log.error("ONTOLOGIES is empty — populate from Step 1b selection.")
        return 1

    # Dedupe to one row per experiment — experiments_to_dataframe returns one
    # row per (experiment × timepoint), so class/table_scope repeat across TPs.
    classified_unique = classified.drop_duplicates("experiment_id").copy()

    # ---- PICKLE STAGE-1 PROBE -------------------------------------------------
    # Run one small pathway_enrichment and verify single-object pickle round-trip
    # + .explain() BEFORE the main loop. Fail-fast per spec §5 Step 2 do.
    probe_row = classified_unique.iloc[0]
    probe_ont = ONTOLOGIES[0]
    log.info(f"Pickle probe: organism={probe_row['organism_name']!r}, "
             f"exp={probe_row['experiment_id']!r}, ontology={probe_ont['ontology']!r}")
    try:
        probe = run_one(
            organism=probe_row["organism_name"],
            ontology=probe_ont["ontology"],
            level=probe_ont["level"],
            experiment_ids=[probe_row["experiment_id"]],
            background="organism",  # safest default for single-experiment probe
            tree=probe_ont.get("tree"),
        )
    except Exception as e:
        log.error(f"Pickle probe enrichment call failed: {e}")
        return 1

    probe_pkl = ANALYSIS_DIR / "data" / "_probe_single.pkl"
    try:
        with open(probe_pkl, "wb") as f:
            pickle.dump(probe, f)
        with open(probe_pkl, "rb") as f:
            loaded_probe = pickle.load(f)
        if not loaded_probe.results.empty:
            first = loaded_probe.results.iloc[0]
            _ = loaded_probe.explain(first["cluster"], first["term_id"])
            log.info("Pickle stage 1 probe (single object): OK (.explain verified)")
        else:
            log.info("Pickle stage 1 probe (single object): OK (empty results, .explain skipped)")
    except Exception as e:
        log.error(f"Pickle stage 1 probe FAILED: {e} — see spec §8 risk 7 for fallbacks")
        probe_pkl.unlink(missing_ok=True)
        return 2
    finally:
        probe_pkl.unlink(missing_ok=True)

    # ---- MAIN LOOP -----------------------------------------------------------
    results: dict[tuple[str, str, str], object] = {}
    all_rows: list[pd.DataFrame] = []

    for org, org_group in classified_unique.groupby("organism_name"):
        is_detected = org_group["table_scope"] == "all_detected_genes"
        detected_ids = org_group.loc[is_detected, "experiment_id"].tolist()
        restricted_ids = org_group.loc[~is_detected, "experiment_id"].tolist()
        log.info(f"Organism {org!r}: {len(detected_ids)} all_detected, "
                 f"{len(restricted_ids)} restricted")

        for ont in ONTOLOGIES:
            for bg_label, exp_ids, bg_value in [
                ("table_scope", detected_ids, "table_scope"),
                ("organism", restricted_ids, "organism"),
            ]:
                if not exp_ids:
                    log.info(f"Skipping {org!r} × {ont['ontology']!r} × {bg_label}: no experiments")
                    continue
                try:
                    result = run_one(
                        organism=org,
                        ontology=ont["ontology"],
                        level=ont["level"],
                        experiment_ids=exp_ids,
                        background=bg_value,
                        tree=ont.get("tree"),
                    )
                except Exception as e:
                    log.error(
                        f"pathway_enrichment failed for ({org!r}, "
                        f"{ont['ontology']!r}, {bg_label}): {e}"
                    )
                    continue

                key = (org, ont["ontology"], bg_label)
                results[key] = result
                df = result.results.copy()
                if df.empty:
                    log.warning(f"{key}: empty results DataFrame; skipping concat")
                    continue
                df["organism"] = org
                df["ontology"] = ont["ontology"]
                df["background_used"] = bg_label
                all_rows.append(df)
                n_sig = (df["p_adjust"] < 0.05).sum() if "p_adjust" in df.columns else 0
                log.info(f"{key}: {len(df)} rows, {df['cluster'].nunique()} clusters, "
                         f"{n_sig} p_adjust<0.05")

    if not all_rows:
        log.error("No enrichment results produced across all (organism, ontology, bg) calls.")
        return 1

    combined = pd.concat(all_rows, ignore_index=True)
    combined.to_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv", index=False)
    log.info(f"Wrote {len(combined)} total rows to enrichment_all.csv "
             f"({combined['cluster'].nunique()} clusters, "
             f"{combined['term_id'].nunique()} unique terms)")

    # ---- PICKLE STAGE-2 DICT ROUND-TRIP --------------------------------------
    pkl_path = ANALYSIS_DIR / "data" / "enrichment_results.pkl"
    try:
        with open(pkl_path, "wb") as f:
            pickle.dump(results, f)
        with open(pkl_path, "rb") as f:
            loaded_dict = pickle.load(f)
    except Exception as e:
        log.error(f"Pickle stage 2 write/load FAILED: {e}")
        pkl_path.unlink(missing_ok=True)
        return 3

    verified = False
    for key, inst in loaded_dict.items():
        if not inst.results.empty:
            first = inst.results.iloc[0]
            try:
                _ = inst.explain(first["cluster"], first["term_id"])
                log.info(f"Pickle stage 2 (dict): OK — .explain verified on {key}")
                verified = True
                break
            except Exception as e:
                log.error(f"Pickle stage 2: .explain failed on loaded {key}: {e}")
                pkl_path.unlink(missing_ok=True)
                return 3
    if not verified:
        log.warning("Pickle stage 2: dict loaded but all values empty; .explain not exercised")

    # ---- DIAGNOSTIC HEATMAP (key-pathway × cluster, saturated at ±10) -------
    key_paths = pd.read_csv(ANALYSIS_DIR / "exploration" / "key_pathways.csv")
    key_term_ids = set(key_paths["term_id"])
    diag = combined[combined["term_id"].isin(key_term_ids)].copy()
    if diag.empty:
        log.warning("No key-pathway terms found in enrichment results — "
                    "diagnostic heatmap not written")
    else:
        pivot = diag.pivot_table(
            index=["ontology", "term_name"],
            columns=["organism", "experiment_id", "timepoint", "direction"],
            values="signed_score",
            aggfunc="first",
        )
        pivot = pivot.clip(lower=-10, upper=10)
        fig, ax = plt.subplots(
            figsize=(max(12, 0.28 * pivot.shape[1]), max(4, 0.35 * pivot.shape[0]))
        )
        sns.heatmap(
            pivot, center=0, cmap="RdBu_r", vmin=-10, vmax=10, ax=ax,
            cbar_kws={"label": "signed_score (capped at ±10)"},
        )
        ax.set_title("Step 2 diagnostic: key pathways × clusters (signed_score)")
        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=7)
        plt.setp(ax.get_yticklabels(), fontsize=8)
        plt.tight_layout()
        qc_path = ANALYSIS_DIR / "exploration" / "qc" / "step2_key_pathway_heatmap.png"
        plt.savefig(qc_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info(f"Wrote diagnostic heatmap to {qc_path} "
                 f"({pivot.shape[0]} key pathways × {pivot.shape[1]} clusters)")

    # ---- SUMMARY TO STDOUT ---------------------------------------------------
    print("\n=== Cluster / test / significance counts per (org, ontology, bg) ===")
    summary = (
        combined.groupby(["organism", "ontology", "background_used"])
        .agg(
            n_clusters=("cluster", "nunique"),
            n_tests=("term_id", "count"),
            n_significant=("p_adjust", lambda x: (x < 0.05).sum()),
        )
    )
    print(summary.to_string())

    print("\n=== Keys in enrichment_results.pkl ===")
    for k in sorted(results):
        n_rows = len(results[k].results) if hasattr(results[k], "results") else 0
        print(f"  {k}: {n_rows} rows")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
