"""Step 5: score all (axis × experiment × TP) cells using two panel kinds.

Two panel kinds per axis (where available):
  - "positive": the step-3 validated positive-control panel
                (3-5 genes per axis, per-gene direction from validation)
  - "cyanorak": the cyanorak hand-curated functional role
                (~25-30 genes per axis, textbook per-axis direction default)

Cell-death has no cyanorak panel because cyanorak.role:D.3 is empty for MED4.

Inputs (read):
  - 3_analysis_framing/data/control_panel.csv           (positive panel + step-3 evidence)
  - 2_kg_selection/data/experiments_pro_med4.csv        (4 trajectory + 1 single-point experiment ids)

Outputs (write):
  - data/all_axes_scores.csv          one row per (axis, panel_kind, experiment, TP)
  - data/genome_de_<short_id>.csv     genome-wide DE per experiment (4 files)
  - data/panel_definitions.csv        the gene lists + per-gene direction used (for transparency)
  - data/01_score_all_axes.log

Run:
  uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/5_analyze/scripts/01_score_all_axes.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

from multiomics_explorer import differential_expression_by_gene, genes_by_ontology
from multiomics_explorer.analysis import to_dataframe

# Local import: stress_score module from the methods step
ANALYSIS_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ANALYSIS_DIR / "4_methods"))
from stress_score import axis_stress_score  # noqa: E402

CONTROL_CSV = ANALYSIS_DIR / "3_analysis_framing" / "data" / "control_panel.csv"
EXP_CSV = ANALYSIS_DIR / "2_kg_selection" / "data" / "experiments_pro_med4.csv"

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = OUT_DIR / "01_score_all_axes.log"


# Direction calibration for the positive panel (from step 3 validation).
# Genes validated DOWN under stress get -1; rest +1. Wrong-handed markers
# (isiB) are calibrated based on observed direction in the data, per the
# step 3 decision to keep them in the panel transparently.
POSITIVE_DIRECTION = {
    # n_stress — all validated UP
    "PMM0246": +1,  # ntcA
    "PMM0263": +1,  # amt1
    "PMM0920": +1,  # glnA
    "PMM0970": +1,  # urtA
    "PMM1463": +1,  # glnB
    # oxidative — prxQ DOWN; rest flat-ish but kept canonical +1
    "PMM0079": -1,  # prxQ — validated DOWN
    "PMM0567": +1,  # gor
    "PMM0856": +1,  # ahpC
    "PMM1061": +1,  # trxA
    "PMM1294": +1,  # sodN
    # proteotoxic — clpX DOWN; rest +1
    "PMM1436": +1,  # groL1
    "PMM1432": +1,  # dnaK1
    "PMM0901": +1,  # htpG
    "PMM1657": -1,  # clpX — validated DOWN
    "PMM1437": +1,  # groES
    # photo — bidirectional
    "PMM0223": -1,  # psbA — DOWN, PSII disassembly
    "PMM1157": -1,  # psbD — DOWN
    "PMM0743": -1,  # ftsH2 — DOWN
    "PMM1404": +1,  # hli — strongly UP
    "PMM0064": +1,  # hli — UP weak
    # cell_death — mostly DOWN (the canonical stress-engaged direction)
    "PMM0191": +1,  # spoT — mixed; canonical +1
    "PMM0400": -1,  # lrtA — validated DOWN under stress
    "PMM1171": -1,  # isiB — validated DOWN (wrong-handed for N-stress per S5; kept with calibrated direction)
}

# Cyanorak panel: textbook per-axis default direction
CYANORAK_AXIS_DIRECTION = {
    "n_stress": +1,       # nitrogen scavenging genes UP under N starvation
    "oxidative": +1,      # oxidative defense UP under oxidative stress
    "proteotoxic": +1,    # chaperones UP under proteotoxic stress
    "photo": -1,          # cyanorak.role:J.8 = PSII; PSII genes DOWN under stress (disassembly)
}

# Cyanorak roles per axis. cell_death has no cyanorak panel because D.3 is empty for MED4.
CYANORAK_ROLES = {
    "n_stress":   ["cyanorak.role:D.1.3", "cyanorak.role:E.4"],
    "oxidative":  ["cyanorak.role:D.1.4"],
    "proteotoxic": ["cyanorak.role:L.3"],
    "photo":      ["cyanorak.role:J.8"],
}


def log(msg: str, *, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def short_exp_id(experiment_id: str) -> str:
    """Make a short id like 'rnaseq_axenic' for filenames."""
    suffix = experiment_id.split("med4_")[-1]
    return suffix


def main() -> None:
    with LOG_PATH.open("w") as fh:
        # Load positive panel
        controls = pd.read_csv(CONTROL_CSV)
        positives = controls[controls["role"] == "positive"].copy()
        log(f"Positive panel (step 3): {len(positives)} genes across "
            f"{positives['axis'].nunique()} axes", fh=fh)

        # Build positive-panel gene sets per axis
        positive_sets = {axis: list(g["locus_tag"]) for axis, g in positives.groupby("axis")}
        for axis, genes in positive_sets.items():
            log(f"  positive[{axis:<12}] n={len(genes)}: {genes}", fh=fh)

        # Build cyanorak gene sets per axis (skip cell_death — no cyanorak set)
        cyanorak_sets: dict[str, list[str]] = {}
        for axis, term_ids in CYANORAK_ROLES.items():
            res = genes_by_ontology(
                ontology="cyanorak_role",
                organism="MED4",
                term_ids=term_ids,
                limit=500,
                min_gene_set_size=1,
            )
            df = to_dataframe(res)
            genes = sorted(df["locus_tag"].unique().tolist())
            cyanorak_sets[axis] = genes
            log(f"  cyanorak[{axis:<12}] n={len(genes)} (terms={term_ids})", fh=fh)

        # Save panel definitions for transparency
        rows = []
        for axis, genes in positive_sets.items():
            for g in genes:
                rows.append({
                    "axis": axis, "panel_kind": "positive",
                    "locus_tag": g, "direction": POSITIVE_DIRECTION[g],
                })
        for axis, genes in cyanorak_sets.items():
            d = CYANORAK_AXIS_DIRECTION[axis]
            for g in genes:
                rows.append({
                    "axis": axis, "panel_kind": "cyanorak",
                    "locus_tag": g, "direction": d,
                })
        panel_def = pd.DataFrame(rows)
        panel_def_path = OUT_DIR / "panel_definitions.csv"
        panel_def.to_csv(panel_def_path, index=False)
        log(f"\nWrote panel definitions to {panel_def_path}", fh=fh)

        # Load experiment list — drop the cross-condition contrast experiment
        experiments = pd.read_csv(EXP_CSV)
        score_experiments = experiments[
            ~experiments["experiment_id"].str.contains("coculture_alteromonas_hot1a3_med4")
        ]["experiment_id"].tolist()
        log(f"\nScoring {len(score_experiments)} experiments:", fh=fh)
        for eid in score_experiments:
            log(f"  - {eid}", fh=fh)

        # Pull genome-wide DE per experiment, save, score
        all_score_rows = []
        for eid in score_experiments:
            log(f"\n[{eid}]", fh=fh)
            result = differential_expression_by_gene(
                experiment_ids=[eid],
                limit=None,
                verbose=True,
            )
            de = to_dataframe(result)
            log(f"  genome-wide rows: {len(de)} (truncated={result['truncated']})", fh=fh)

            # Drop pooled days 60+89 per step 2 decision
            before = len(de)
            de = de[de["timepoint"].astype(str) != "days 60+89"].copy()
            if len(de) < before:
                log(f"  dropped {before - len(de)} pooled 'days 60+89' rows", fh=fh)

            short = short_exp_id(eid)
            out_path = OUT_DIR / f"genome_de_{short}.csv"
            de.to_csv(out_path, index=False)
            log(f"  wrote {out_path}", fh=fh)

            # Determine condition + omics from metadata
            condition = "coculture" if "coculture" in short else "axenic"
            omics = "RNASEQ" if "rnaseq" in short else "PROTEOMICS"

            # For each axis × panel × TP, score
            for axis in positive_sets.keys():
                # Determine TPs in this experiment
                tps = de["timepoint"].dropna().unique().tolist()
                if not tps:
                    # Single-point experiment (axenic-RNA)
                    tps = [None]

                for panel_kind in ("positive", "cyanorak"):
                    if panel_kind == "positive":
                        gene_set = positive_sets[axis]
                        direction = {g: POSITIVE_DIRECTION[g] for g in gene_set}
                    else:
                        if axis not in cyanorak_sets:
                            continue  # no cyanorak panel for this axis (cell_death)
                        gene_set = cyanorak_sets[axis]
                        d = CYANORAK_AXIS_DIRECTION[axis]
                        direction = {g: d for g in gene_set}

                    for tp in tps:
                        if tp is None:
                            tp_de = de
                            tp_label = "single"
                            tp_phase = de["growth_phase"].dropna().iloc[0] if not de["growth_phase"].dropna().empty else None
                            tp_hours = None
                        else:
                            tp_de = de[de["timepoint"] == tp]
                            tp_label = tp
                            tp_phase = tp_de["growth_phase"].dropna().iloc[0] if "growth_phase" in tp_de.columns and not tp_de["growth_phase"].dropna().empty else None
                            tp_hours = tp_de["timepoint_hours"].dropna().iloc[0] if "timepoint_hours" in tp_de.columns and not tp_de["timepoint_hours"].dropna().empty else None

                        res = axis_stress_score(
                            tp_de[["locus_tag", "log2fc"]],
                            axis_genes=gene_set,
                            direction=direction,
                        )
                        all_score_rows.append({
                            "axis": axis,
                            "panel_kind": panel_kind,
                            "experiment_id": eid,
                            "experiment_short": short,
                            "condition": condition,
                            "omics": omics,
                            "timepoint": tp_label,
                            "timepoint_hours": tp_hours,
                            "growth_phase": tp_phase,
                            "axis_score": res["axis_score"],
                            "axis_mean_signed_lfc": res["axis_mean"],
                            "background_mean_lfc": res["background_mean"],
                            "background_sd_lfc": res["background_sd"],
                            "n_axis_with_data": res["n_axis"],
                            "n_axis_total": len(gene_set),
                            "n_background": res["n_background"],
                        })

        scores = pd.DataFrame(all_score_rows)
        scores_path = OUT_DIR / "all_axes_scores.csv"
        scores.to_csv(scores_path, index=False)
        log(f"\nWrote {len(scores)} score rows to {scores_path}", fh=fh)
        log(f"  axes:        {sorted(scores['axis'].unique())}", fh=fh)
        log(f"  panel_kinds: {sorted(scores['panel_kind'].unique())}", fh=fh)
        log(f"  conditions:  {sorted(scores['condition'].unique())}", fh=fh)
        log(f"  omics:       {sorted(scores['omics'].unique())}", fh=fh)

        # Coverage summary: how many axis genes had data per (axis × panel × experiment)
        log("\nCoverage (n_axis_with_data / n_axis_total) per (axis × panel × experiment_short × first-TP):", fh=fh)
        first_tp_per = scores.groupby(["axis", "panel_kind", "experiment_short"]).first().reset_index()
        for _, r in first_tp_per.iterrows():
            log(f"  {r['axis']:<12} {r['panel_kind']:<10} {r['experiment_short']:<25} "
                f"{int(r['n_axis_with_data']):>2}/{int(r['n_axis_total']):>2}", fh=fh)


if __name__ == "__main__":
    main()
