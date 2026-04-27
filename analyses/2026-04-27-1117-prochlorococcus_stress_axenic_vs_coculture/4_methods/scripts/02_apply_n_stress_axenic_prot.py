"""Driving example for step 4: N-stress axis × axenic-proteomics trajectory.

Pulls genome-wide DE for the axenic-proteomics experiment, applies
axis_stress_score per TP using the 5 validated N-stress positives, saves
the per-TP scores plus a trajectory figure.

Inputs (read):
  - 3_analysis_framing/data/control_panel.csv (for the n_stress positives)

Outputs (write):
  - data/n_stress_axenic_prot_scores.csv          one row per TP, score + diagnostics
  - data/n_stress_axenic_prot_de.csv              full DE used (for transparency)
  - figures/n_stress_axenic_prot_trajectory.png   bar plot of axis_score per TP
  - data/02_apply_n_stress_axenic_prot.log

Run from the multiomics_research repo root:
  uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/4_methods/scripts/02_apply_n_stress_axenic_prot.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd

from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe

# Local import: stress_score module from this analysis step's folder
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from stress_score import axis_stress_score  # noqa: E402

ANALYSIS_DIR = Path(__file__).resolve().parents[2]
CONTROL_CSV = ANALYSIS_DIR / "3_analysis_framing" / "data" / "control_panel.csv"

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = OUT_DIR / "02_apply_n_stress_axenic_prot.log"

EXPERIMENT_ID = "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic"


def log(msg: str, *, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def main() -> None:
    with LOG_PATH.open("w") as fh:
        # Load control panel and pick the n_stress positives
        controls = pd.read_csv(CONTROL_CSV)
        n_stress = controls[(controls["axis"] == "n_stress") & (controls["role"] == "positive")]
        axis_genes = n_stress["locus_tag"].tolist()
        # All n_stress positives are validated UP — direction = +1
        direction = {g: 1 for g in axis_genes}

        log(f"Driving example: N-stress × axenic-proteomics", fh=fh)
        log(f"Experiment: {EXPERIMENT_ID}", fh=fh)
        log(f"Axis gene set: {len(axis_genes)} genes (validated positives, all direction = +1)", fh=fh)
        for _, r in n_stress.iterrows():
            log(f"  {r['locus_tag']:>10}  {r['gene_name']}", fh=fh)

        # Pull all genes DE for the experiment
        log("\nFetching genome-wide DE for the experiment ...", fh=fh)
        result = differential_expression_by_gene(
            experiment_ids=[EXPERIMENT_ID],
            limit=None,
            verbose=True,
        )
        de = to_dataframe(result)
        log(f"  total_matching = {result['total_matching']}", fh=fh)
        log(f"  rows returned  = {len(de)}", fh=fh)
        log(f"  truncated      = {result['truncated']}", fh=fh)

        # Drop pooled "days 60+89" rows (not present in axenic-proteomics anyway, but defensive)
        before = len(de)
        de = de[de["timepoint"].astype(str) != "days 60+89"].copy()
        if len(de) < before:
            log(f"  dropped {before - len(de)} pooled 'days 60+89' rows", fh=fh)

        # Save full DE for transparency
        de_path = OUT_DIR / "n_stress_axenic_prot_de.csv"
        de.to_csv(de_path, index=False)
        log(f"\nWrote full DE to {de_path}", fh=fh)
        log(f"  unique genes:      {de['locus_tag'].nunique()}", fh=fh)
        log(f"  unique timepoints: {sorted(de['timepoint'].dropna().unique())}", fh=fh)

        # Score per TP
        rows = []
        for tp, group in de.groupby("timepoint", sort=False):
            res = axis_stress_score(
                group[["locus_tag", "log2fc"]],
                axis_genes=axis_genes,
                direction=direction,
            )
            growth_phase = group["growth_phase"].dropna().iloc[0] if "growth_phase" in group.columns and not group["growth_phase"].dropna().empty else None
            tp_hours = group["timepoint_hours"].dropna().iloc[0] if "timepoint_hours" in group.columns and not group["timepoint_hours"].dropna().empty else None
            rows.append({
                "timepoint": tp,
                "timepoint_hours": tp_hours,
                "growth_phase": growth_phase,
                "axis_score": res["axis_score"],
                "axis_mean_signed_lfc": res["axis_mean"],
                "background_mean_lfc": res["background_mean"],
                "background_sd_lfc": res["background_sd"],
                "n_axis_genes_with_data": res["n_axis"],
                "n_background_genes_with_data": res["n_background"],
                "axis_genes_with_data": ";".join(res["axis_genes_with_data"]),
                "axis_genes_missing_data": ";".join(res["axis_genes_missing_data"]),
            })
        scores = pd.DataFrame(rows).sort_values("timepoint_hours").reset_index(drop=True)
        scores_path = OUT_DIR / "n_stress_axenic_prot_scores.csv"
        scores.to_csv(scores_path, index=False)
        log(f"\nWrote per-TP scores to {scores_path}", fh=fh)
        log(scores.to_string(index=False), fh=fh)

        # Worked example for the smallest-magnitude TP (proves the formula end-to-end)
        log("\nWorked example for first TP (day 14, nutrient_limited):", fh=fh)
        first = scores.iloc[0]
        first_tp_de = de[de["timepoint"] == first["timepoint"]]
        log("  Per-axis-gene log2fc:", fh=fh)
        for g in axis_genes:
            r = first_tp_de[first_tp_de["locus_tag"] == g]
            if r.empty:
                log(f"    {g}: missing", fh=fh)
            else:
                log(f"    {g}: log2fc={float(r.iloc[0]['log2fc']):+.4f}", fh=fh)
        log(f"  axis_mean_signed_lfc      = {first['axis_mean_signed_lfc']:.4f}", fh=fh)
        log(f"  background_mean_lfc       = {first['background_mean_lfc']:.4f}", fh=fh)
        log(f"  background_sd_lfc         = {first['background_sd_lfc']:.4f}", fh=fh)
        log(f"  axis_score                = ({first['axis_mean_signed_lfc']:.4f} - {first['background_mean_lfc']:.4f}) / {first['background_sd_lfc']:.4f} = {first['axis_score']:.4f}", fh=fh)

        # Two-panel trajectory: raw axis mean (left) and signed-Z (right).
        # Showing both surfaces the methodology insight that the axis can stay
        # raw-strong while losing distinctiveness as the global response widens.
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
        x = (scores["timepoint_hours"].astype(float) / 24.0).to_numpy()
        colors = ["#5b9bd5" if p == "nutrient_limited" else "#d04545"
                  for p in scores["growth_phase"]]
        bar_width = 5  # days

        # Left panel: raw axis mean signed log2FC
        ax = axes[0]
        y_raw = scores["axis_mean_signed_lfc"].to_numpy()
        ax.bar(x, y_raw, width=bar_width, color=colors, edgecolor="black", linewidth=0.6)
        for xi, yi, lbl, ph in zip(x, y_raw, scores["timepoint"], scores["growth_phase"]):
            ax.text(xi, yi + 0.05, f"{lbl}\n({ph[:4]})",
                    ha="center", va="bottom", fontsize=8)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xlabel("timepoint (days)")
        ax.set_ylabel("axis_mean signed log2FC\n(mean of d_g·log2fc over axis genes)")
        ax.set_title("Raw axis response (not normalized)")
        ax.set_ylim(0, max(y_raw.max() * 1.25, 1.6))

        # Right panel: signed-Z axis score
        ax = axes[1]
        y_z = scores["axis_score"].to_numpy()
        ax.bar(x, y_z, width=bar_width, color=colors, edgecolor="black", linewidth=0.6)
        for xi, yi, lbl, ph in zip(x, y_z, scores["timepoint"], scores["growth_phase"]):
            ax.text(xi, yi + 0.05, f"{lbl}\n({ph[:4]})",
                    ha="center", va="bottom", fontsize=8)
        ax.axhline(0, color="black", lw=0.8)
        ax.axhline(2, color="gray", lw=0.5, ls="--")
        ax.text(x[-1], 2.05, "z=+2", color="gray", fontsize=8, ha="right")
        ax.set_xlabel("timepoint (days)")
        ax.set_ylabel("axis_score (signed-Z vs genome background)")
        ax.set_title("Distinctiveness vs genome (z-score)")
        ax.set_ylim(0, max(y_z.max() * 1.25, 2.4))

        fig.suptitle(
            "Driving example: N-stress axis × axenic proteomics\n"
            "Pro MED4, Weissberg 2025 — within-condition trajectory (vs each condition's own exponential)",
            fontsize=11,
            y=1.02,
        )
        plt.tight_layout()
        out_png = FIG_DIR / "n_stress_axenic_prot_trajectory.png"
        out_pdf = FIG_DIR / "n_stress_axenic_prot_trajectory.pdf"
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.savefig(out_pdf, bbox_inches="tight")
        log(f"\nWrote {out_png}", fh=fh)
        log(f"Wrote {out_pdf}", fh=fh)


if __name__ == "__main__":
    main()
