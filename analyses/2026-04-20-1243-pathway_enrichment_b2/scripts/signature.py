"""Shared signature-derivation primitive.

Pure functions — no I/O, no logging. Imported by 04_derive_signature.py (main
run) and 05_compute_scores.py (LOO-R stability re-derivation). Kept as a
single-file module, not a package, to match spec §9's "no _utils/ package"
discipline while still avoiding cross-script importlib hacks.
"""
from __future__ import annotations

import json
from collections import defaultdict

import pandas as pd

# Temporal filter: drop R clusters with timepoint_hours <= cutoff.
# Keep-mask at call sites: `hours.isna() | (hours > TIMEPOINT_HOURS_CUTOFF)`
# per decisions.md D1 (strict greater-than: excludes 0h AND 3h).
TIMEPOINT_HOURS_CUTOFF = 3.0
CORE_SUPPORT = 3
FALLBACK_SUPPORT = 2
NOTABLE_SIGNED_SCORE = 3.0


def derive_for_ontology(
    df: pd.DataFrame, support_threshold: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (signature_df, dropped_df) for one ontology.

    Caller is responsible for upstream filtering (restricting to R clusters of
    one ontology, applying any temporal filter). Required columns on df:
    term_id, term_name, cluster, direction, experiment_id, p_adjust,
    signed_score.
    """
    sig_rows = []
    dropped_rows = []
    for (term_id, term_name), grp in df.groupby(["term_id", "term_name"]):
        sig = grp[grp["p_adjust"] < 0.05]
        n_up = (sig["direction"] == "up").sum()
        n_down = (sig["direction"] == "down").sum()
        up_clusters = sig.loc[sig["direction"] == "up", "cluster"].tolist()
        down_clusters = sig.loc[sig["direction"] == "down", "cluster"].tolist()
        max_abs = grp["signed_score"].abs().max() if len(grp) else 0.0

        per_exp = defaultdict(lambda: {"up": 0, "down": 0})
        for _, r in sig.iterrows():
            per_exp[r["experiment_id"]][r["direction"]] += 1
        per_exp_json = json.dumps(dict(per_exp))

        is_up = n_up >= support_threshold and n_down < support_threshold
        is_down = n_down >= support_threshold and n_up < support_threshold
        is_bidirectional = n_up >= support_threshold and n_down >= support_threshold

        if is_up:
            sig_rows.append({
                "term_id": term_id, "term_name": term_name, "direction": "up",
                "n_clusters_supporting": n_up, "n_up": n_up, "n_down": n_down,
                "contributing_clusters": "|".join(up_clusters),
                "per_experiment_breakdown": per_exp_json,
            })
        elif is_down:
            sig_rows.append({
                "term_id": term_id, "term_name": term_name, "direction": "down",
                "n_clusters_supporting": n_down, "n_up": n_up, "n_down": n_down,
                "contributing_clusters": "|".join(down_clusters),
                "per_experiment_breakdown": per_exp_json,
            })
        elif is_bidirectional:
            dropped_rows.append({
                "term_id": term_id, "term_name": term_name,
                "n_up": n_up, "n_down": n_down,
                "drop_reason": "bidirectional",
                "max_signed_score": max_abs,
                "contributing_clusters_up": "|".join(up_clusters),
                "contributing_clusters_down": "|".join(down_clusters),
                "per_experiment_breakdown": per_exp_json,
            })
        elif max(n_up, n_down) < support_threshold and (
            max_abs >= NOTABLE_SIGNED_SCORE or (n_up + n_down) >= 2
        ):
            dropped_rows.append({
                "term_id": term_id, "term_name": term_name,
                "n_up": n_up, "n_down": n_down,
                "drop_reason": "below_threshold_notable",
                "max_signed_score": max_abs,
                "contributing_clusters_up": "|".join(up_clusters),
                "contributing_clusters_down": "|".join(down_clusters),
                "per_experiment_breakdown": per_exp_json,
            })

    return pd.DataFrame(sig_rows), pd.DataFrame(dropped_rows)
