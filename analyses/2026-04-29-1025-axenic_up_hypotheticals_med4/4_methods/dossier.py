"""
Dossier methodology module — builds and renders one card per candidate gene.

Implements the 7-axis dossier surface locked at step 3:
  1. Identity & DE evidence (this experiment)
  2. Cluster memberships (with curated functional/dynamics/temporal text)
  3. Ortholog groups (per group at all 4 specificity ranks)
  4. Ortholog group signal probe (cached per group_id; envelope-only)
  5. Cross-study response profile (2-call: locked-only + full; adjust+flag)
  6. Ontology terms (leaves only; rollup deferred)
  7. Derived metrics

Public surface:
    build_dossier(locus_tag, locked_experiment_id, organism, group_cache) -> dict
    render_card_markdown(card) -> str
    GroupProbeCache(cache_path=None)

Each call to build_dossier issues:
  - 1 gene_overview (verbose)
  - 1 differential_expression_by_gene (this experiment only)
  - 1 gene_clusters_by_gene (verbose)
  - 1 gene_homologs (verbose)
  - 0–4 group probes (cached by group_id; each probe = 1 genes_by_homolog_group + 1 gene_overview summary)
  - 1 gene_ontology_terms (leaves only)
  - 1 gene_derived_metrics (verbose)
  - 2 gene_response_profile (locked-only + full)
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from multiomics_explorer import (
    differential_expression_by_gene,
    gene_clusters_by_gene,
    gene_derived_metrics,
    gene_homologs,
    gene_ontology_terms,
    gene_overview,
    gene_response_profile,
    genes_by_homolog_group,
)

# ---------- Group signal probe cache ----------


@dataclass
class GroupProbeCache:
    """Cache `genes_by_homolog_group` + `gene_overview(summary=True)` envelopes by group_id.

    Stored value per group_id:
      - dict (envelope summary) when the group has members beyond the candidate.
      - None when the group is a singleton (only the candidate itself).
    """

    cache_path: Path | None = None
    _cache: dict[str, dict | None] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.cache_path and self.cache_path.exists():
            try:
                self._cache = json.loads(self.cache_path.read_text())
            except json.JSONDecodeError:
                self._cache = {}

    def get_or_compute(self, group_id: str, candidate_locus_tag: str) -> dict | None:
        if group_id in self._cache:
            return self._cache[group_id]

        members_resp = genes_by_homolog_group(group_ids=[group_id], limit=None)
        member_lts = sorted({
            r["locus_tag"]
            for r in members_resp.get("results", [])
            if r.get("locus_tag") and r["locus_tag"] != candidate_locus_tag
        })

        if not member_lts:
            self._cache[group_id] = None
            return None

        env = gene_overview(locus_tags=member_lts, summary=True)

        probe = {
            "n_members_excl_candidate": len(member_lts),
            "by_organism": {d["organism_name"]: d["count"] for d in env.get("by_organism", [])},
            "by_category": {d["category"]: d["count"] for d in env.get("by_category", [])},
            "by_annotation_type": {
                d["annotation_type"]: d["count"] for d in env.get("by_annotation_type", [])
            },
            "has_expression": int(env.get("has_expression", 0) or 0),
            "has_significant_expression": int(env.get("has_significant_expression", 0) or 0),
            "has_clusters": int(env.get("has_clusters", 0) or 0),
            "has_derived_metrics": int(env.get("has_derived_metrics", 0) or 0),
            "has_orthologs": int(env.get("has_orthologs", 0) or 0),
        }
        self._cache[group_id] = probe
        return probe

    def save(self) -> None:
        if self.cache_path is not None:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            self.cache_path.write_text(json.dumps(self._cache, indent=2, sort_keys=True))


# ---------- Section builders ----------


def _identity(ov_row: dict) -> dict:
    return {
        "locus_tag": ov_row.get("locus_tag"),
        "gene_name": ov_row.get("gene_name"),
        "product": ov_row.get("product"),
        "gene_category": ov_row.get("gene_category"),
        "annotation_quality": ov_row.get("annotation_quality"),
        "organism_name": ov_row.get("organism_name"),
        "all_identifiers": ov_row.get("all_identifiers") or [],
        "annotation_types": ov_row.get("annotation_types") or [],
        "function_description": ov_row.get("function_description"),
    }


def _de_evidence(de_resp: dict, locked_experiment_id: str) -> dict:
    rows = de_resp.get("results") or []
    if not rows:
        return {"experiment_id": locked_experiment_id, "no_data": True}
    r = rows[0]
    return {
        "experiment_id": r.get("experiment_id"),
        "timepoint": r.get("timepoint"),
        "timepoint_hours": r.get("timepoint_hours"),
        "log2fc": r.get("log2fc"),
        "padj": r.get("padj"),
        "rank": r.get("rank"),
        "rank_up": r.get("rank_up"),
        "rank_down": r.get("rank_down"),
        "expression_status": r.get("expression_status"),
    }


def _clusters(clu_resp: dict) -> list[dict]:
    out = []
    for r in clu_resp.get("results") or []:
        out.append({
            "cluster_id": r.get("cluster_id"),
            "cluster_name": r.get("cluster_name"),
            "cluster_type": r.get("cluster_type"),
            "cluster_method": r.get("cluster_method"),
            "member_count": r.get("member_count"),
            "cluster_functional_description": r.get("cluster_functional_description"),
            "cluster_expression_dynamics": r.get("cluster_expression_dynamics"),
            "cluster_temporal_pattern": r.get("cluster_temporal_pattern"),
            "treatment": r.get("treatment"),
            "analysis_id": r.get("analysis_id"),
            "analysis_name": r.get("analysis_name"),
        })
    return out


def _ortholog_groups(hom_resp: dict, locus_tag: str, cache: GroupProbeCache) -> list[dict]:
    out = []
    for r in hom_resp.get("results") or []:
        group_id = r.get("group_id")
        probe = cache.get_or_compute(group_id, locus_tag) if group_id else None
        out.append({
            "group_id": group_id,
            "source": r.get("source"),
            "taxonomic_level": r.get("taxonomic_level"),
            "specificity_rank": r.get("specificity_rank"),
            "consensus_gene_name": r.get("consensus_gene_name"),
            "consensus_product": r.get("consensus_product"),
            "member_count": r.get("member_count"),
            "organism_count": r.get("organism_count"),
            "genera": r.get("genera") or [],
            "has_cross_genus_members": r.get("has_cross_genus_members"),
            "description": r.get("description"),
            "functional_description": r.get("functional_description"),
            "consensus_cyanorak_roles": r.get("cyanorak_roles") or [],
            "consensus_cog_categories": r.get("cog_categories") or [],
            "signal_probe": probe,  # None if singleton; dict otherwise
        })
    out.sort(key=lambda g: (g.get("specificity_rank") is None, g.get("specificity_rank") if g.get("specificity_rank") is not None else 99))
    return out


def _treatment_summary(tname: str, ts: dict) -> dict:
    """Flatten one entry from response_summary."""
    return {
        "treatment": tname,
        "experiments_total": ts.get("experiments_total", 0),
        "experiments_tested": ts.get("experiments_tested", 0),
        "experiments_up": ts.get("experiments_up", 0),
        "experiments_down": ts.get("experiments_down", 0),
        "timepoints_total": ts.get("timepoints_total", 0),
        "timepoints_tested": ts.get("timepoints_tested", 0),
        "timepoints_up": ts.get("timepoints_up", 0),
        "timepoints_down": ts.get("timepoints_down", 0),
        "up_best_rank": ts.get("up_best_rank"),
        "up_median_rank": ts.get("up_median_rank"),
        "up_max_log2fc": ts.get("up_max_log2fc"),
        "down_best_rank": ts.get("down_best_rank"),
        "down_median_rank": ts.get("down_median_rank"),
        "down_max_log2fc": ts.get("down_max_log2fc"),
    }


def _adjust_nitrogen_row(full_summary: dict, locked_summary: dict) -> dict | None:
    """Subtract the locked experiment's contribution from full nitrogen counts.

    Cleanly subtractable: experiments_*, timepoints_*.
    Not subtractable: up_best_rank, up_max_log2fc, down_*. Flagged.
    """
    if not full_summary:
        return None
    if not locked_summary:
        # No locked nitrogen contribution — no adjustment needed.
        return {
            "other_experiments": dict(full_summary),
            "extremes_may_include_locked": False,
            "notes": "no locked experiment contribution to subtract",
        }

    other = {
        "experiments_total": full_summary["experiments_total"] - locked_summary["experiments_total"],
        "experiments_tested": full_summary["experiments_tested"] - locked_summary["experiments_tested"],
        "experiments_up": full_summary["experiments_up"] - locked_summary["experiments_up"],
        "experiments_down": full_summary["experiments_down"] - locked_summary["experiments_down"],
        "timepoints_total": full_summary["timepoints_total"] - locked_summary["timepoints_total"],
        "timepoints_tested": full_summary["timepoints_tested"] - locked_summary["timepoints_tested"],
        "timepoints_up": full_summary["timepoints_up"] - locked_summary["timepoints_up"],
        "timepoints_down": full_summary["timepoints_down"] - locked_summary["timepoints_down"],
    }
    # Sanity check
    for k, v in other.items():
        assert v >= 0, f"adjustment produced negative {k}={v}"

    return {
        "other_experiments": other,
        "extremes_may_include_locked": True,
        "extremes_kept_from_full": {
            "up_best_rank": full_summary.get("up_best_rank"),
            "up_max_log2fc": full_summary.get("up_max_log2fc"),
            "down_best_rank": full_summary.get("down_best_rank"),
            "down_max_log2fc": full_summary.get("down_max_log2fc"),
        },
        "locked_contribution": locked_summary,
    }


def _response_profile(rp_locked: dict, rp_full: dict) -> dict:
    """Build the response section: locked-only profile + full cross-study + adjust+flag for nitrogen row."""
    locked_results = rp_locked.get("results") or []
    full_results = rp_full.get("results") or []

    locked_per_treatment = {}
    if locked_results:
        rs = locked_results[0].get("response_summary") or {}
        locked_per_treatment = {t: _treatment_summary(t, ts) for t, ts in rs.items()}

    full_per_treatment = {}
    full_groups = {
        "groups_responded": [],
        "groups_not_responded": [],
        "groups_tested_not_responded": [],
        "groups_not_known": [],
    }
    if full_results:
        r = full_results[0]
        rs = r.get("response_summary") or {}
        full_per_treatment = {t: _treatment_summary(t, ts) for t, ts in rs.items()}
        for k in full_groups:
            full_groups[k] = r.get(k) or []

    nitrogen_adjusted = _adjust_nitrogen_row(
        full_per_treatment.get("nitrogen"), locked_per_treatment.get("nitrogen")
    )

    return {
        "locked": {
            "per_treatment": locked_per_treatment,
        },
        "full": {
            "per_treatment": full_per_treatment,
            **full_groups,
        },
        "nitrogen_adjusted": nitrogen_adjusted,
    }


def _ontology_terms(ont_resp: dict) -> list[dict]:
    out = []
    for r in ont_resp.get("results") or []:
        out.append({
            "term_id": r.get("term_id"),
            "term_name": r.get("term_name"),
            "ontology_type": r.get("ontology_type"),
            "level": r.get("level"),
            "tree": r.get("tree"),
            "tree_code": r.get("tree_code"),
        })
    out.sort(key=lambda t: (t.get("ontology_type") or "", t.get("level") or 0))
    return out


def _derived_metrics(dm_resp: dict) -> list[dict]:
    out = []
    for r in dm_resp.get("results") or []:
        out.append({
            "derived_metric_id": r.get("derived_metric_id"),
            "name": r.get("name"),
            "value": r.get("value"),
            "value_kind": r.get("value_kind"),
            "rank_by_metric": r.get("rank_by_metric"),
            "metric_percentile": r.get("metric_percentile"),
            "metric_bucket": r.get("metric_bucket"),
            "metric_type": r.get("metric_type"),
            "compartment": r.get("compartment"),
            "treatment_type": r.get("treatment_type") or [],
            "publication_doi": r.get("publication_doi"),
            "field_description": r.get("field_description"),
            "unit": r.get("unit"),
        })
    return out


# ---------- Top-level builder ----------


def build_dossier(
    locus_tag: str,
    locked_experiment_id: str,
    organism: str,
    group_cache: GroupProbeCache,
) -> dict:
    """Build a dossier card for one candidate gene.

    Returns a nested dict matching the 7-axis structure. Use render_card_markdown
    to emit a human-readable card.
    """
    ov = gene_overview(locus_tags=[locus_tag], verbose=True)
    if not ov.get("results"):
        return {
            "locus_tag": locus_tag,
            "error": f"{locus_tag} not found in KG (gene_overview returned 0 rows)",
        }
    ov_row = ov["results"][0]

    de = differential_expression_by_gene(
        locus_tags=[locus_tag], experiment_ids=[locked_experiment_id], verbose=False
    )
    clu = gene_clusters_by_gene(locus_tags=[locus_tag], verbose=True, limit=None)
    hom = gene_homologs(locus_tags=[locus_tag], verbose=True, limit=None)
    ont = gene_ontology_terms(locus_tags=[locus_tag], organism=organism, limit=None)
    dm = gene_derived_metrics(locus_tags=[locus_tag], verbose=True, limit=None)
    rp_locked = gene_response_profile(
        locus_tags=[locus_tag], experiment_ids=[locked_experiment_id], group_by="treatment_type"
    )
    rp_full = gene_response_profile(locus_tags=[locus_tag], group_by="treatment_type")

    return {
        "locus_tag": locus_tag,
        "identity": _identity(ov_row),
        "de_evidence": _de_evidence(de, locked_experiment_id),
        "clusters": _clusters(clu),
        "ortholog_groups": _ortholog_groups(hom, locus_tag, group_cache),
        "response": _response_profile(rp_locked, rp_full),
        "ontology": _ontology_terms(ont),
        "derived_metrics": _derived_metrics(dm),
    }


# ---------- Markdown renderer ----------


def _fmt(v: Any, default: str = "—") -> str:
    if v is None or v == "" or v == [] or v == {}:
        return default
    return str(v)


def _fmt_na(v: Any) -> str:
    """Render N/A / None / empty as em-dash; otherwise the value."""
    if v is None or v == "" or (isinstance(v, str) and v.strip().upper() == "N/A"):
        return "—"
    return str(v)


def render_card_markdown(card: dict) -> str:
    if "error" in card:
        return f"# Dossier — {card['locus_tag']}\n\n**Error:** {card['error']}\n"

    lines: list[str] = []
    ident = card["identity"]
    title = ident["locus_tag"]
    if ident.get("gene_name"):
        title = f"{title} ({ident['gene_name']})"
    lines.append(f"# Dossier — {title}\n")

    # Identity
    lines.append("## Identity\n")
    lines.append(f"- locus_tag: `{ident['locus_tag']}`")
    lines.append(f"- gene_name: {_fmt(ident.get('gene_name'))}")
    lines.append(f"- product: {_fmt(ident.get('product'))}")
    lines.append(f"- gene_category: {_fmt(ident.get('gene_category'))}")
    lines.append(f"- annotation_quality: {ident.get('annotation_quality')}")
    lines.append(f"- organism_name: {ident.get('organism_name')}")
    lines.append(f"- all_identifiers: {', '.join(ident.get('all_identifiers') or []) or '—'}")
    lines.append(f"- annotation_types: {', '.join(ident.get('annotation_types') or []) or '—'}")
    if ident.get("function_description"):
        lines.append(f"- function_description: {ident['function_description']}")
    lines.append("")

    # DE evidence
    de = card["de_evidence"]
    lines.append("## DE evidence — this experiment\n")
    if de.get("no_data"):
        lines.append(f"No DE row returned for `{de['experiment_id']}` — gene is not in the experiment's DE table.")
    else:
        lines.append(f"- experiment_id: `{de['experiment_id']}`")
        lines.append(f"- timepoint: {_fmt(de.get('timepoint'))} ({_fmt(de.get('timepoint_hours'))} h)")
        lines.append(f"- log2fc: {de.get('log2fc'):.3f}" if de.get("log2fc") is not None else "- log2fc: —")
        lines.append(f"- padj: {de.get('padj'):.3e}" if de.get("padj") is not None else "- padj: —")
        lines.append(f"- rank (within experiment, by |log2fc|): {de.get('rank')}")
        lines.append(f"- rank_up: {de.get('rank_up')}")
        lines.append(f"- expression_status: **{de.get('expression_status')}**")
    lines.append("")

    # Clusters
    clusters = card["clusters"]
    lines.append(f"## Cluster memberships ({len(clusters)})\n")
    if not clusters:
        lines.append("No cluster memberships in any clustering analysis.\n")
    else:
        for c in clusters:
            lines.append(f"### `{c['cluster_id']}`")
            lines.append(f"- name: {c['cluster_name']}")
            lines.append(f"- type: {c['cluster_type']} / method: {c.get('cluster_method') or '—'}")
            lines.append(f"- member_count: {c.get('member_count')}")
            lines.append(f"- cluster_functional_description: {_fmt_na(c.get('cluster_functional_description'))}")
            lines.append(f"- cluster_expression_dynamics: {_fmt_na(c.get('cluster_expression_dynamics'))}")
            lines.append(f"- cluster_temporal_pattern: {_fmt_na(c.get('cluster_temporal_pattern'))}")
            lines.append(f"- treatment: {_fmt(c.get('treatment'))}")
            lines.append(f"- analysis: {c['analysis_name']} (`{c['analysis_id']}`)")
            lines.append("")

    # Ortholog groups + signal probe
    ogs = card["ortholog_groups"]
    lines.append(f"## Ortholog groups ({len(ogs)})\n")
    if not ogs:
        lines.append("No homolog groups in any source. (See `gaps_and_friction.md` F3 if RefSeq-only locus tag.)\n")
    else:
        for og in ogs:
            lines.append(f"### `{og['group_id']}` — {og['source']}, {og['taxonomic_level']} (rank {og['specificity_rank']})")
            lines.append(f"- consensus_gene_name: {_fmt(og.get('consensus_gene_name'))}")
            lines.append(f"- consensus_product: {_fmt(og.get('consensus_product'))}")
            lines.append(f"- member_count: {og['member_count']} / organism_count: {og['organism_count']} / has_cross_genus_members: {og.get('has_cross_genus_members')}")
            lines.append(f"- genera: {', '.join(og.get('genera') or []) or '—'}")
            if og.get("description"):
                lines.append(f"- description: {og['description']}")
            if og.get("functional_description"):
                lines.append(f"- functional_description: {og['functional_description']}")
            cy = og.get("consensus_cyanorak_roles") or []
            if cy:
                lines.append(f"- consensus cyanorak_roles: {'; '.join(r.get('name', '?') for r in cy)}")
            cg = og.get("consensus_cog_categories") or []
            if cg:
                lines.append(f"- consensus cog_categories: {'; '.join(r.get('name', '?') for r in cg)}")
            probe = og.get("signal_probe")
            if probe is None:
                lines.append("- **signal probe:** singleton group (only the candidate gene itself) — no orthologs to probe")
            else:
                lines.append(f"- **signal probe** (envelope from `gene_overview(locus_tags=[{probe['n_members_excl_candidate']} ortholog members], summary=True)`):")
                lines.append(f"    - by_organism: {', '.join(f'{org}={n}' for org, n in probe['by_organism'].items())}")
                lines.append(f"    - by_category: {', '.join(f'{c}={n}' for c, n in probe['by_category'].items()) or '—'}")
                lines.append(f"    - by_annotation_type: {', '.join(f'{a}={n}' for a, n in probe['by_annotation_type'].items()) or '—'}")
                lines.append(f"    - has_expression: {probe['has_expression']} / {probe['n_members_excl_candidate']}")
                lines.append(f"    - has_significant_expression: {probe['has_significant_expression']} / {probe['n_members_excl_candidate']}")
                lines.append(f"    - has_clusters: {probe['has_clusters']} / {probe['n_members_excl_candidate']}")
                lines.append(f"    - has_derived_metrics: {probe['has_derived_metrics']} / {probe['n_members_excl_candidate']}")
            lines.append("")

    # Response profile
    resp = card["response"]
    lines.append("## Cross-study response profile\n")

    locked_pt = resp["locked"].get("per_treatment", {})
    if not locked_pt:
        lines.append("- This experiment: no DE row found in `gene_response_profile(experiment_ids=[locked])`.\n")
    else:
        lines.append("### This experiment (locked)\n")
        for tname, ts in locked_pt.items():
            up_max = ts.get("up_max_log2fc")
            up_max_str = f"{up_max:.3f}" if isinstance(up_max, (int, float)) else "—"
            lines.append(
                f"- {tname}: experiments_tested={ts['experiments_tested']}/{ts['experiments_total']}, "
                f"timepoints_up={ts['timepoints_up']}/{ts['timepoints_tested']}, "
                f"up_best_rank={ts.get('up_best_rank') or '—'}, up_max_log2fc={up_max_str}"
            )
        lines.append("")

    full = resp["full"]
    lines.append("### Full cross-study (all MED4 experiments — INCLUDES this experiment)\n")
    lines.append(f"- groups_responded: {full.get('groups_responded') or '—'}")
    lines.append(f"- groups_not_responded: {full.get('groups_not_responded') or '—'}")
    lines.append(f"- groups_tested_not_responded: {full.get('groups_tested_not_responded') or '—'}")
    lines.append(f"- groups_not_known: {full.get('groups_not_known') or '—'}\n")

    full_pt = full.get("per_treatment", {})
    if full_pt:
        lines.append("Per-treatment summary:\n")
        for tname, ts in full_pt.items():
            up_max = ts.get("up_max_log2fc")
            up_max_str = f"{up_max:.3f}" if isinstance(up_max, (int, float)) else "—"
            down_max = ts.get("down_max_log2fc")
            down_max_str = f"{down_max:.3f}" if isinstance(down_max, (int, float)) else "—"
            lines.append(
                f"- **{tname}** — exps_tested={ts['experiments_tested']}/{ts['experiments_total']}, "
                f"up_exps={ts['experiments_up']}, down_exps={ts['experiments_down']}, "
                f"timepoints_up={ts['timepoints_up']}/{ts['timepoints_tested']}, "
                f"timepoints_down={ts['timepoints_down']}, "
                f"up_best_rank={ts.get('up_best_rank') or '—'}, up_max_log2fc={up_max_str}, "
                f"down_max_log2fc={down_max_str}"
            )
        lines.append("")

    # Nitrogen adjusted
    adj = resp.get("nitrogen_adjusted")
    if adj:
        lines.append("### Nitrogen — adjust+flag for the locked experiment overlap\n")
        if adj.get("notes"):
            lines.append(f"_{adj['notes']}_\n")
        lc = adj.get("locked_contribution") or {}
        if lc:
            lines.append("- locked experiment contributed: "
                         f"experiments_total={lc['experiments_total']}, "
                         f"experiments_up={lc['experiments_up']}, "
                         f"experiments_down={lc['experiments_down']}, "
                         f"timepoints_total={lc['timepoints_total']}, "
                         f"timepoints_up={lc['timepoints_up']}, "
                         f"timepoints_down={lc['timepoints_down']}")
        oth = adj.get("other_experiments") or {}
        if oth:
            lines.append("- **other-than-locked** (cross-study minus locked): "
                         f"experiments_total={oth['experiments_total']}, "
                         f"experiments_tested={oth['experiments_tested']}, "
                         f"experiments_up={oth['experiments_up']}, "
                         f"experiments_down={oth['experiments_down']}, "
                         f"timepoints_total={oth['timepoints_total']}, "
                         f"timepoints_tested={oth['timepoints_tested']}, "
                         f"timepoints_up={oth['timepoints_up']}, "
                         f"timepoints_down={oth['timepoints_down']}")
        if adj.get("extremes_may_include_locked"):
            ex = adj.get("extremes_kept_from_full") or {}
            lines.append(f"- ⚠ extremes (up_best_rank={ex.get('up_best_rank')}, "
                         f"up_max_log2fc={ex.get('up_max_log2fc')}) "
                         f"may include or be from this analysis's locked experiment — "
                         f"compare with the locked-only stats above")
        lines.append("")

    # Ontology
    ont = card["ontology"]
    lines.append(f"## Ontology terms ({len(ont)} — leaves only)\n")
    if not ont:
        lines.append("No ontology terms in any source. (Member of the F1 17-gene no-ontology subset.)\n")
    else:
        for t in ont:
            extra = ""
            if t.get("tree"):
                extra = f" / tree={t['tree']}"
            lines.append(f"- **{t['ontology_type']}** (level {t.get('level')}): {t['term_name']} `{t['term_id']}`{extra}")
        lines.append("")

    # Derived metrics
    dms = card["derived_metrics"]
    lines.append(f"## Derived metrics ({len(dms)})\n")
    if not dms:
        lines.append("No derived metrics.\n")
    else:
        for d in dms:
            v = d.get("value")
            v_str = f"{v:.3f}" if isinstance(v, (int, float)) else str(v)
            lines.append(f"### {d['name']}")
            lines.append(f"- derived_metric_id: `{d['derived_metric_id']}`")
            lines.append(f"- value: {v_str} ({d.get('value_kind')}); rank={d.get('rank_by_metric')}, "
                         f"percentile={d.get('metric_percentile')}, bucket={d.get('metric_bucket')}")
            lines.append(f"- metric_type: {d.get('metric_type')} / compartment: {d.get('compartment')}")
            lines.append(f"- treatment_type: {', '.join(d.get('treatment_type') or [])} / publication_doi: {d.get('publication_doi')}")
            if d.get("field_description"):
                lines.append(f"- description: {d['field_description']}")
            lines.append("")

    return "\n".join(lines)
