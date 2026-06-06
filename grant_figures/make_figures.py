#!/usr/bin/env python
"""Grant-proposal figures for the multi-omics KG + MCP + AI-research project.

Generates three vector figures (PDF + SVG):
  fig1_architecture  — the 3-layer system: KG -> MCP server -> AI research agent
  fig2_workflow      — the 6-step research methodology with do/show/explore/decide
  fig3_exemplar      — concrete science win (MED4 axenic vs coculture, N-limitation)

All numbers are grounded in the live KG (kg_release_info / concepts doc) and the
committed analyses under analyses/. Run from repo root:

    uv run grant_figures/make_figures.py
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

# --------------------------------------------------------------------------- #
# Palette — marine / cyanobacteria theme
# --------------------------------------------------------------------------- #
INK = "#1A2238"          # near-black text
MUTED = "#5A6478"        # secondary text
PAGE = "#FFFFFF"

KG = "#0E7C7B"           # knowledge-graph layer (teal)
KG_L = "#D6ECEC"
MCP = "#2A6F97"          # MCP/API layer (blue)
MCP_L = "#D8E7F1"
RES = "#5B4B8A"          # research/agent layer (indigo)
RES_L = "#E5E0F0"
GOLD = "#E0A100"         # highlight / accent
GOLD_L = "#FBEFC9"
CORAL = "#D7553B"        # contrast accent
GREEN = "#2E8B57"        # "alive" / positive
GREY_L = "#EEF1F5"

# DejaVu Sans is bundled with matplotlib, complete, and embeds cleanly in PDF.
# (System "Helvetica" on this box is a partial font that breaks PDF subsetting.)
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["svg.fonttype"] = "none"  # keep text as text in SVG (editable)
plt.rcParams["pdf.fonttype"] = 42       # embed TrueType in PDF (editable text)

OUT = Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- #
# Drawing helpers (axes in 0..100 user coords)
# --------------------------------------------------------------------------- #
def new_canvas(w_in, h_in):
    fig, ax = plt.subplots(figsize=(w_in, h_in))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    return fig, ax


def box(ax, x, y, w, h, *, fc, ec, lw=1.4, r=0.025, z=2, alpha=1.0):
    # r is a fraction; rounding is in data units, capped so it can't exceed a
    # half-dimension (which would balloon the corner into a chevron blob).
    rs = min(r * 100, 0.5 * min(w, h))
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={rs}",
        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z, alpha=alpha,
        mutation_aspect=1,
    )
    ax.add_patch(p)
    return p


def text(ax, x, y, s, *, size=10, color=INK, weight="normal", ha="center",
         va="center", style="normal", z=5, wrap=False, family=None, lh=1.15):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha,
            va=va, fontstyle=style, zorder=z, family=family, linespacing=lh,
            wrap=wrap)


def arrow(ax, p0, p1, *, color=MUTED, lw=2.0, style="-|>", ms=14, z=1,
          rad=0.0, ls="-"):
    a = FancyArrowPatch(
        p0, p1, arrowstyle=style, mutation_scale=ms, lw=lw, color=color,
        zorder=z, connectionstyle=f"arc3,rad={rad}", linestyle=ls,
        shrinkA=0, shrinkB=0,
    )
    ax.add_patch(a)


def chip(ax, x, y, w, h, label, *, fc, ec, tc=INK, size=8.2, weight="bold"):
    box(ax, x, y, w, h, fc=fc, ec=ec, lw=1.1, r=0.04, z=3)
    text(ax, x + w / 2, y + h / 2, label, size=size, color=tc, weight=weight)


def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200, bbox_inches="tight",
                    facecolor=PAGE)
    plt.close(fig)
    print(f"  wrote {name}.pdf / .svg / .png")


def inset(fig, x, y, w, h):
    """Add a real plotting axes at (x, y, w, h) given in the 0..100 user coords
    used by the schematic. The host ax spans figure fraction [0.01, 0.99]."""
    fl = 0.01 + (x / 100) * 0.98
    fb = 0.01 + (y / 100) * 0.98
    return fig.add_axes([fl, fb, (w / 100) * 0.98, (h / 100) * 0.98], zorder=6)


# Real signature scores (rank_score, core 189-gene tier) from
# analyses/2026-04-08-1038-n_limitation_signature_v2/results/scores_all.csv
SCORE_CONDS = ["Axenic\nRNA-seq", "Coculture\nRNA-seq", "Coculture\nproteomics"]
SCORE_VALS = [0.58, 0.05, 0.22]
SCORE_COLS = [CORAL, GREEN, GOLD]
SCORE_VERDICT = ["strong\nN-stress", "signal\nerased", "protein signal\npersists"]


def draw_score_bars(iax, *, base=9.0, verdicts=True):
    """Real bar chart of the 189-gene N-limitation signature score by condition."""
    xs = range(len(SCORE_VALS))
    iax.bar(xs, SCORE_VALS, color=SCORE_COLS, width=0.66, zorder=3,
            edgecolor="white", linewidth=1.2)
    iax.set_ylim(0, 0.86)
    iax.set_xlim(-0.6, len(SCORE_VALS) - 0.4)
    # "no signal" reference band near zero
    iax.axhspan(0, 0.07, color=MUTED, alpha=0.12, zorder=1)
    iax.text(len(SCORE_VALS) - 0.5, 0.085, "≈ no signal", ha="right",
             va="bottom", fontsize=base - 2.5, color=MUTED, style="italic")
    for xi, v, vd, c in zip(xs, SCORE_VALS, SCORE_VERDICT, SCORE_COLS):
        iax.text(xi, v + 0.02, f"{v:.2f}", ha="center", va="bottom",
                 fontsize=base, fontweight="bold", color=c)
        if verdicts:
            iax.text(xi, v + 0.115, vd, ha="center", va="bottom",
                     fontsize=base - 2, color=c, fontweight="bold")
    iax.set_xticks(list(xs))
    iax.set_xticklabels(SCORE_CONDS, fontsize=base - 0.5, color=INK)
    iax.tick_params(axis="x", length=0, pad=2)
    iax.tick_params(axis="y", labelsize=base - 2, color=MUTED, length=2)
    iax.set_yticks([0, 0.2, 0.4, 0.6])
    iax.set_ylabel("signature score", fontsize=base - 1, color=MUTED)
    for s in ("top", "right"):
        iax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        iax.spines[s].set_color(MUTED)
    iax.set_facecolor("none")


# Real analysis artifacts, read live from the committed N-limitation analysis.
ANALYSIS = OUT.parent / "analyses" / "2026-04-08-1038-n_limitation_signature_v2"
SCORES_CSV = ANALYSIS / "results" / "scores_all.csv"
SIG_CSV = ANALYSIS / "data" / "core_signature.csv"
RNA_COC_CSV = ANALYSIS / "data" / "applied_de_weissberg_rnaseq_coculture.csv"
PRO_COC_CSV = ANALYSIS / "data" / "applied_de_weissberg_proteomics_coculture.csv"

# (csv label, legend label, color, marker, linestyle, markersize)
TRAJ_STYLE = [
    ("Weissberg RNA-seq axenic",      "Axenic · RNA-seq",       CORAL,     "*", "none", 17),
    ("Weissberg RNA-seq coculture",   "Coculture · RNA-seq",    MCP,       "o", "-",     5.5),
    ("Weissberg proteomics coculture", "Coculture · proteomics", "#9A7400", "s", "-",    5.0),
    ("Weissberg proteomics axenic",   "Axenic · proteomics",    "#9AA3B2", "^", "--",    4.5),
]


def _load_target_trajectories():
    series = {}
    try:
        with open(SCORES_CSV) as f:
            for r in csv.DictReader(f):
                if r.get("tier") != "core" or r.get("role") != "target":
                    continue
                tph = (r.get("timepoint_hours") or "").strip()
                if not tph:           # skip pooled "days 60+89" summary rows
                    continue
                series.setdefault(r["label"], []).append(
                    (float(tph) / 24.0, float(r["rank_score"])))
    except FileNotFoundError:
        return {}
    for k in series:
        series[k].sort()
    return series


def draw_trajectory(iax, *, base=9.0, mini=False):
    """Real signature-score trajectory over the 90-day course, from scores_all.csv."""
    s = _load_target_trajectories()
    iax.axhspan(0.18, 0.50, color=GREEN, alpha=0.11, zorder=1)
    if not mini:
        iax.text(7, 0.255, "reference\nN-limited range", ha="left", va="bottom",
                 fontsize=base - 2.5, color="#2E8B57", style="italic",
                 linespacing=1.1)
    iax.axhline(0, color=MUTED, lw=0.7, ls=":", zorder=1)
    for key, label, col, mk, ls, ms in TRAJ_STYLE:
        pts = s.get(key, [])
        if not pts:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        iax.plot(xs, ys, color=col, marker=mk,
                 markersize=ms * (0.8 if mini else 1.0),
                 linestyle=ls, linewidth=1.3 if mini else 1.8, label=label,
                 zorder=4, markeredgecolor="white", markeredgewidth=0.5)
    iax.set_xlim(5, 95)
    iax.set_ylim(-0.05, 0.66)
    iax.set_xticks([14, 31, 60, 89])
    iax.set_yticks([0, 0.2, 0.4, 0.6])
    iax.tick_params(labelsize=base - 2, color=MUTED)
    if mini:
        iax.set_xlabel("days", fontsize=base - 1, color=MUTED, labelpad=1)
    else:
        iax.set_xlabel("time (days)", fontsize=base - 0.5, color=MUTED)
        iax.set_ylabel("signature score", fontsize=base - 0.5, color=MUTED)
    for sp in ("top", "right"):
        iax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        iax.spines[sp].set_color(MUTED)
    iax.set_facecolor("none")
    if not mini:
        leg = iax.legend(loc="upper right", fontsize=base - 2, frameon=True,
                         framealpha=0.92, edgecolor="#DDDDDD", handlelength=1.6,
                         borderpad=0.5, labelspacing=0.35)
        leg.get_frame().set_linewidth(0.6)


def _style_axes(iax, base):
    for sp in ("top", "right"):
        iax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        iax.spines[sp].set_color(MUTED)
    iax.tick_params(labelsize=base - 2, color=MUTED)
    iax.set_facecolor("none")


UP_C, DOWN_C = CORAL, MCP   # shared up/down encoding across panels


def draw_discordance(iax, *, base=9.0, mini=False):
    """Real per-gene RNA vs protein log2FC for the 147 signature genes,
    coculture day 31 (proteomics peak). Visualizes the discordance."""
    sig = {r["locus_tag"]: r["direction"]
           for r in csv.DictReader(open(SIG_CSV))}

    def day31(path):
        out = {}
        for r in csv.DictReader(open(path)):
            if r.get("timepoint") != "day 31":
                continue
            try:
                out[r["locus_tag"]] = float(r["log2fc"])
            except (ValueError, KeyError):
                pass
        return out

    rna, pro = day31(RNA_COC_CSV), day31(PRO_COC_CSV)
    for lt, d in sig.items():
        if lt not in rna or lt not in pro:
            continue
        iax.scatter(rna[lt], pro[lt], s=10 if mini else 18, alpha=0.8, zorder=3,
                    color=UP_C if d == "up" else DOWN_C,
                    edgecolor="white", linewidth=0.3)
    lim = 5.2
    iax.plot([-lim, lim], [-lim, lim], ls="--", color=MUTED, lw=0.9, zorder=1)
    iax.axhline(0, color=MUTED, lw=0.6, zorder=1)
    iax.axvline(0, color=MUTED, lw=0.6, zorder=1)
    # the discordant quadrant: protein up while RNA flat / down
    iax.axhspan(0, lim, xmin=0, xmax=0.5, color=GOLD, alpha=0.08, zorder=0)
    if not mini:
        iax.text(-lim + 0.3, lim - 0.4, "protein ↑,\nRNA flat / ↓",
                 fontsize=base - 2.5, color="#9A7400", va="top", ha="left",
                 style="italic", linespacing=1.05)
    iax.set_xlim(-lim, lim)
    iax.set_ylim(-lim, lim)
    iax.set_xticks([-4, 0, 4])
    iax.set_yticks([-4, 0, 4])
    if mini:
        iax.set_xlabel("RNA log₂FC", fontsize=base - 1, color=MUTED, labelpad=1)
        iax.set_ylabel("protein", fontsize=base - 1, color=MUTED, labelpad=1)
    else:
        iax.set_xlabel("RNA-seq log₂FC", fontsize=base - 0.5, color=MUTED)
        iax.set_ylabel("proteomics log₂FC", fontsize=base - 0.5, color=MUTED)
    _style_axes(iax, base)
    if not mini:
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], marker="o", ls="none", color=UP_C, label="up gene",
                          markeredgecolor="white", markersize=5),
                   Line2D([0], [0], marker="o", ls="none", color=DOWN_C, label="down gene",
                          markeredgecolor="white", markersize=5)]
        leg = iax.legend(handles=handles, loc="lower right", fontsize=base - 2.5,
                         frameon=True, framealpha=0.92, edgecolor="#DDDDDD")
        leg.get_frame().set_linewidth(0.6)


def draw_experiment_matrix(iax, *, base=9.0):
    """Compact view of the experiment selection: omics × condition for MED4.
    Axenic RNA-seq is single-timepoint; the rest are time courses."""
    rows = ["Proteomics", "RNA-seq"]          # bottom-up on the y-axis
    cols = ["Axenic", "Coculture"]
    # (row, col): "full" time course vs "partial" (single timepoint)
    status = {(1, 0): "partial", (1, 1): "full",
              (0, 0): "full", (0, 1): "full"}
    for (r, c), st in status.items():
        fc = KG if st == "full" else KG_L
        tc = "white" if st == "full" else KG
        iax.add_patch(Rectangle((c + 0.08, r + 0.08), 0.84, 0.84, facecolor=fc,
                                edgecolor=KG, lw=1.2, zorder=3))
        check_y = r + (0.62 if st == "partial" else 0.5)
        iax.text(c + 0.5, check_y, "✓", ha="center", va="center",
                 fontsize=base + 2, color=tc, fontweight="bold", zorder=4)
        if st == "partial":
            iax.text(c + 0.5, r + 0.27, "1 tp", ha="center", va="center",
                     fontsize=base - 3, color=KG, zorder=5, style="italic")
    iax.set_xlim(0, 2)
    iax.set_ylim(0, 2)
    iax.set_xticks([0.5, 1.5])
    iax.set_xticklabels(cols, fontsize=base - 0.5)
    iax.set_yticks([0.5, 1.5])
    iax.set_yticklabels(rows, fontsize=base - 0.5)
    iax.tick_params(length=0)
    for sp in iax.spines.values():
        sp.set_visible(False)
    iax.set_facecolor("none")


def draw_scorecard(ax, rx, ry, rw, rh, *, base=7.0):
    """Evaluate step: assess the result against the step-3 framing.
    Drawn directly in ax coords inside the region (rx, ry, rw, rh)."""
    items = [
        ("✓", GREEN, "Prediction held: RNA signature erased"),
        ("★", GOLD, "New finding: protein signature persists"),
        ("✓", GREEN, "Specific: passes negative controls"),
        ("▲", CORAL, "Caveat: axenic RNA = 1 timepoint"),
    ]
    n = len(items)
    for i, (ic, c, txt) in enumerate(items):
        yy = ry + rh - (i + 0.5) * (rh / n)
        text(ax, rx + 1.8, yy, ic, size=base + 1.5, color=c, weight="bold",
             ha="left", va="center", z=4)
        text(ax, rx + 4.8, yy, txt, size=base - 0.5, color=INK, ha="left",
             va="center", z=4)


def draw_controls(iax, *, base=9.0, mini=False):
    """Real signature score by experiment role — references score high,
    negative controls score low (specificity check). From scores_all.csv."""
    groups = {"reference": [], "negative_control": []}
    flagged = None
    for r in csv.DictReader(open(SCORES_CSV)):
        if r.get("tier") != "core":
            continue
        role = r.get("role")
        if role in groups:
            val = float(r["rank_score"])
            groups[role].append(val)
            if role == "negative_control" and val > 0.3:
                flagged = val   # high-light: a known cross-responsive control
    iax.axhspan(0.18, 0.50, color=GREEN, alpha=0.11, zorder=1)
    cols = {"reference": KG, "negative_control": "#9AA3B2"}
    if mini:
        labels = {"reference": "refs", "negative_control": "controls"}
    else:
        labels = {"reference": "reference\nN-dep studies",
                  "negative_control": "negative\ncontrols"}
    for xi, role in enumerate(("reference", "negative_control")):
        vals = groups[role]
        for j, v in enumerate(vals):
            jx = xi + ((j % 3) - 1) * 0.11
            is_flag = (role == "negative_control" and flagged is not None
                       and abs(v - flagged) < 1e-9)
            iax.scatter(jx, v, s=28 if mini else 46, zorder=4,
                        color=CORAL if is_flag else cols[role],
                        edgecolor="white", linewidth=0.6)
    if flagged is not None and not mini:
        iax.annotate("high-light\n(cross-responsive)", xy=(1, flagged),
                     xytext=(1.18, flagged + 0.04), fontsize=base - 3,
                     color=CORAL, ha="left", va="center", linespacing=1.0,
                     arrowprops=dict(arrowstyle="-", color=CORAL, lw=0.7))
    iax.axhline(0, color=MUTED, lw=0.6, ls=":", zorder=1)
    iax.set_xlim(-0.6, 1.9)
    iax.set_ylim(-0.08, 0.62)
    iax.set_xticks([0, 1])
    iax.set_xticklabels([labels["reference"], labels["negative_control"]],
                        fontsize=base - 1)
    iax.tick_params(axis="x", length=0)
    iax.set_yticks([0, 0.2, 0.4, 0.6])
    if mini:
        iax.set_ylabel("score", fontsize=base - 1, color=MUTED, labelpad=1)
    else:
        iax.set_ylabel("signature score", fontsize=base - 0.5, color=MUTED)
        iax.text(-0.5, 0.345, "N-limited\nrange", fontsize=base - 3,
                 color="#2E8B57", style="italic", va="center", linespacing=1.0)
    _style_axes(iax, base)


def draw_composition(iax, *, base=9.0):
    """The 189-gene signature: up/down split across top functional categories.
    From core_signature.csv."""
    import collections
    rows = list(csv.DictReader(open(SIG_CSV)))
    cats = collections.Counter(r["gene_category"] for r in rows)
    top = [c for c, _ in cats.most_common(7)]
    up = {c: 0 for c in top}
    dn = {c: 0 for c in top}
    for r in rows:
        c = r["gene_category"]
        if c in up:
            (up if r["direction"] == "up" else dn)[c] += 1
    order = list(reversed(top))            # largest at top of the barh
    ys = range(len(order))
    up_v = [up[c] for c in order]
    dn_v = [dn[c] for c in order]
    iax.barh(ys, dn_v, color=DOWN_C, zorder=3, height=0.7, label="down")
    iax.barh(ys, up_v, left=dn_v, color=UP_C, zorder=3, height=0.7, label="up")
    iax.set_yticks(list(ys))
    iax.set_yticklabels([c.replace(" and ", " &\n") for c in order],
                        fontsize=base - 2)
    iax.set_xlabel("genes", fontsize=base - 0.5, color=MUTED)
    iax.set_xlim(0, max(u + d for u, d in zip(up_v, dn_v)) + 3)
    iax.text(0.97, 0.06, "189 genes · 74 ↑ / 115 ↓", transform=iax.transAxes,
             ha="right", va="bottom", fontsize=base - 1.5, color=INK,
             fontweight="bold")
    _style_axes(iax, base)
    leg = iax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.13),
                     fontsize=base - 2.5, frameon=True, framealpha=0.92,
                     edgecolor="#DDDDDD", ncol=2, handlelength=1.2,
                     columnspacing=1.0)
    leg.get_frame().set_linewidth(0.6)


def kg_node(ax, cx, cy, w, h, title, sub, fc, ec, *, tsize=8.0, ssize=5.8):
    box(ax, cx - w / 2, cy - h / 2, w, h, fc=fc, ec=ec, lw=1.3, r=0.08, z=4)
    if sub:
        text(ax, cx, cy + h * 0.16, title, size=tsize, weight="bold", z=5)
        text(ax, cx, cy - h * 0.22, sub, size=ssize, color=MUTED, z=5)
    else:
        text(ax, cx, cy, title, size=tsize, weight="bold", z=5)


def kg_edge(ax, p0, p1, label, *, color=MUTED):
    arrow(ax, p0, p1, color=color, lw=1.3, ms=9, z=2)
    mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
    text(ax, mx, my + 0.9, label, size=5.6, color=color, style="italic", z=3)


def draw_kg_model(ax, ox, oy):
    """Native, on-brand recreation of the KG data model (Gene-centered hub).
    ox, oy = center of the Gene node. Adapted from
    multiomics_biocypher_kg/docs/kg_concept_figure.svg."""
    gw, gh = 13, 6.6
    # satellites: (dx, dy, w, h, title, sub, fc, ec)
    sats = [
        (-18, 7.0, 14, 5.2, "Publication", "studies", MCP_L, MCP),
        (0, 8.0, 15, 5.2, "Experiment", "RNA-seq · proteomics", MCP_L, MCP),
        (-20, 0, 13, 5.2, "Organism", "strain / taxonomy", GOLD_L, GOLD),
        (20, 0, 15, 5.2, "Ortholog group", "cross-strain", GOLD_L, GOLD),
        (-15, -8.5, 14, 5.2, "Ontology", "GO · KEGG · Pfam", RES_L, RES),
        (14, -8.5, 14, 5.2, "Metabolite", "reactions · transport", "#E2F1E9", GREEN),
    ]
    pos = {}
    for dx, dy, w, h, title, sub, fc, ec in sats:
        cx, cy = ox + dx, oy + dy
        pos[title] = (cx, cy, w, h)
    # edges (drawn under nodes)
    g = (ox, oy)
    kg_edge(ax, pos["Publication"][:2], pos["Experiment"][:2], "reports")
    kg_edge(ax, (pos["Experiment"][0], pos["Experiment"][1] - 2.7),
            (ox, oy + gh / 2), "measures")
    kg_edge(ax, (ox - gw / 2, oy), (pos["Organism"][0] + 6.5, pos["Organism"][1]), "in")
    kg_edge(ax, (ox + gw / 2, oy), (pos["Ortholog group"][0] - 7.5, pos["Ortholog group"][1]), "conserved")
    kg_edge(ax, (ox - gw / 4, oy - gh / 2), pos["Ontology"][:2], "annotated")
    kg_edge(ax, (ox + gw / 4, oy - gh / 2), pos["Metabolite"][:2], "acts on")
    # satellite nodes
    for dx, dy, w, h, title, sub, fc, ec in sats:
        kg_node(ax, ox + dx, oy + dy, w, h, title, sub, fc, ec)
    # Gene hub (drawn last, emphasized)
    box(ax, ox - gw / 2, oy - gh / 2, gw, gh, fc=KG, ec=KG, lw=1.5, r=0.08, z=6)
    text(ax, ox, oy + 1.0, "Gene", size=11, weight="bold", color="white", z=7)
    text(ax, ox, oy - 1.6, "central entity", size=6.0, color="#CDE7E7", z=7)


# --------------------------------------------------------------------------- #
# FIGURE 1 — System architecture (3-layer stack)
# --------------------------------------------------------------------------- #
def fig_architecture():
    fig, ax = new_canvas(12.4, 9.6)

    L, R = 5, 81          # content band for the three layer boxes
    W = R - L
    cx_mid = (L + R) / 2

    text(ax, 50, 97.5, "An AI research assistant for multi-omics discovery in marine microbes",
         size=16, weight="bold")
    text(ax, 50, 93.8,
         "Integrated published data  →  ask questions in plain language  →  reproducible analyses",
         size=11, color=MUTED, style="italic")

    # ---- Data types feeding the KG ---------------------------------------- #
    src_y, src_h = 87.0, 5.6
    text(ax, L, 91.6, "Published multi-omics data", size=10, weight="bold",
         color=KG, ha="left")
    sources = ["Genomes", "Transcriptomes", "Proteomes", "Metabolomes",
               "Annotations\n& orthologs"]
    sw = 12.5
    gap = (W - sw) / (len(sources) - 1)
    for i, s in enumerate(sources):
        sx = L + i * gap
        chip(ax, sx, src_y, sw, src_h, s, fc=KG_L, ec=KG, tc=INK, size=8.2,
             weight="normal")
        arrow(ax, (sx + sw / 2, src_y), (sx + sw / 2, 85.5), color=KG, lw=1.3,
              ms=9)

    # ---- Layer 1: the Knowledge Graph (with native data-model diagram) ---- #
    box(ax, L, 52.0, W, 33.0, fc=KG_L, ec=KG, lw=2.0, r=0.02)
    text(ax, L + 3.5, 82.6, "1   KNOWLEDGE GRAPH", size=13, weight="bold",
         color=KG, ha="left")
    text(ax, L + 3.5, 79.6,
         "Prochlorococcus & Alteromonas — published multi-omics, integrated into one graph",
         size=9.5, color=INK, ha="left")

    # data model: Gene-centered hub on the left
    draw_kg_model(ax, ox=29, oy=64.0)

    # divider + headline counts on the right
    ax.plot([60, 60], [55, 75.5], color=KG, lw=0.8, alpha=0.4, zorder=2)
    text(ax, 71, 76.0, "at a glance", size=8.5, weight="bold", color=KG,
         style="italic")
    kg_stats = [("120,000+", "genes"), ("43", "published\nstudies"),
                ("197", "experiments"), ("45", "organisms")]
    sx_cols, sy_rows = [66.0, 76.5], [70.0, 60.0]
    for i, (val, lab) in enumerate(kg_stats):
        cx = sx_cols[i % 2]
        cy = sy_rows[i // 2]
        text(ax, cx, cy, val, size=14.5, weight="bold", color=KG)
        text(ax, cx, cy - 3.4, lab, size=8.0, color=MUTED)

    arrow(ax, (cx_mid, 52.0), (cx_mid, 49.0), color=INK, lw=2.6, ms=18)

    # ---- Layer 2: ask questions ------------------------------------------- #
    box(ax, L, 30.0, W, 19.0, fc=MCP_L, ec=MCP, lw=2.0, r=0.02)
    text(ax, L + 3.5, 45.6, "2   ASK QUESTIONS", size=13, weight="bold",
         color=MCP, ha="left")
    text(ax, L + 3.5, 42.0,
         "In plain language or in Python — no query language to learn",
         size=9.5, color=INK, ha="left")

    qtypes = ["Expression\n& DE", "Function\n& pathways",
              "Orthologs &\ncross-organism", "Metabolism\n& chemistry",
              "Sequence &\ngenome context"]
    qw, qh = 13.5, 7.4
    qgap = (W - qw) / (len(qtypes) - 1)
    for i, q in enumerate(qtypes):
        qx = L + i * qgap
        chip(ax, qx, 31.6, qw, qh, q, fc="#FFFFFF", ec=MCP, tc=MCP, size=8.4,
             weight="bold")

    arrow(ax, (cx_mid, 30.0), (cx_mid, 27.0), color=INK, lw=2.6, ms=18)

    # ---- Layer 3: AI research assistant ----------------------------------- #
    box(ax, L, 4.0, W, 22.0, fc=RES_L, ec=RES, lw=2.0, r=0.02)
    text(ax, L + 3.5, 22.6, "3   AI RESEARCH ASSISTANT", size=13,
         weight="bold", color=RES, ha="left")
    text(ax, L + 3.5, 19.0,
         "Plans the analysis, writes and runs the code, and reports results with caveats",
         size=9.5, color=INK, ha="left")

    steps = ["Plans the\nanalysis", "Runs reproducible\ncode",
             "Reports results\n+ caveats"]
    aw = 22.5
    agap = (W - aw) / (len(steps) - 1)
    for i, head in enumerate(steps):
        axp = L + i * agap
        box(ax, axp, 6.0, aw, 8.4, fc="#FFFFFF", ec=RES, lw=1.3, r=0.04, z=3)
        text(ax, axp + aw / 2, 10.2, head, size=9.2, weight="bold", color=RES)
        if i < len(steps) - 1:
            arrow(ax, (axp + aw, 10.2), (L + (i + 1) * agap, 10.2),
                  color=MUTED, lw=1.6, ms=11)

    # ---- Researcher (right rail spanning the full stack) ----------------- #
    rx, rw = 84.0, 13.0
    box(ax, rx, 4.0, rw, 81.0, fc=GOLD_L, ec=GOLD, lw=2.0, r=0.04)
    ax.text(rx + rw / 2, 56.0, "RESEARCHER", rotation=90, ha="center",
            va="center", fontsize=14, fontweight="bold", color="#9A7400",
            zorder=5)
    # bidirectional exchange with the assistant layer
    arrow(ax, (rx, 18.0), (R, 18.0), color=GOLD, lw=2.2, ms=14)
    text(ax, (R + rx) / 2, 19.7, "question", size=7.4, color="#9A7400")
    arrow(ax, (R, 10.5), (rx, 10.5), color=RES, lw=2.2, ms=14)
    text(ax, (R + rx) / 2, 12.2, "answer", size=7.4, color=RES)

    save(fig, "fig1_architecture")


# --------------------------------------------------------------------------- #
# FIGURE 2 — Research workflow (6-step methodology)
# --------------------------------------------------------------------------- #
def fig_workflow():
    fig, ax = new_canvas(12.0, 6.2)

    text(ax, 50, 95, "From question to finished result, in six steps",
         size=16, weight="bold")
    text(ax, 50, 88.5,
         "A researcher reviews and approves each step before the analysis moves on",
         size=11, color=MUTED, style="italic")

    steps = [
        ("1", "Research\nquestion", "Clarify and\nlock the question", KG),
        ("2", "Find the\ndata", "Relevant studies,\nexperiments, organisms", KG),
        ("3", "Frame the\nanalysis", "Hypothesis, controls,\nexpected outcome", KG),
        ("4", "Build the\nmethod", "Analysis code from\na worked example", MCP),
        ("5", "Run &\nanalyze", "Scored results,\nfigures, tables", MCP),
        ("6", "Evaluate", "Assess, harvest\ncaveats, finalize", MCP),
    ]
    sw, sh = 13.5, 23.0
    y = 50.0
    margin = 3.5
    sgap = (100 - 2 * margin - sw) / (len(steps) - 1)
    for i, (num, title, body, col) in enumerate(steps):
        x = margin + i * sgap
        cx = x + sw / 2
        fc = KG_L if col == KG else MCP_L
        box(ax, x, y, sw, sh, fc=fc, ec=col, lw=2.0, r=0.05, z=3)
        # number badge
        box(ax, cx - 3.0, y + sh - 6.4, 6.0, 6.0, fc=col, ec=col, r=0.5, z=4)
        text(ax, cx, y + sh - 3.4, num, size=14, weight="bold", color="white",
             z=5)
        text(ax, cx, y + 11.0, title, size=11, weight="bold", color=INK)
        text(ax, cx, y + 4.8, body, size=8.2, color=MUTED)
        if i < len(steps) - 1:
            nx = margin + (i + 1) * sgap
            arrow(ax, (x + sw, y + sh / 2), (nx, y + sh / 2), color=MUTED,
                  lw=2.2, ms=14)

    # ---- phase brackets: proposal (1-3) vs execute (4-6) ------------------ #
    def bracket(x0, x1, yb, label, color):
        ax.plot([x0, x0, x1, x1], [yb + 1.6, yb, yb, yb + 1.6], color=color,
                lw=1.8, zorder=2)
        text(ax, (x0 + x1) / 2, yb - 3.2, label, size=10.5, weight="bold",
             color=color)

    bracket(margin, margin + 2 * sgap + sw, 46.0, "PLAN THE STUDY", KG)
    bracket(margin + 3 * sgap, 100 - margin, 46.0, "RUN IT", MCP)

    # ---- light rigor note ------------------------------------------------- #
    box(ax, 10, 6.0, 80, 26.0, fc=GOLD_L, ec=GOLD, lw=2.0, r=0.03)
    text(ax, 50, 27.5, "Built for rigor and reproducibility", size=12,
         weight="bold", color="#9A7400")
    notes = [
        "The researcher reviews and approves each step before it continues",
        "Hypotheses, controls and caveats are written down as the data demands them",
        "Every result is scripted and reproducible from start to finish",
    ]
    for i, it in enumerate(notes):
        gy = 21.0 - i * 4.6
        text(ax, 17, gy, "✓", size=12, weight="bold", color=GREEN, ha="left")
        text(ax, 20.5, gy, it, size=10, color=INK, ha="left")

    save(fig, "fig2_workflow")


# --------------------------------------------------------------------------- #
# FIGURE 3 — Science exemplar / impact
# --------------------------------------------------------------------------- #
def fig_exemplar():
    fig, ax = new_canvas(12.0, 7.6)

    text(ax, 50, 96.5,
         "Exemplar: how a heterotroph rescues Prochlorococcus from nitrogen starvation",
         size=14.5, weight="bold")
    text(ax, 50, 92.3,
         "KG-driven re-analysis reproduces and extends Weissberg et al. 2025 — a 90-day "
         "axenic-vs-coculture multi-omics time course",
         size=9.5, color=MUTED, style="italic")

    # ---- the biological system (left) ------------------------------------ #
    box(ax, 4, 60.0, 44, 26.0, fc=GREY_L, ec=MUTED, lw=1.3, r=0.03)
    text(ax, 6, 83.0, "The system", size=10, weight="bold", color=INK, ha="left")
    text(ax, 6, 79.2,
         "Prochlorococcus MED4 in N-poor seawater (PRO99-lowN),\n"
         "continuous light, 24 °C — grown two ways:",
         size=8.4, color=INK, ha="left", va="center")

    # axenic — dies
    box(ax, 6, 67.5, 19, 8.5, fc="#FBE7E2", ec=CORAL, lw=1.5, r=0.05)
    text(ax, 15.5, 73.4, "Axenic", size=9.2, weight="bold", color=CORAL)
    text(ax, 15.5, 69.8, "dies in ~2 weeks", size=7.8, color=CORAL)
    # coculture — survives
    box(ax, 27, 67.5, 19, 8.5, fc="#E2F1E9", ec=GREEN, lw=1.5, r=0.05)
    text(ax, 36.5, 73.4, "+ Alteromonas", size=8.8, weight="bold", color=GREEN)
    text(ax, 36.5, 69.8, "survives 90+ days", size=7.8, color=GREEN)
    text(ax, 26, 63.2,
         "Same genotype, same medium — a heterotrophic partner is the only difference.",
         size=7.6, color=MUTED, ha="center", style="italic")

    # ---- the KG analysis pipeline (right) -------------------------------- #
    box(ax, 52, 60.0, 44, 26.0, fc=MCP_L, ec=MCP, lw=1.3, r=0.03)
    text(ax, 54, 82.5, "The analysis, built on the KG", size=10.5,
         weight="bold", color=MCP, ha="left")
    pipe = [
        "Pulled RNA-seq + proteomics for MED4 in both conditions",
        "Built a 189-gene nitrogen-limitation signature from two\nindependent reference studies",
        "Scored each condition and compared RNA vs protein on\ngenes seen by both",
    ]
    for i, p in enumerate(pipe):
        py = 77.0 - i * 6.0
        text(ax, 55, py, "•", size=11, color=MCP, ha="left", weight="bold")
        text(ax, 57, py, p, size=8.4, color=INK, ha="left", va="center")

    arrow(ax, (50, 60.0), (50, 55.5), color=INK, lw=2.4, ms=15)

    # ---- result: real embedded trajectory plot --------------------------- #
    text(ax, 50, 53.5,
         "Result: 189-gene N-stress signature score over the 90-day course",
         size=10.5, weight="bold")
    box(ax, 17, 24.5, 66, 27.0, fc="#FFFFFF", ec=MUTED, lw=1.0, r=0.02, z=2)
    iax = inset(fig, 23, 28.5, 55, 19.5)
    draw_trajectory(iax, base=8.5)

    # ---- takeaway -------------------------------------------------------- #
    box(ax, 8, 6.0, 84, 16.5, fc=RES_L, ec=RES, lw=1.8, r=0.03)
    text(ax, 50, 18.8, "The finding", size=10.5, weight="bold", color=RES)
    text(ax, 50, 13.8,
         "Coculture abolishes the transcriptional nitrogen-stress response, yet a protein-level\n"
         "signature persists — a genuine RNA/protein discordance (not a coverage artifact).",
         size=9.4, color=INK, va="center")
    text(ax, 50, 8.4,
         "The KG + agent reproduced the published rescue phenotype and surfaced a new "
         "post-transcriptional layer — in days, scripted and reproducible, with caveats attached.",
         size=8.6, color=RES, weight="bold", va="center")

    save(fig, "fig3_exemplar")


# --------------------------------------------------------------------------- #
# FIGURE 3a — science exemplar, four real embedded data panels
# --------------------------------------------------------------------------- #
def _panel(ax, x, y, w, h, tag, title, color):
    """Draw a panel frame with a tag badge + title; return inset region coords."""
    box(ax, x, y, w, h, fc="#FFFFFF", ec=color, lw=1.5, r=0.02, z=2)
    box(ax, x + 1.6, y + h - 5.4, 5.2, 4.0, fc=color, ec=color, r=0.4, z=3)
    text(ax, x + 4.2, y + h - 3.4, tag, size=10, weight="bold", color="white",
         z=4)
    text(ax, x + 8.4, y + h - 3.4, title, size=9.2, weight="bold", color=INK,
         ha="left", z=4)
    return (x + 6.5, y + 5.0, w - 11, h - 13.5)   # inset: l, b, w, h


def fig_exemplar_multi():
    fig, ax = new_canvas(13.6, 8.8)

    text(ax, 50, 97.2,
         "How a heterotroph rescues Prochlorococcus from nitrogen starvation",
         size=15.5, weight="bold")
    text(ax, 50, 93.4,
         "MED4 dies axenically under N-starvation but survives 90+ days with Alteromonas — "
         "KG re-analysis of Weissberg et al. 2025 shows why",
         size=9.8, color=MUTED, style="italic")

    PW, PH = 44.5, 38.0
    xL, xR = 4.0, 51.5
    yT, yB = 48.5, 7.5

    # Panel A — trajectory
    l, b, w, h = _panel(ax, xL, yT, PW, PH, "A",
                        "Signature score over 90 days", MCP)
    draw_trajectory(inset(fig, l, b, w, h), base=8.0)

    # Panel B — RNA vs protein discordance
    l, b, w, h = _panel(ax, xR, yT, PW, PH, "B",
                        "RNA vs protein, day 31 (147 genes)", GOLD)
    draw_discordance(inset(fig, l, b, w, h), base=8.0)

    # Panel C — method validation / controls
    l, b, w, h = _panel(ax, xL, yB, PW, PH, "C",
                        "Specificity: references vs controls", KG)
    draw_controls(inset(fig, l, b, w, h), base=8.0)

    # Panel D — signature composition
    l, b, w, h = _panel(ax, xR, yB, PW, PH, "D",
                        "What the 189-gene signature is", RES)
    draw_composition(inset(fig, l, b, w, h), base=8.0)

    # finding strip
    text(ax, 50, 3.0,
         "Coculture erases the transcriptional N-stress signature (A), yet a protein-level "
         "signature persists (A, B) — a genuine RNA/protein discordance. The signature is "
         "specific (C) and dominated by stress, translation & photosynthesis genes (D).",
         size=8.4, color=INK, weight="bold")

    save(fig, "fig3a_exemplar_panels")


# --------------------------------------------------------------------------- #
# FIGURE 2a — the 6-step process, worked through one real study
# --------------------------------------------------------------------------- #
def fig_workflow_example():
    fig, ax = new_canvas(12.6, 8.2)

    text(ax, 50, 96.8, "The six steps, worked through one study",
         size=16, weight="bold")
    text(ax, 50, 92.8,
         "Each step produces a concrete, inspectable artifact — shown here for the "
         "Alteromonas–Prochlorococcus N-stress question (Weissberg 2025)",
         size=10, color=MUTED, style="italic")

    # kind: "text" → body string; otherwise a draw_* thumbnail + caption
    steps = [
        ("1", "Research question", KG, "text",
         "Does coculture with\nAlteromonas lower the\nN-stress response of\nProchlorococcus MED4?"),
        ("2", "Find the data", KG, "text",
         "Weissberg 2025 — MED4\nacross RNA-seq + proteomics,\naxenic and coculture"),
        ("3", "Frame the analysis", KG, "experiments",
         "select experiments · define\nthe 189-gene signature"),
        ("4", "Build the method", MCP, "controls",
         "QC: the score separates\nN-stress from controls"),
        ("5", "Run & analyze", MCP, "trajectory",
         "signature score by condition,\nover the 90-day course"),
        ("6", "Evaluate", MCP, "scorecard",
         "assess vs the prediction;\nharvest caveats"),
    ]

    bw, bh = 29.5, 31.0
    xs = [4, 35.4, 66.8]
    ys = [49.0, 7.5]
    positions = [(xs[c], ys[r]) for r in range(2) for c in range(3)]

    for idx, (num, title, col, kind, body) in enumerate(steps):
        x, y = positions[idx]
        cx = x + bw / 2
        fc = KG_L if col == KG else MCP_L
        box(ax, x, y, bw, bh, fc=fc, ec=col, lw=2.0, r=0.04, z=3)
        box(ax, x + 2.0, y + bh - 7.0, 5.6, 5.6, fc=col, ec=col, r=0.5, z=4)
        text(ax, x + 4.8, y + bh - 4.2, num, size=12.5, weight="bold",
             color="white", z=5)
        text(ax, x + 10.5, y + bh - 4.2, title, size=10.5, weight="bold",
             color=INK, ha="left")
        ax.plot([x + 2.5, x + bw - 2.5], [y + bh - 8.4, y + bh - 8.4],
                color=col, lw=0.8, alpha=0.5, zorder=4)
        if kind == "text":
            text(ax, cx, y + (bh - 8.4) / 2 + 0.5, body, size=9.0, color=INK,
                 va="center")
            continue
        # caption + embedded thumbnail
        text(ax, cx, y + bh - 11.0, body, size=7.2, color=MUTED, va="center")
        # white plot backing so the thumbnail reads on the tinted card
        bxx, bxy, bxw, bxh = x + 2.6, y + 2.4, bw - 5.2, 12.6
        box(ax, bxx, bxy, bxw, bxh, fc="#FFFFFF", ec=col, lw=0.8, r=0.04, z=3)
        if kind == "scorecard":
            draw_scorecard(ax, bxx, bxy, bxw, bxh, base=7.0)
            continue
        iax = inset(fig, x + 5.5, y + 4.2, bw - 9.5, 8.6)
        if kind == "experiments":
            draw_experiment_matrix(iax, base=6.6)
        elif kind == "controls":
            draw_controls(iax, base=6.6, mini=True)
        elif kind == "trajectory":
            draw_trajectory(iax, base=6.6, mini=True)
        elif kind == "discordance":
            draw_discordance(iax, base=6.6, mini=True)

    # flow arrows: 1→2→3 (top, L→R), down to 4, 4→5→6 wait order is row2 L→R = 4,5,6
    a = MUTED
    # top row left→right
    arrow(ax, (xs[0] + bw, ys[0] + bh / 2), (xs[1], ys[0] + bh / 2), color=a, lw=2.2, ms=14)
    arrow(ax, (xs[1] + bw, ys[0] + bh / 2), (xs[2], ys[0] + bh / 2), color=a, lw=2.2, ms=14)
    # wrap down from step 3 (top-right) to step 4 (bottom-left)
    arrow(ax, (xs[2] + bw / 2, ys[0]), (xs[0] + bw / 2, ys[1] + bh),
          color=a, lw=2.0, ms=13, rad=-0.18)
    # bottom row left→right
    arrow(ax, (xs[0] + bw, ys[1] + bh / 2), (xs[1], ys[1] + bh / 2), color=a, lw=2.2, ms=14)
    arrow(ax, (xs[1] + bw, ys[1] + bh / 2), (xs[2], ys[1] + bh / 2), color=a, lw=2.2, ms=14)

    save(fig, "fig2a_workflow_example")


if __name__ == "__main__":
    print("Generating grant figures →", OUT)
    fig_architecture()
    fig_workflow()
    fig_workflow_example()
    fig_exemplar()
    fig_exemplar_multi()
    print("Done.")
