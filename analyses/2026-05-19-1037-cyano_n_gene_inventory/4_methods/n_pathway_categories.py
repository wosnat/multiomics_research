"""Map an ortholog group (by consensus gene name) to a coarse N pathway category.

Categories used in the inventory heatmap as the column annotation strip:
- regulation
- assimilation_core
- ammonium_uptake
- urea_uptake
- urease
- cyanate
- nitrate_nitrite
- mo_cofactor
- met_biosynthesis
- other

Order matters: when displayed, categories are color-coded in this canonical order.
"""
from __future__ import annotations

CATEGORY_ORDER = [
    "regulation",
    "assimilation_core",
    "ammonium_uptake",
    "urea_uptake",
    "urease",
    "cyanate",
    "nitrate_nitrite",
    "mo_cofactor",
    "met_biosynthesis",
    "other",
]

# Gene-name → category. Case-insensitive exact match.
_GENE_TO_CATEGORY = {
    # regulation
    "ntcA": "regulation",
    "glnB": "regulation",
    "ntrB": "regulation",
    "ntrX": "regulation",
    "ptrA": "regulation",            # Crp/FNR paralog of NtcA in MIT9303
    "merR": "regulation",            # LysR/MerR family transcriptional regulator
    # assimilation core
    "glnA": "assimilation_core",
    "glnN1": "assimilation_core",    # alternative glutamine synthetase
    "glsF": "assimilation_core",
    "carA": "assimilation_core",
    # ammonium uptake
    "amt1": "ammonium_uptake",
    "amt2": "ammonium_uptake",
    # urea uptake
    "urtA": "urea_uptake",
    "urtB": "urea_uptake",
    "urtC": "urea_uptake",
    "urtD": "urea_uptake",
    "urtE": "urea_uptake",
    # urease
    "ureA": "urease",
    "ureB": "urease",
    "ureC": "urease",
    "ureD": "urease",
    "ureE": "urease",
    "ureF": "urease",
    "ureG": "urease",
    # cyanate
    "cynA": "cyanate",
    "cynB": "cyanate",
    "cynD": "cyanate",
    "cynS": "cyanate",
    "cynH": "cyanate",
    # nitrate / nitrite
    "nrtP": "nitrate_nitrite",
    "narB": "nitrate_nitrite",
    "narM": "nitrate_nitrite",
    "nirA": "nitrate_nitrite",
    "focA": "nitrate_nitrite",       # nitrite transporter FNT family
    # Mo cofactor biosynthesis
    "moaA": "mo_cofactor",
    "moaB": "mo_cofactor",
    "moaC": "mo_cofactor",
    "moaD": "mo_cofactor",
    "moaE": "mo_cofactor",
    "moeA": "mo_cofactor",
    "moeB": "mo_cofactor",
    "mobA": "mo_cofactor",
    # Met / sulfur amino acid biosynthesis
    "cysK": "met_biosynthesis",
    "metB": "met_biosynthesis",
    "metC": "met_biosynthesis",
    "metY": "met_biosynthesis",
    "MetB-like": "met_biosynthesis",
}


def categorize(gene_name) -> str:
    """Return the pathway category for a consensus_gene_name. NaN or unknown → 'other'."""
    if gene_name is None:
        return "other"
    try:
        if str(gene_name) == "nan":
            return "other"
    except Exception:
        return "other"
    name = str(gene_name).strip()
    if not name:
        return "other"
    return _GENE_TO_CATEGORY.get(name, "other")


# Palette: distinguishable colors in CATEGORY_ORDER order (tab10 + custom for "other")
CATEGORY_PALETTE = {
    "regulation":        "#1f77b4",   # blue
    "assimilation_core": "#2ca02c",   # green
    "ammonium_uptake":   "#17becf",   # cyan
    "urea_uptake":       "#ff7f0e",   # orange
    "urease":            "#d62728",   # red
    "cyanate":           "#9467bd",   # purple
    "nitrate_nitrite":   "#8c564b",   # brown
    "mo_cofactor":       "#e377c2",   # pink
    "met_biosynthesis":  "#bcbd22",   # olive
    "other":             "#cccccc",   # gray
}
