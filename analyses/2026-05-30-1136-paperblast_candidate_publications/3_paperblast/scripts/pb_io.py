"""pb_io.py — shared loader for saved PaperBLAST pages.

Handles both save formats the researcher used:
  - "Webpage, HTML Only"  -> .html (plain)
  - "Webpage, Single File" -> .mhtml (MIME multipart, quoted-printable HTML)
and recursive subfolders (results/<pool>/<seed>.{html,mhtml}).

Use load_pages(RESULTS) to iterate (path, html_text) over every saved page,
skipping the _quarantine/ folder.
"""
from __future__ import annotations

import quopri
from email import policy
from email.parser import BytesParser
from pathlib import Path


def _read_mhtml(path: Path) -> str:
    """Return the decoded text/html body of an MHTML (single-file) save."""
    msg = BytesParser(policy=policy.default).parsebytes(path.read_bytes())
    for part in msg.walk():
        if part.get_content_type() == "text/html":
            payload = part.get_payload(decode=True)
            if payload:
                return payload.decode(part.get_content_charset() or "utf-8", "replace")
    # fallback: decode whole file as quoted-printable
    return quopri.decodestring(path.read_bytes()).decode("utf-8", "replace")


def read_page(path: Path) -> str:
    if path.suffix.lower() == ".mhtml":
        return _read_mhtml(path)
    return path.read_text(encoding="utf-8", errors="replace")


def load_pages(results_dir: Path):
    """Yield (path, html_text) for every saved page under results_dir
    (recursive), skipping _quarantine/. Sorted for deterministic order."""
    for p in sorted(results_dir.rglob("*")):
        if p.suffix.lower() not in (".html", ".mhtml"):
            continue
        if "_quarantine" in p.parts:
            continue
        yield p, read_page(p)
