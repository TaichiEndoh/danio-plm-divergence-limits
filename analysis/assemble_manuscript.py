#!/usr/bin/env python3
"""
Assemble the per-section drafts in docs/ into one manuscript file.

The section files are the source of truth; this script only concatenates them, so it can be
re-run after any edit. Output: docs/MANUSCRIPT_FULL_en.md — the single file to hand to the
English editor or to paste into a journal template.
"""
from __future__ import annotations

import datetime
import re
from pathlib import Path

DOCS = Path("docs")


# Headings that are working notes for the authors, not manuscript content. Everything from
# such a heading up to the next heading of the same level is dropped from the assembled file;
# it stays in the section file, which remains the working document.
INTERNAL_SECTIONS = {
    "Figure production checklist",
    "Open decisions for the corresponding author",
    "Currently uncited entries",
    "Placement of the authors' own work",
}


# A blockquote block opening with one of these is a working note for the authors, not
# manuscript content, and is dropped wherever it appears. Figure legends are also
# blockquotes, so the test has to be on the opening marker, never on the quoting alone.
NOTE_OPENER = re.compile(r"^>\s*\*\*(\[internal\]|Structure note|Revision note|Note for)")


def body(path: Path) -> str:
    """Section text, minus the H1 title, the leading editorial note-block, any blockquote
    marked as a working note, and any author-facing section listed in INTERNAL_SECTIONS."""
    lines = path.read_text(encoding="utf-8").split("\n")
    out: list[str] = []
    in_lead = True          # still above the section's first real content
    skip_until_level = None
    i = 0
    while i < len(lines):
        line = lines[i]

        m = re.match(r"^(#{2,6})\s+(.*)$", line)
        if m:
            level, text = len(m.group(1)), m.group(2).strip()
            if skip_until_level is not None and level <= skip_until_level:
                skip_until_level = None
            if text in INTERNAL_SECTIONS:
                skip_until_level = level
                i += 1
                continue
        if skip_until_level is not None or line.startswith("# "):
            i += 1
            continue

        if line.startswith(">"):
            block = []
            while i < len(lines) and lines[i].startswith(">"):
                block.append(lines[i])
                i += 1
            if in_lead or NOTE_OPENER.match(block[0]):
                if i < len(lines) and not lines[i].strip():
                    i += 1          # and the blank line the dropped block leaves behind
                continue
            out.extend(block)
            in_lead = False
            continue

        if line.strip():
            in_lead = False
        out.append(line)
        i += 1

    # collapse the trailing separator a dropped section may leave behind
    text = "\n".join(out).strip()
    return re.sub(r"\n-{3,}\s*$", "", text).strip()


def main() -> None:
    ab = (DOCS / "ABSTRACT_en.md").read_text(encoding="utf-8")
    abstract = re.search(r"## Structured version.*?\n\n(.*?)\n---", ab, re.S).group(1).strip()
    keywords = re.search(r"## Keywords\n\n(.*?)\n\n---", ab, re.S).group(1).strip()
    title = re.search(r"^1\. \*\*(.*?)\*\*", ab, re.M | re.S).group(1).replace("\n   ", " ").strip()

    doc = f"""# {title}

**Authors.** Taichi Endoh, Gerry Amor Camer, Kotetsu Kayama, Daiji Endoh, Hiroki Teraoka
*(author list and order to be confirmed)*

**Target journal.** *Proteomes* (MDPI) — the follow-up to Endoh et al. 2026, *Proteomes* 14(3):36.

**Draft assembled.** {datetime.date.today().isoformat()} · Repository:
https://github.com/TaichiEndoh/zebrafish-aescallii-esm3-alphafold

> **Note for the English editor.** This is a complete first draft, assembled from the section
> files in `docs/`. Every numeric value was recomputed from the committed analysis outputs in
> `reports/` during a pre-submission audit, and all 27 references were verified against
> publisher records. Bracketed numbers such as [ref 17] point to the reference list at the end.
> Figure legends are in the Figures section; the figure files are in `reports/figures/final/`
> at 300 dpi. Regenerate this file with `python analysis/assemble_manuscript.py`.

---

## Abstract

{abstract}

**Keywords:** {keywords}

---

## 1. Introduction

{body(DOCS / 'INTRODUCTION_en.md')}

---

## 2. Materials and Methods

{body(DOCS / 'METHODS_en.md')}

---

## 3. Results

{body(DOCS / 'RESULTS_en.md')}

---

## 4. Discussion

{body(DOCS / 'DISCUSSION_en.md')}

---

## 5. Conclusions

{body(DOCS / 'CONCLUSIONS_en.md')}

---

{body(DOCS / 'BACKMATTER_en.md')}

---

## Figures and Tables

{body(DOCS / 'FIGURES_en.md')}

---

## References

{body(DOCS / 'REFERENCES_en.md')}
"""
    out = DOCS / "MANUSCRIPT_FULL_en.md"
    out.write_text(doc, encoding="utf-8")
    print(f"Wrote {out}")
    print(f"  characters {len(doc):,}  words ~{len(doc.split()):,}")
    print(f"  unresolved [CITE: markers: {doc.count('[CITE:')}")
    # [confirm] flags a decision the corresponding author still owes; they are useful to
    # co-authors but must not survive to submission, so count them on every build.
    confirms = doc.count("[confirm]")
    print(f"  unresolved [confirm] markers: {confirms}"
          + ("  <- clear these before submitting" if confirms else ""))
    leaked = doc.count("[internal]")
    if leaked:
        raise SystemExit(f"ERROR: {leaked} internal note block(s) leaked into the manuscript")


if __name__ == "__main__":
    main()
