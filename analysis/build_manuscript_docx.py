#!/usr/bin/env python3
"""
Word (.docx) build of the full manuscript, with the composed figures embedded inline.

The Markdown and PDF copies do not put the figures where the reader meets the text, which
makes them of little use to an editor working through the manuscript. This produces the file
that can actually be marked up: one document, figures sitting directly above their legends,
tracked-changes ready.

The Figures section of MANUSCRIPT_FULL_en.md names the *source panel* files that went into
each composed figure, which is what a reader reproducing the analysis needs but not what an
editor needs. Here those provenance lines are replaced by the composed figure itself; the
provenance survives in docs/FIGURES_en.md and in reports/figures/final/README.md.

Requires pandoc. Writes docs/MANUSCRIPT_FULL_en.docx
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

SRC = Path("docs/MANUSCRIPT_FULL_en.md")
OUT = Path("docs/MANUSCRIPT_FULL_en.docx")
FIGDIR = Path("reports/figures/final")

# manuscript figure number -> composed file in reports/figures/final/
FIGURES = {
    "1": "Figure1_study_design.png",
    "2": "Figure2_scale_robustness.png",
    "3": "Figure3_prioritization.png",
    "4": "Figure4_hotspots_and_axes.png",
    "5": "Figure5_structural_grounding.png",
    "6": "Figure6_feature_annotation.png",
    "S1": "FigureS1_rmsd_control.png",
    "S2": "FigureS2_variant_rediscovery.png",
    "S3": "FigureS3_baseline_distribution.png",
    "S4": "FigureS4_ahr2_residue_profile.png",
}

HEADING = re.compile(r"^### Figure (S?\d+)\.")
PROVENANCE = re.compile(r"^\*\*(Source|Panels):\*\*")


def transform(md: str) -> tuple[str, list[str]]:
    """Swap each figure's provenance lines for the composed image; unquote legends."""
    out: list[str] = []
    embedded: list[str] = []
    lines = md.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)

        m = HEADING.match(line)
        if not m:
            i += 1
            continue

        num = m.group(1)
        i += 1

        # drop the provenance block (a Source:/Panels: line plus its continuations)
        if i < len(lines) and PROVENANCE.match(lines[i]):
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(">"):
                i += 1

        path = FIGDIR / FIGURES[num]
        if not path.exists():
            raise SystemExit(f"missing composed figure for Figure {num}: {path}")
        out.append("")
        out.append(f"![]({path})")
        embedded.append(f"Figure {num} <- {path.name}")

        # legends are blockquoted in the Markdown so they read as captions; in Word they
        # should be ordinary paragraphs the editor can edit in place
        while i < len(lines) and (lines[i].startswith(">") or not lines[i].strip()):
            if lines[i].startswith(">"):
                out.append(re.sub(r"^> ?", "", lines[i]))
            else:
                out.append("")
                if i + 1 < len(lines) and not lines[i + 1].startswith(">"):
                    i += 1
                    break
            i += 1

    return "\n".join(out), embedded


def main() -> None:
    if not shutil.which("pandoc"):
        sys.exit("pandoc not found; install it (brew install pandoc)")

    md, embedded = transform(SRC.read_text())

    staged = Path(".manuscript_docx_src.md")
    staged.write_text(md)
    try:
        subprocess.run(
            ["pandoc", str(staged), "-o", str(OUT),
             "--from", "markdown+pipe_tables+tex_math_dollars",
             "--resource-path", ".",
             "--toc", "--toc-depth=2"],
            check=True,
        )
    finally:
        staged.unlink(missing_ok=True)

    for line in embedded:
        print(f"  {line}")
    print(f"\nWrote {OUT}  ({OUT.stat().st_size / 1e6:.1f} MB, {len(embedded)} figures embedded)")


if __name__ == "__main__":
    main()
