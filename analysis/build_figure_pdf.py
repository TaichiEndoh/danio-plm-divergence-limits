#!/usr/bin/env python3
"""
Bundle the final figures into (a) one downloadable PDF with legends, one figure per page,
and (b) a GitHub-rendered Markdown page that displays them inline.

Legends are pulled from docs/FIGURES_en.md so there is a single source of truth.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

FIN = Path("reports/figures/final")
DOCS = Path("docs")

ORDER = [
    ("Figure1_study_design.png", "Figure 1"),
    ("Figure2_scale_robustness.png", "Figure 2"),
    ("Figure3_prioritization.png", "Figure 3"),
    ("Figure4_hotspots_and_axes.png", "Figure 4"),
    ("Figure5_structural_grounding.png", "Figure 5"),
    ("Figure6_feature_annotation.png", "Figure 6"),
    ("FigureS1_rmsd_control.png", "Figure S1"),
    ("FigureS2_variant_rediscovery.png", "Figure S2"),
    ("FigureS3_baseline_distribution.png", "Figure S3"),
    ("FigureS4_ahr2_residue_profile.png", "Figure S4"),
]


def legends() -> dict[str, str]:
    """Extract each '> **Figure N.** …' block from docs/FIGURES_en.md."""
    text = (DOCS / "FIGURES_en.md").read_text(encoding="utf-8")
    out: dict[str, str] = {}
    for m in re.finditer(r"^> \*\*(Figure S?\d+)\.\*\*(.*?)(?=\n\n|\n#|\Z)", text, re.M | re.S):
        key = m.group(1)
        body = m.group(2)
        body = re.sub(r"^> ?", "", body, flags=re.M)          # strip quote markers
        body = re.sub(r"\*\*(.*?)\*\*", r"\1", body)          # strip bold
        body = re.sub(r"\*(.*?)\*", r"\1", body)              # strip italics
        body = re.sub(r"`(.*?)`", r"\1", body)
        out[key] = re.sub(r"\s+", " ", body).strip()
    return out


def wrap(text: str, width: int = 108) -> str:
    words, line, lines = text.split(), "", []
    for w in words:
        if len(line) + len(w) + 1 > width:
            lines.append(line)
            line = w
        else:
            line = f"{line} {w}".strip()
    if line:
        lines.append(line)
    return "\n".join(lines)


def main() -> None:
    leg = legends()
    out_pdf = FIN / "ALL_FIGURES.pdf"
    with PdfPages(out_pdf) as pdf:
        for fname, key in ORDER:
            path = FIN / fname
            if not path.is_file():
                print(f"  missing: {fname}")
                continue
            img = mpimg.imread(path)
            h, w = img.shape[0], img.shape[1]
            # Portrait page for tall figures, landscape for wide ones, so every figure fills
            # as much of its page as possible. The legend gets a proportionate strip below.
            tall = h / w > 1.0
            page = (8.3, 11.7) if tall else (11.7, 8.3)
            body = leg.get(key, "(legend not found)")
            width_chars = 96 if tall else 118
            wrapped = wrap(body, width_chars)
            n_lines = wrapped.count("\n") + 1
            # legend strip height as a fraction of the page
            legend_frac = min(0.34, 0.035 + n_lines * 0.0135 * (11.7 / page[1]))
            fig = plt.figure(figsize=page)
            ax = fig.add_axes([0.045, legend_frac + 0.035, 0.91, 0.925 - legend_frac])
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(key, fontsize=13, weight="bold", loc="left", pad=8)
            fig.text(0.045, legend_frac + 0.012, wrapped, fontsize=7.4, va="top", ha="left",
                     family="DejaVu Sans", linespacing=1.45)
            fig.text(0.96, 0.02, f"{fname}  ({w}×{h} px)", fontsize=6,
                     color="#666", ha="right")
            pdf.savefig(fig)
            plt.close(fig)
            print(f"  added {key}: {fname}")
        d = pdf.infodict()
        d["Title"] = "Figures — PLM ortholog divergence in Danio (follow-up study)"
        d["Subject"] = "All main and supplementary figures with legends"
    size = out_pdf.stat().st_size / 1048576
    print(f"\nWrote {out_pdf}  ({size:.1f} MB)")

    # GitHub-rendered preview page
    lines = [
        "# Figures (preview)",
        "",
        "All figures for the follow-up manuscript, displayed inline so they can be viewed on",
        "GitHub without downloading anything. Legends are the same text as in",
        "[`FIGURES_en.md`](../../docs/FIGURES_en.md).",
        "",
        "**Download everything as one file:** [`ALL_FIGURES.pdf`](ALL_FIGURES.pdf)",
        "(right-click → Save link as, or open it and use the download button).",
        "",
        "---",
        "",
    ]
    for fname, key in ORDER:
        if not (FIN / fname).is_file():
            continue
        lines += [f"## {key}", "", f"![{key}]({fname})", "",
                  f"> {leg.get(key, '(legend not found)')}", "", "---", ""]
    (FIN / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {FIN / 'README.md'}")


if __name__ == "__main__":
    main()
