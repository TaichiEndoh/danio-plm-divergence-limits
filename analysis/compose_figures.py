#!/usr/bin/env python3
"""
Compose the multi-panel main figures described in docs/FIGURES_en.md from the
single-panel plots in reports/figures/.

Panels are laid out on a white canvas, scaled to a common width per row, and labelled
(a), (b), (c)... in the top-left corner of each panel. Output is written to
reports/figures/final/ at publication resolution (>= 300 dpi at the intended print width).

Usage:  python analysis/compose_figures.py
"""
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

SRC = Path("reports/figures")
OUT = SRC / "final"
MARGIN = 26          # outer margin (px)
GAP = 22             # gap between panels (px)
LABEL_PAD = 10       # inset of the panel label from the panel corner
LABEL_PT = 58        # panel label font size
BG = "white"

FONT_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
    "/Library/Fonts/Arial Bold.ttf",
]


def load_font(size: int):
    for p in FONT_CANDIDATES:
        if Path(p).is_file():
            try:
                return ImageFont.truetype(p, size)
            except Exception:  # noqa: BLE001
                continue
    return ImageFont.load_default()


def scaled(img: Image.Image, width: int) -> Image.Image:
    if img.width == width:
        return img
    h = round(img.height * width / img.width)
    return img.resize((width, h), Image.LANCZOS)


def compose(name: str, rows: list[list[str]], target_w: int = 2400) -> Path:
    """rows: list of rows; each row is a list of source filenames (without .png)."""
    font = load_font(LABEL_PT)
    inner_w = target_w - 2 * MARGIN

    # scale each row to inner_w, splitting width between its panels
    laid: list[list[Image.Image]] = []
    for row in rows:
        n = len(row)
        each_w = (inner_w - GAP * (n - 1)) // n
        imgs = []
        for fn in row:
            p = SRC / f"{fn}.png"
            if not p.is_file():
                raise FileNotFoundError(p)
            im = Image.open(p).convert("RGB")
            imgs.append(scaled(im, each_w))
        laid.append(imgs)

    total_h = 2 * MARGIN + sum(max(i.height for i in r) for r in laid) + GAP * (len(laid) - 1)
    canvas = Image.new("RGB", (target_w, total_h), BG)
    draw = ImageDraw.Draw(canvas)

    label = iter("abcdefghij")
    y = MARGIN
    for r in laid:
        x = MARGIN
        rh = max(i.height for i in r)
        for im in r:
            canvas.paste(im, (x, y))
            ch = next(label)
            # white halo behind the label so it stays readable over plot content
            for dx, dy in ((-2, 0), (2, 0), (0, -2), (0, 2)):
                draw.text((x + LABEL_PAD + dx, y + LABEL_PAD + dy), f"({ch})",
                          font=font, fill="white")
            draw.text((x + LABEL_PAD, y + LABEL_PAD), f"({ch})", font=font, fill="black")
            x += im.width + GAP
        y += rh + GAP

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / f"{name}.png"
    # 300 dpi metadata for a ~8 inch print width
    canvas.save(dest, dpi=(300, 300))
    print(f"  {dest.name}: {canvas.width}x{canvas.height}px, {len([i for r in laid for i in r])} panels")
    return dest


def copy_single(name: str, src: str) -> None:
    im = Image.open(SRC / f"{src}.png").convert("RGB")
    OUT.mkdir(parents=True, exist_ok=True)
    im.save(OUT / f"{name}.png", dpi=(300, 300))
    print(f"  {name}.png: {im.width}x{im.height}px, 1 panel (as-is)")


def main() -> None:
    print("Composing main figures ->", OUT)

    # Fig 1 — concept diagram, single panel, no label needed
    copy_single("Figure1_study_design", "fig1_concept")

    # Fig 2 — ranking preserved across scale
    compose("Figure2_scale_robustness",
            [["fig2_esm2_vs_esm3_subset", "esm2_ladder_to_esm3"]])

    # Fig 3 — candidate prioritization vs chance + identity control
    compose("Figure3_prioritization",
            [["fig3_ranking_overlap_esm3"],
             ["fig_ranking_stats"]])

    # Fig 4 — residue hotspots + ortholog/paralog axis
    compose("Figure4_hotspots_and_axes",
            [["esm2_vs_esm3_divergence"],
             ["ahr_pocket_ortholog_vs_paralog", "ahr_pocket_axis_esm2_esm3"]])

    # Fig 5 — structural grounding
    compose("Figure5_structural_grounding",
            [["ahr2_plddt_vs_divergence", "kcnj13_plddt_vs_divergence"],
             ["ahr2_structure_divergence", "kcnj13_structure_divergence"],
             ["tier_stratification"]])

    # Fig 6 — feature annotation (already two panels in one image)
    copy_single("Figure6_feature_annotation", "feature_annotation")

    print("\nComposing supplementary figures")
    compose("FigureS1_rmsd_control",
            [["ahr2_plm_vs_structure", "kcnj13_plm_vs_structure"]])
    copy_single("FigureS2_variant_rediscovery", "fig1_variant_rediscovery_scaling")
    copy_single("FigureS3_baseline_distribution", "fig2_baseline_divergence_distribution")
    copy_single("FigureS4_ahr2_residue_profile", "fig3_ahr2_residue_divergence")

    print(f"\nDone. {len(list(OUT.glob('*.png')))} files in {OUT}")


if __name__ == "__main__":
    main()
