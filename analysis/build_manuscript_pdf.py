#!/usr/bin/env python3
"""
Render docs/MANUSCRIPT_FULL_en.md as a paginated PDF for co-authors and the English editor.

Markdown -> HTML (python-markdown, with tables) -> PDF (headless Chrome). Styled to read like a
journal manuscript: serif body, numbered running heads, tables that fit the page.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import markdown

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
SRC = Path("docs/MANUSCRIPT_FULL_en.md")
OUT = Path("docs/MANUSCRIPT_FULL_en.pdf")

CSS = """
@page { size: A4; margin: 20mm 18mm 18mm; }
* { box-sizing: border-box; }
body { font-family:"Hiragino Mincho ProN","Times New Roman",serif; font-size:10pt;
       line-height:1.65; color:#111; margin:0; }
h1 { font-size:16pt; line-height:1.4; text-align:center; margin:0 0 6mm;
     padding-bottom:4mm; border-bottom:1.5px solid #333; }
h2 { font-size:12.5pt; margin:8mm 0 3mm; padding-bottom:1mm;
     border-bottom:.6pt solid #999; page-break-after:avoid; }
h3 { font-size:11pt; margin:5mm 0 2mm; color:#1a3d6d; page-break-after:avoid; }
p { margin:0 0 2.6mm; text-align:justify; }
ul,ol { margin:0 0 3mm; padding-left:6mm; } li { margin-bottom:1.2mm; }
blockquote { margin:3mm 0; padding:2.5mm 4mm; background:#f4f6f9;
             border-left:3px solid #1a3d6d; font-size:9.4pt; }
blockquote p { margin:0 0 1.5mm; }
table { width:100%; border-collapse:collapse; margin:3mm 0 4mm; font-size:8.6pt;
        page-break-inside:auto; }
tr { page-break-inside:avoid; }
th,td { border:.4pt solid #999; padding:1.4mm 1.8mm; vertical-align:top; text-align:left; }
th { background:#e8edf4; font-weight:600; }
code { font-family:"SF Mono",Menlo,monospace; font-size:8.6pt; background:#f2f2f2;
       padding:0 2px; border-radius:2px; }
strong { font-weight:600; }
hr { border:0; border-top:.5pt solid #ccc; margin:6mm 0; }
a { color:#1a3d6d; word-break:break-all; }
"""


def main() -> None:
    if not SRC.is_file():
        sys.exit(f"missing {SRC}")
    html_body = markdown.markdown(
        SRC.read_text(encoding="utf-8"),
        extensions=["tables", "sane_lists", "attr_list"],
    )
    tmp = Path("/private/tmp/claude-501/manuscript_render.html")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    tmp.write_text(
        f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        f"<title>Manuscript</title><style>{CSS}</style></head><body>{html_body}</body></html>",
        encoding="utf-8",
    )
    subprocess.run(
        [CHROME, "--headless", "--disable-gpu", f"--print-to-pdf={OUT}",
         "--no-pdf-header-footer", str(tmp)],
        capture_output=True, check=False,
    )
    if not OUT.is_file():
        sys.exit("Chrome did not produce a PDF")
    print(f"Wrote {OUT}  ({OUT.stat().st_size/1048576:.1f} MB)")


if __name__ == "__main__":
    main()
