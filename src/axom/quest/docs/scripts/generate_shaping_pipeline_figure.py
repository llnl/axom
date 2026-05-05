#!/usr/bin/env python3
"""Generate an SVG figure for Quest's shaping pipeline documentation."""

from __future__ import annotations

import argparse
import html
from pathlib import Path


WIDTH = 1600
HEIGHT = 900


def esc(text: str) -> str:
    return html.escape(text, quote=True)


def rect(x: int, y: int, w: int, h: int, cls: str, rx: int = 28) -> str:
    return (
        f"<rect class='{cls}' x='{x}' y='{y}' width='{w}' height='{h}' "
        f"rx='{rx}' ry='{rx}' />"
    )


def line(x1: int, y1: int, x2: int, y2: int, cls: str, marker: bool = True) -> str:
    marker_end = " marker-end='url(#arrow)'" if marker else ""
    return (
        f"<line class='{cls}' x1='{x1}' y1='{y1}' x2='{x2}' y2='{y2}'{marker_end} />"
    )


def text_block(x: int, y: int, lines: list[str], cls: str, line_gap: int = 34) -> str:
    parts = [f"<text class='{cls}' x='{x}' y='{y}'>"]
    for idx, item in enumerate(lines):
        dy = 0 if idx == 0 else line_gap
        parts.append(f"<tspan x='{x}' dy='{dy}'>{esc(item)}</tspan>")
    parts.append("</text>")
    return "".join(parts)


def bullet_list(x: int, y: int, items: list[str], cls: str, line_gap: int = 30) -> str:
    parts: list[str] = []
    for idx, item in enumerate(items):
        yy = y + idx * line_gap
        parts.append(f"<circle class='bullet-dot' cx='{x}' cy='{yy - 6}' r='4' />")
        parts.append(f"<text class='{cls}' x='{x + 18}' y='{yy}'>{esc(item)}</text>")
    return "".join(parts)


def stage_box(x: int, y: int, w: int, h: int, label: str, items: list[str], kind: str) -> str:
    return "".join(
        [
            rect(x, y, w, h, f"panel {kind}"),
            text_block(x + 28, y + 52, [label], "stage-title"),
            bullet_list(x + 32, y + 110, items, "body-text"),
        ]
    )


def step_box(
    x: int, y: int, w: int, h: int, title: list[str], subtitle: list[str]
) -> str:
    return "".join(
        [
            rect(x, y, w, h, "step"),
            text_block(x + 20, y + 36, title, "step-title", line_gap=22),
            text_block(x + 20, y + 74, subtitle, "step-text", line_gap=22),
        ]
    )


def dashed_loop(x: int, y: int, w: int, h: int) -> str:
    return (
        f"<rect class='loop-frame' x='{x}' y='{y}' width='{w}' height='{h}' "
        f"rx='34' ry='34' />"
    )


def build_svg() -> str:
    out: list[str] = []
    out.append(
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{WIDTH}' height='{HEIGHT}' "
        f"viewBox='0 0 {WIDTH} {HEIGHT}' role='img' "
        f"aria-labelledby='title desc'>"
    )
    out.append("<title id='title'>Quest shaping pipeline</title>")
    out.append(
        "<desc id='desc'>Diagram showing Klee shape input flowing through Quest "
        "shaper setup, a per-shape processing loop, and output material volume "
        "fraction fields on the target mesh.</desc>"
    )
    out.append(
        """
<defs>
  <linearGradient id="bg" x1="0" y1="0" x2="1" y2="1">
    <stop offset="0%" stop-color="#f7f9fb" />
    <stop offset="100%" stop-color="#eef3f7" />
  </linearGradient>
  <linearGradient id="meshGlow" x1="0" y1="0" x2="1" y2="1">
    <stop offset="0%" stop-color="#dbe8f4" stop-opacity="0.85" />
    <stop offset="100%" stop-color="#f7fafc" stop-opacity="0.35" />
  </linearGradient>
  <filter id="shadow" x="-20%" y="-20%" width="140%" height="140%">
    <feDropShadow dx="0" dy="10" stdDeviation="18" flood-color="#2d3742" flood-opacity="0.12" />
  </filter>
  <marker id="arrow" markerWidth="8" markerHeight="8" refX="6" refY="4" orient="auto">
    <path d="M0,0 L8,4 L0,8 z" fill="#48667d" />
  </marker>
</defs>
<style>
  .bg { fill: url(#bg); }
  .halo-a { fill: #d9e7f3; opacity: 0.55; }
  .halo-b { fill: #e4e8ec; opacity: 0.72; }
  .panel { filter: url(#shadow); stroke-width: 2.5; }
  .input { fill: #f4f7fa; stroke: #5f7385; }
  .setup { fill: #e9f2f9; stroke: #2980b9; }
  .output { fill: #eef3f7; stroke: #607d94; }
  .step { fill: #ffffff; stroke: #70818f; stroke-width: 2; }
  .loop-frame { fill: none; stroke: #6b8193; stroke-width: 3; stroke-dasharray: 12 12; }
  .connector { stroke: #48667d; stroke-width: 3; stroke-linecap: round; }
  .connector-soft { stroke: #7a90a3; stroke-width: 2.5; stroke-linecap: round; stroke-dasharray: 8 8; }
  .stage-title { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 30px; font-weight: 700; fill: #343131; }
  .loop-title { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 28px; font-weight: 700; fill: #3b4d5d; }
  .step-title { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 20px; font-weight: 700; fill: #344654; }
  .step-text { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 15px; fill: #556977; }
  .body-text { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 18px; fill: #445867; }
  .small-label { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 18px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; fill: #627583; }
  .caption { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 20px; fill: #586874; }
  .mono { font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', monospace; font-size: 20px; fill: #344654; }
  .bullet-dot { fill: #2980b9; }
</style>
"""
    )
    out.append(f"<rect class='bg' x='0' y='0' width='{WIDTH}' height='{HEIGHT}' />")
    out.append("<circle class='halo-a' cx='1280' cy='150' r='180' />")
    out.append("<circle class='halo-b' cx='240' cy='760' r='170' />")

    out.append(text_block(72, 88, ["Quest shaping pipeline"], "stage-title"))
    out.append(
        text_block(
            74,
            130,
            ["Klee describes shapes. Quest loads them, queries the target mesh,", "and writes material volume fractions."],
            "caption",
            line_gap=28,
        )
    )

    out.append(text_block(100, 214, ["Klee input"], "small-label"))
    out.append(
        stage_box(
            72,
            234,
            340,
            250,
            "Shape set",
            [
                "material name",
                "geometry path and format",
                "transforms",
                "replacement rules",
            ],
            "input",
        )
    )

    out.append(text_block(490, 214, ["Quest setup"], "small-label"))
    out.append(
        stage_box(
            462,
            234,
            340,
            250,
            "Target mesh + shaper",
            [
                "MFEM or Blueprint mesh",
                "sampling or intersection shaper",
                "runtime policy and tolerances",
                "optional initial vol_frac_*",
            ],
            "setup",
        )
    )

    loop_x = 842
    loop_y = 180
    loop_w = 500
    loop_h = 540
    out.append(text_block(loop_x + 32, 214, ["Per-shape loop"], "small-label"))
    out.append(dashed_loop(loop_x, loop_y + 38, loop_w, loop_h))
    out.append(text_block(loop_x + 34, loop_y + 98, ["For each shape"], "loop-title"))

    step_x = loop_x + 30
    step_w = 212
    step_h = 100
    row1_y = loop_y + 146
    row2_y = loop_y + 314

    col2_x = step_x + 238
    bottom_x = step_x + 120

    out.append(
        step_box(
            step_x,
            row1_y,
            step_w,
            step_h,
            ["loadShape()"],
            ["read geometry", "and apply transforms"],
        )
    )
    out.append(
        step_box(
            col2_x,
            row1_y,
            step_w,
            step_h,
            ["prepare", "ShapeQuery()"],
            ["build spatial index", "or intersection data"],
        )
    )
    out.append(
        step_box(
            step_x,
            row2_y,
            step_w,
            step_h,
            ["runShapeQuery()"],
            ["measure overlap", "with target zones"],
        )
    )
    out.append(
        step_box(
            col2_x,
            row2_y,
            step_w,
            step_h,
            ["applyReplacement", "Rules()"],
            ["merge with existing", "materials"],
        )
    )
    out.append(
        step_box(
            bottom_x,
            loop_y + 442,
            step_w + 52,
            step_h,
            ["finalizeShapeQuery()"],
            ["clear temporary state", "before the next shape"],
        )
    )

    out.append(line(412, 360, 462, 360, "connector"))
    out.append(line(802, 360, 842, 360, "connector"))
    out.append(line(step_x + step_w, row1_y + 50, col2_x, row1_y + 50, "connector"))
    out.append(line(step_x + 106, row1_y + step_h, step_x + 106, row2_y, "connector"))
    out.append(line(col2_x + 106, row1_y + step_h, col2_x + 106, row2_y, "connector"))
    out.append(line(step_x + step_w, row2_y + 50, col2_x, row2_y + 50, "connector"))
    out.append(line(step_x + 106, row2_y + step_h, bottom_x + 110, loop_y + 442, "connector"))
    out.append(line(col2_x + 106, row2_y + step_h, bottom_x + 154, loop_y + 442, "connector"))
    out.append(line(bottom_x + 132, loop_y + 542, bottom_x + 132, loop_y + 584, "connector-soft"))
    out.append(text_block(bottom_x + 154, loop_y + 578, ["next shape"], "step-text"))

    out.append(text_block(1406, 214, ["Quest output"], "small-label"))
    out.append(
        "".join(
            [
                rect(1378, 234, 180, 250, "panel output"),
                text_block(1406, 286, ["Volume", "fractions"], "stage-title", line_gap=30),
                bullet_list(
                    1410,
                    346,
                    [
                        "vol_frac_matA",
                        "vol_frac_matB",
                        "vol_frac_free",
                        "ready for output",
                    ],
                    "body-text",
                ),
            ]
        )
    )
    out.append(line(1342, 360, 1378, 360, "connector"))

    out.append(rect(1354, 552, 240, 140, "setup", rx=24))
    out.append(text_block(1380, 598, ["Target mesh state"], "step-title"))
    out.append(text_block(1380, 638, ["materials now encode", "shape overlap per zone"], "step-text", line_gap=26))
    out.append(line(1468, 484, 1468, 552, "connector"))

    out.append("<path d='M 632 486 C 612 540, 662 612, 842 612' class='connector-soft' fill='none' />")
    out.append(line(632, 484, 632, 486, "connector-soft", marker=False))
    out.append(text_block(500, 592, ["optional: import existing", "material fields before shaping"], "step-text", line_gap=26))

    out.append("</svg>")
    return "".join(out)


def parse_args() -> argparse.Namespace:
    default_output = (
        Path(__file__).resolve().parent.parent / "sphinx" / "figs" / "shaping_pipeline.svg"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=default_output,
        help=f"Path to the output SVG file (default: {default_output})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(build_svg(), encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
