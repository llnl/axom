#!/bin/sh
"exec" "python3" "-u" "-B" "$0" "$@"

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""
 file: svg2contours.py

 description: 
  Reads in an SVG document and outputs an MFEM NURBS mesh.
  Depends on the svgpathtools module

 notes:
  - By default, all SVG curve segments are output as cubic NURBS (degree 3) for
    compatibility with older MFEM/VisIt workflows.
  - This script optionally supports MFEM's newer "patches" NURBS mesh format
    for 1D NURBS segments embedded in 2D. Support for reading patch-based 1D
    NURBS meshes was added to MFEM after MFEM 4.9.0.
"""

import sys
import os
import json
from typing import Optional
from svgpathtools import (
    Document,
    Path,
    Line,
    QuadraticBezier,
    CubicBezier,
    Arc,
    is_bezier_path,
    is_bezier_segment,
    is_path_segment,
    svg2paths,
    bpoints2bezier,
)
import numpy as np
import re
import argparse


def parse_viewbox(view_box: str):
    """Parse an SVG viewBox attribute string.

    Returns a 4-tuple (min_x, min_y, width, height) as floats, or None if
    view_box is not provided or is invalid.
    """

    if not view_box:
        return None

    # SVG spec allows whitespace and/or comma separators.
    parts = [p for p in re.split(r"[,\s]+", view_box.strip()) if p]
    if len(parts) != 4:
        return None

    try:
        return tuple(map(float, parts))
    except Exception:
        return None


def parse_svg_length(length: str):
    """Parse an SVG length attribute and return its numeric value as a float.

    If the value is missing or expressed as a percentage, returns None.
    Units (e.g. 'mm', 'px') are ignored and the leading numeric portion is used.
    """

    if not length:
        return None

    length = length.strip()
    if length.endswith("%"):
        return None

    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", length)
    if not match:
        return None

    try:
        return float(match.group(0))
    except Exception:
        return None


def get_root_transform(doc: Document):
    """Create transform to convert 2D coordinate system from y pointing down to y pointing up"""

    attr = doc.tree.getroot().attrib
    width = attr.get("width", None)
    height = attr.get("height", None)
    viewBox = attr.get("viewBox", None)

    print(f"""SVG dimensions: {width=} {height=} {viewBox=}""")

    tf = np.identity(3)

    # Prefer viewBox since it defines the internal user-coordinate system used by path data.
    parsed_viewbox = parse_viewbox(viewBox)
    if parsed_viewbox is not None:
        _, min_y, _, vb_height = parsed_viewbox
        tf[1, 1] = -1
        tf[1, 2] = 2.0 * min_y + vb_height
        return tf

    # Fallback to root height attribute when viewBox is unavailable/invalid.
    height_number = parse_svg_length(height)
    if height_number is not None:
        tf[1, 1] = -1
        tf[1, 2] = height_number

    # Last-resort fallback: flip about y=0 when we cannot infer a reference height.
    if tf[1, 1] == 1:
        tf[1, 1] = -1

    return tf


def transform_cubic(cubic: CubicBezier, tf):
    """Apply transformation `tf` to control points of cubic Bezier.
    Adapted from internal `transform()` function in svgpathtools library
    """

    def to_point(p):
        return np.array([[p.real], [p.imag], [1.0]])

    def to_complex(v):
        return v.item(0) + 1j * v.item(1)

    return bpoints2bezier([to_complex(tf.dot(to_point(p))) for p in cubic.bpoints()])


def transform_segment(seg, tf):
    """Apply transformation `tf` to a Line/QuadraticBezier/CubicBezier segment."""

    def to_point(p):
        return np.array([[p.real], [p.imag], [1.0]])

    def to_complex(v):
        return v.item(0) + 1j * v.item(1)

    return bpoints2bezier([to_complex(tf.dot(to_point(p))) for p in seg.bpoints()])


def lerp(a, b, t):
    """linear interpolation from a to b with parameter t, typically between 0 and 1"""
    return (1 - t) * a + t * b


def line_to_cubic(line: Line):
    """Converts from an svgpathtools Line to a (rational) CubicBezier (with unit weights)"""
    q_0, q_1 = line.bpoints()
    return (CubicBezier(q_0, lerp(q_0, q_1, 1.0 / 3), lerp(q_0, q_1, 2.0 / 3), q_1), [1, 1, 1, 1])


def quadratic_to_cubic(quad: QuadraticBezier):
    """Converts from an svgpathtools QuadraticBezier to a (rational) CubicBezier (with unit weights)"""
    q_0, q_1, q_2 = quad.bpoints()
    return (CubicBezier(q_0, lerp(q_0, q_1, 2.0 / 3), lerp(q_1, q_2, 1.0 / 3), q_2), [1, 1, 1, 1])


def arc_to_cubic(arc: Arc):
    """Convertes from an svgpathtools Arc to a rational CubicBezier"""
    q_0 = arc.start
    q_3 = arc.end

    def area(p1, p2, p3):
        """Computes the area of a triangle defined by vertices p1, p2 and p3"""
        v21, v31 = p2 - p1, p3 - p1
        return 0.5 * np.abs(v21.real * v31.imag - v21.imag * v31.real)

    # Notes:
    # (1) We have control point positions 0 and 3 and their tangents on the ellipse
    # as well as the midpoint; we need to find control points 1 and 2
    # (2) We are computing these as the intersections of the tangent lines, with lines connecting
    # the opposite endpoint to the midpoints
    # (3) The weights are then derived through an isometry with the semicirle,
    # whose weights are proportional to [3,1,1,3]
    try:
        # Use a scaling factor to extend lines from endpoints to shoulder
        # these lines contain the internal control points c_1 and c_2
        scale_fac = 10

        shoulder = arc.point(0.5)
        d_0 = arc.derivative(0)
        d_1 = arc.derivative(1)
        # print(f"""For arc {arc} w/ {d_0=} and {d_1=};\n {arc.theta=}, {arc.phi=}, {arc.rotation=}, {arc.delta=}""")

        # extend the line segment from q3 to shoulder point
        # and find the intersection w/ tangent line @ control point 0
        l_3_1 = Line(q_3, q_3 + scale_fac * (shoulder - q_3))
        l_0_1 = Line(q_0, q_0 + d_0)
        ints_31_01 = l_3_1.intersect(l_0_1)
        # print(f"""Finding intersection @ control point 1\n  {shoulder=}\n  {l_3_1=}\n  {l_0_1=}\n  {ints_31_01=}""")

        # extend the line segment from q0 to shoulder point
        # and find the intersection w/ tangent line @ control point 3
        l_0_2 = Line(q_0, q_0 + scale_fac * (shoulder - q_0))
        l_3_2 = Line(q_3, q_3 - d_1)
        ints_02_32 = l_0_2.intersect(l_3_2)
        # print(f"""Finding intersection @ control point 2\n  {shoulder=}\n  {l_0_2=}\n  {l_3_2=}\n  {ints_02_32=}""")

        c_1 = l_0_1.point(ints_31_01[0][1])
        c_2 = l_3_2.point(ints_02_32[0][1])
        # print(f"""Control points from intersections: CP 1: {c_1}; CP 2 {c_2}""")

        # print(f"""\t<path d="M {q_3.real} {q_3.imag} {c_2.real} {c_2.imag} {c_1.real} {c_1.imag} {q_0.real} {q_0.imag}" />""")
        # print(f"""\t<circle cx="{shoulder.real}" cy="{shoulder.imag}" r="5" />""")

        reversed = not arc.sweep
        curve = CubicBezier(q_3, c_2, c_1, q_0) if reversed else CubicBezier(q_0, c_1, c_2, q_3)

        # compute the rational weights based on areas of triangles that skip the current index
        # formula is from "Shape factors and shoulder points for shape control of rational Bezier curves"
        #                  https://doi.org/10.1016/j.cad.2023.103477
        b0, b1, b2, b3 = curve.bpoints()

        areas = [area(b1, b2, b3), area(b0, b2, b3) / 3.0, area(b0, b1, b3) / 3.0, area(b1, b2, b3)]
        shape_fac = [areas[1] ** 2 / (areas[0] * areas[2]), areas[2] ** 2 / (areas[1] * areas[3])]
        weights = [1, 0, 0, 1]
        weights[1] = np.cbrt(shape_fac[0] ** 2 * shape_fac[1])
        weights[2] = weights[1] ** 2 / shape_fac[0]
        # print(f""" -- {weights=} and {shape_fac=} for {curve=}""")

        return (curve, weights)

    except Exception as err:
        print(f"Exception: {err}")
        print(f"""*** Problem with arc {arc}:\n\t {arc.theta=}; {arc.delta=}; {arc.phi=} ***""")
        # as a fall-back, use as_cubic_curves function from svgpathtools
        # which approximates rational curve
        # note, we're currently only taking the first cubic; there might be more.
        for c in arc.as_cubic_curves():
            q_0, q_1 = c.start, c.end
            c_1, c_2 = c.control1, c.control2
            return (CubicBezier(q_0, c_1, c_2, q_1), [3, 1, 1, 3])


def arc_to_quadratic_beziers(arc: Arc, *, max_sweep_angle_rad: float = np.pi / 2):
    """Convert an svgpathtools Arc to rational QuadraticBezier segments.

    Notes:
      - Rational quadratics represent conic sections exactly.
      - We subdivide the arc into pieces (default: <= 90 degrees) to keep the
        interior weight bounded away from zero.
      - This function avoids relying on svgpathtools' internal angle units for
        `arc.theta` / `arc.delta` (which may be degrees depending on version).
        Instead, it reconstructs the start angle and sweep angle from the arc's
        start/end points in the ellipse's local parameter space, then selects
        the correct branch using `arc.point(0.5)` (and `arc.large_arc` as a hint).
    Returns:
      List[(QuadraticBezier, weights)] where weights has length 3 and the degree is 2.
    """

    if max_sweep_angle_rad <= 0:
        raise ValueError("max_sweep_angle_rad must be positive")

    center = getattr(arc, "center", None)
    radius = getattr(arc, "radius", None)

    if center is None or radius is None:
        raise Exception("Arc is missing required attributes (center/radius)")

    rx = float(np.abs(radius.real))
    ry = float(np.abs(radius.imag))

    if rx == 0.0 or ry == 0.0:
        # Degenerate arc; let caller fall back (svgpathtools generally models these as Lines).
        return []

    # We build rational quadratic pieces by using the conic (tangent intersection) construction.
    # This avoids any ambiguity about svgpathtools' internal angle units and tends to be robust.

    def cross(a: complex, b: complex) -> float:
        return float(a.real * b.imag - a.imag * b.real)

    def dot(a: complex, b: complex) -> float:
        return float(a.real * b.real + a.imag * b.imag)

    def try_quad_for_arc_piece(arc_piece: Arc):
        p0 = arc_piece.start
        p2 = arc_piece.end
        m = arc_piece.point(0.5)

        d0 = arc_piece.derivative(0.0)
        d2 = arc_piece.derivative(1.0)

        # Compute intersection of tangents at endpoints.
        denom = cross(d0, d2)
        if np.abs(denom) < 1e-14:
            return None

        # p0 + s*d0 = p2 + t*d2
        s = cross((p2 - p0), d2) / denom
        p1 = p0 + s * d0

        # Solve for the middle weight w s.t. the rational quadratic matches the midpoint:
        # M = (P0 + 2 w P1 + P2) / (2(1+w))  =>  2 w (M - P1) = P0 + P2 - 2 M
        u = m - p1
        uu = dot(u, u)
        if uu <= 0.0:
            return None

        rhs = p0 + p2 - 2.0 * m
        w = dot(rhs, u) / (2.0 * uu)

        if not np.isfinite(w) or w <= 0.0:
            return None

        # Consistency check (least-squares residual; should be near zero for true conics).
        res = (p0 + p2 + 2.0 * w * p1) - (2.0 * (1.0 + w) * m)
        scale2 = max(rx, ry) ** 2
        if dot(res, res) > 1e-8 * max(1.0, scale2):
            return None

        return (QuadraticBezier(p0, p1, p2), [1.0, float(w), 1.0])

    # Split until each piece's w is at least cos(max_sweep/2) (matches the circle-arc weight bound).
    w_min = float(np.cos(0.5 * max_sweep_angle_rad))
    max_pieces = 1024

    pieces = []
    work = [arc]
    while work:
        if len(pieces) + len(work) > max_pieces:
            return []

        a = work.pop(0)
        quad = try_quad_for_arc_piece(a)
        if quad is None:
            a1, a2 = a.split(0.5)
            work.insert(0, a2)
            work.insert(0, a1)
            continue

        _, wts = quad
        if wts[1] < w_min - 1e-12:
            a1, a2 = a.split(0.5)
            work.insert(0, a2)
            work.insert(0, a1)
            continue

        pieces.append(quad)

    # Match legacy `arc_to_cubic` orientation handling: reverse when sweep flag is not set.
    if not getattr(arc, "sweep", True):
        pieces = [(q.reversed(), list(reversed(w))) for (q, w) in reversed(pieces)]

    return pieces


def segment_as_cubic(seg, reverse_paths: bool):
    if isinstance(seg, Line):
        cubic, weights = line_to_cubic(seg)
    elif isinstance(seg, QuadraticBezier):
        cubic, weights = quadratic_to_cubic(seg)
    elif isinstance(seg, CubicBezier):
        cubic, weights = seg, [1, 1, 1, 1]
    elif isinstance(seg, Arc):
        cubic, weights = arc_to_cubic(seg)
    else:
        raise Exception(f"'{type(seg)}' type not supported yet")

    if reverse_paths:
        cubic = cubic.reversed()
        weights.reverse()

    return (cubic, weights)


def segment_as_native_nurbs_segments(seg, reverse_paths: bool):
    """Convert an svgpathtools segment to Bezier (NURBS) segments with native degree where possible.

    Returns a list of (segment, weights, degree).
      - Line: 1 segment, degree 1
      - QuadraticBezier: 1 segment, degree 2
      - CubicBezier: 1 segment, degree 3
      - Arc: 1+ rational QuadraticBezier segments (degree 2), subdivided for robustness
    """

    out = []
    if isinstance(seg, Line):
        out = [(seg, [1, 1], 1)]
    elif isinstance(seg, QuadraticBezier):
        out = [(seg, [1, 1, 1], 2)]
    elif isinstance(seg, CubicBezier):
        out = [(seg, [1, 1, 1, 1], 3)]
    elif isinstance(seg, Arc):
        quads = arc_to_quadratic_beziers(seg)
        if not quads:
            # Fall back to a rational cubic if the arc conversion is ill-conditioned.
            cubic, weights = arc_to_cubic(seg)
            out = [(cubic, weights, 3)]
        else:
            out = [(q, w, 2) for (q, w) in quads]
    else:
        raise Exception(f"'{type(seg)}' type not supported yet")

    if reverse_paths:
        out = [(s.reversed(), list(reversed(w)), d) for (s, w, d) in reversed(out)]

    return out


def dist_to_ellipse(center, radius, angle, pt):
    cx, cy = center.real, center.imag
    rx, ry = radius.real, radius.imag

    rot = np.exp(-1j * np.radians(angle))
    transformed_pt = rot * complex(pt.real - cx, pt.imag - cy)
    return transformed_pt.real ** 2 / rx ** 2 + transformed_pt.imag ** 2 / ry ** 2 - 1


class MFEMData:

    def __init__(self):
        self.elem_cnt = 0
        self.vert_cnt = 0
        self.elems = []
        self.edges = []
        self.knots = []

        # mfem format lists the endpoints and then the interiors
        self.wgts_ends = []
        self.wgts_ints = []
        self.dof_ends = []
        self.dof_ints = []

    def add_cubic_bezier(self, cubic, weights, attrib):

        self.elems.append(" ".join(map(str, [attrib, 1, self.vert_cnt, self.vert_cnt + 1])))
        self.vert_cnt += 2

        self.edges.append(f"{self.elem_cnt} 0 1")
        self.elem_cnt += 1

        # Assume for now that the order is always 3
        self.knots.append("3 4 0 0 0 0 1 1 1 1")
        self.wgts_ends.append(f"{weights[0]} {weights[3]}")
        self.wgts_ints.append(f"{weights[2]} {weights[1]}")
        self.dof_ends.append(" ".join(
            map(str, [cubic.start.real, cubic.start.imag, cubic.end.real, cubic.end.imag])))
        self.dof_ints.append(" ".join(
            map(str, [
                cubic.control2.real, cubic.control2.imag, cubic.control1.real, cubic.control1.imag
            ])))

    def write_file(self, filename):
        mfem_file = []

        mfem_file.extend([
            "MFEM NURBS mesh v1.0",
            "",
            "# MFEM Geometry Types (see fem/geom.hpp):",
            "#",
            "# SEGMENT = 1 | SQUARE = 3 | CUBE = 5",
            "#",
            "# element: <attr> 1 <v0> <v1>",
            "# edge: <idx++> 0 1  <-- idx increases by one each time",
            "# knotvector: <order> <num_ctrl_pts> [knots]; sizeof(knots) is 1+order+num_ctrl_pts",
            "# weights: array of weights corresponding to the NURBS element",
            "# FES: list of control points; vertex control points at top, then interior control points",
            "",
        ])

        mfem_file.extend(["dimension", "1", ""])

        mfem_file.extend(["elements", f"{self.elem_cnt}", "\n".join(self.elems), ""])

        mfem_file.extend(["boundary", "0", ""])

        mfem_file.extend(["edges", f"{self.elem_cnt}", "\n".join(self.edges), ""])

        mfem_file.extend(["vertices", f"{self.vert_cnt}", ""])

        mfem_file.extend(["knotvectors", f"{self.elem_cnt}", "\n".join(self.knots), ""])

        mfem_file.extend(["weights", "\n".join(self.wgts_ends), "\n".join(self.wgts_ints), ""])

        mfem_file.extend([
            "FiniteElementSpace",
            "FiniteElementCollection: NURBS",
            "VDim: 2",
            "Ordering: 1",
            "",
            "\n".join(self.dof_ends),
            "\n".join(self.dof_ints),
            "",
        ])

        with open(filename, mode="w") as f:
            f.write("\n".join(mfem_file))


class MFEMPatchesData:

    def __init__(self):
        self.elem_cnt = 0
        self.vert_cnt = 0
        self.elems = []
        self.edges = []
        self.patches = []

    @staticmethod
    def _bezier_knotvector(degree: int):
        # Bezier knot vector on [0,1]
        return [0] * (degree + 1) + [1] * (degree + 1)

    @staticmethod
    def _quadratic_multispan_knotvector(num_spans: int, *, degree: int = 2):
        if degree != 2:
            raise ValueError("Only quadratic multispan knotvectors are supported here")
        if num_spans < 1:
            raise ValueError("num_spans must be >= 1")

        knots = [0.0] * (degree + 1)
        if num_spans > 1:
            # Use multiplicity=degree at internal knots so each span is a Bezier segment.
            for i in range(1, num_spans):
                t = float(i) / float(num_spans)
                knots.extend([t] * degree)
        knots.extend([1.0] * (degree + 1))
        return knots

    @staticmethod
    def quadratic_beziers_to_multispan(quads_with_weights):
        """Merge quadratic rational Bezier spans into a single quadratic multi-span NURBS patch.

        This uses internal knot multiplicity=degree (2), so each span remains a Bezier segment.
        Returns (cps, weights, knots).
        """

        num_spans = len(quads_with_weights)
        if num_spans < 1:
            raise ValueError("Expected at least one span")

        cps = []
        weights = []
        for span_idx, (quad, wts) in enumerate(quads_with_weights):
            bpts = quad.bpoints()
            if len(bpts) != 3 or len(wts) != 3:
                raise Exception(
                    "Expected quadratic Bezier spans with 3 control points and 3 weights")
            if span_idx == 0:
                cps.extend(bpts)
                weights.extend(wts)
            else:
                cps.extend(bpts[1:])
                weights.extend(wts[1:])

        knots = MFEMPatchesData._quadratic_multispan_knotvector(num_spans)
        return cps, weights, knots

    def add_nurbs_patch(self,
                        *,
                        cps,
                        degree: int,
                        weights,
                        knots,
                        attrib: int,
                        patch_comment: Optional[str] = None):
        if len(cps) != len(weights):
            raise Exception(f"Expected {len(cps)} weights, got {len(weights)}")
        expected_knots = 1 + degree + len(cps)
        if len(knots) != expected_knots:
            raise Exception(f"Expected {expected_knots} knots, got {len(knots)}")

        v0 = self.vert_cnt
        v1 = self.vert_cnt + 1
        self.vert_cnt += 2

        self.elems.append(" ".join(map(str, [attrib, 1, v0, v1])))
        self.edges.append(f"{self.elem_cnt} {v0} {v1}")
        self.elem_cnt += 1

        patch_lines = []
        patch_lines.append("")
        if patch_comment:
            patch_lines.append(f"# Patch {self.elem_cnt - 1}: {patch_comment}")
        else:
            patch_lines.append(
                f"# Patch {self.elem_cnt - 1}: degree {degree} ({len(cps)} control points)")
        patch_lines.append("knotvectors")
        patch_lines.append("1")
        patch_lines.append("{} {}  {}".format(
            degree,
            len(cps),
            " ".join(str(k) for k in knots),
        ))
        patch_lines.append("")
        patch_lines.append("dimension")
        patch_lines.append("2")
        patch_lines.append("")
        patch_lines.append("controlpoints")
        for (cp, w) in zip(cps, weights):
            # MFEM expects NURBS control points in homogeneous form: (x*w, y*w, w).
            ww = float(w)
            patch_lines.append(f"{cp.real * ww}  {cp.imag * ww}  {ww}")
        patch_lines.append("")

        self.patches.append("\n".join(patch_lines))

    def add_bezier(self, seg, degree: int, weights, attrib: int):
        cps = seg.bpoints()
        if len(cps) != degree + 1:
            raise Exception(
                f"Expected {degree + 1} control points for degree {degree}, got {len(cps)}")
        self.add_nurbs_patch(
            cps=cps,
            degree=degree,
            weights=weights,
            knots=self._bezier_knotvector(degree),
            attrib=attrib,
        )

    def write_file(self, filename):
        mfem_file = []

        mfem_file.extend([
            "MFEM NURBS mesh v1.0",
            "",
            "#",
            "# Patch-based 1D NURBS segments embedded in 2D.",
            "# NOTE: MFEM support for reading patch-based 1D NURBS meshes was added after MFEM 4.9.0.",
            "#",
            "",
        ])

        mfem_file.extend(["dimension", "1", ""])
        mfem_file.extend(["elements", f"{self.elem_cnt}", "\n".join(self.elems), ""])
        mfem_file.extend(["boundary", "0", ""])
        mfem_file.extend(["edges", f"{self.elem_cnt}", "\n".join(self.edges), ""])
        mfem_file.extend(["vertices", f"{self.vert_cnt}", ""])
        mfem_file.extend(["patches", "\n".join(self.patches)])

        with open(filename, mode="w") as f:
            f.write("\n".join(mfem_file))


def compute_svg_path_stats(paths):
    stats = {
        "paths_total": len(paths),
        "curves_order_1": 0,
        "curves_order_2": 0,
        "curves_order_3": 0,
        "elliptical_arcs": 0,
        "unknown_segments": 0,
    }

    for p in paths:
        for seg in p:
            if isinstance(seg, Line):
                stats["curves_order_1"] += 1
            elif isinstance(seg, QuadraticBezier):
                stats["curves_order_2"] += 1
            elif isinstance(seg, CubicBezier):
                stats["curves_order_3"] += 1
            elif isinstance(seg, Arc):
                stats["elliptical_arcs"] += 1
            else:
                stats["unknown_segments"] += 1

    return stats


def parse_args():

    parser = argparse.ArgumentParser(
        description="svg2contours: Convert the curves in an SVG to MFEM NURBS mesh")

    parser.add_argument(
        "-i",
        "--input",
        dest="inputfile",
        required=True,
        type=argparse.FileType("r", encoding="UTF-8"),
        help="Input SVG image (*.svg)",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="outputfile",
        default="drawing.mesh",
        help="Output file in mfem NURBS mesh format (*.mesh)",
    )

    parser.add_argument(
        "--mfem-patches",
        dest="mfem_patches",
        default=False,
        action="store_true",
        help=
        ("Write the newer MFEM NURBS 'patches' mesh format for 1D segments embedded in 2D "
         "(requires MFEM > 4.9.0 to read). When enabled, Lines/Quadratic/Cubic segments are "
         "written using degree 1/2/3 respectively; elliptical arcs are written as rational quadratics."
         ),
    )

    parser.add_argument(
        "--stats",
        dest="statsfile",
        default=None,
        help="Optional output JSON file with SVG curve statistics (use '-' for stdout)",
    )

    parser.add_argument(
        "-r",
        "--reverse",
        dest="reverse_paths",
        default=False,
        action="store_true",
        help=
        "reverses paths (can be helpful during coordinate system transformation from y-axis pointing down to up)",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        default=False,
        action="store_true",
        help="verbose output flag",
    )

    opts = parser.parse_args()
    return vars(opts)


def main():
    opts = parse_args()
    verbose = opts["verbose"]

    if verbose:
        print(f"Running from '{os.getcwd()}' with arguments")
        for k, v in opts.items():
            print(f"\t{k}: {v}")

    ## Load the SVG document
    input_file = opts["inputfile"].name
    doc = Document(input_file)
    paths = doc.paths()

    if verbose:
        print(f"Attributes at root of '{input_file}':")
        for k, v in doc.tree.getroot().attrib.items():
            print(f"\t{k}: {v}")

    coordinate_transform = get_root_transform(doc)

    ## Optionally, generate a json file w/ stats about the number of curves/paths
    stats_file = opts.get("statsfile", None)
    if stats_file:
        stats = compute_svg_path_stats(paths)
        stats["input_svg"] = os.path.relpath(input_file, os.getcwd())

        root_attr = doc.tree.getroot().attrib
        parsed_viewbox = parse_viewbox(root_attr.get("viewBox", None))
        if parsed_viewbox is not None:
            min_x, min_y, vb_width, vb_height = parsed_viewbox
            stats["bounding_box"] = {
                "min_x": min_x,
                "min_y": min_y,
                "max_x": min_x + vb_width,
                "max_y": min_y + vb_height,
            }
        else:
            width_number = parse_svg_length(root_attr.get("width", None))
            height_number = parse_svg_length(root_attr.get("height", None))
            if width_number is not None and height_number is not None:
                stats["bounding_box"] = {
                    "min_x": 0.0,
                    "min_y": 0.0,
                    "max_x": width_number,
                    "max_y": height_number,
                }

        stats_json = json.dumps(stats, indent=2, sort_keys=True)
        if stats_file == "-":
            print(stats_json)
        else:
            with open(stats_file, mode="w", encoding="utf-8") as f:
                f.write(stats_json + "\n")

    ## Process SVG paths
    if verbose:
        print("SVG paths: \n", paths)

    mfem_data = MFEMData()
    mfem_patches = opts.get("mfem_patches", False)
    mfem_patches_data = MFEMPatchesData() if mfem_patches else None

    for p_idx, p in enumerate(paths):
        # print(f"""reading {p_idx=} {p=} \n w/ {p.d()=}""")

        is_d_path = "d" in p.element.keys()
        attrib = p_idx + 1

        reverse_paths = True if np.linalg.det(coordinate_transform) < 0 else False
        if opts["reverse_paths"]:
            reverse_paths = not reverse_paths

        if not all(map(is_path_segment, p)):
            continue

        for seg_idx, seg in enumerate(p):
            # print(f"""processing {seg_idx=} {seg=}""")

            if (not mfem_patches) and isinstance(seg, Arc) and seg.large_arc and is_d_path:
                # split large elliptical arcs for easier processing
                # this simplifies the derivation of the internal control points
                # in `arc_to_cubic` algorithm
                arc1, arc2 = seg.split(0.5)

                cubic, weights = segment_as_cubic(arc1, reverse_paths)
                xformed_cubic = transform_cubic(cubic, coordinate_transform)
                mfem_data.add_cubic_bezier(xformed_cubic, weights, attrib)

                cubic, weights = segment_as_cubic(arc2, reverse_paths)
                xformed_cubic = transform_cubic(cubic, coordinate_transform)
                mfem_data.add_cubic_bezier(xformed_cubic, weights, attrib)
            else:
                if mfem_patches:
                    out_segments = segment_as_native_nurbs_segments(seg, reverse_paths)
                    if isinstance(seg, Arc) and len(out_segments) > 1 and all(
                            d == 2 for (_, _, d) in out_segments):
                        quads_with_weights = []
                        for out_seg, weights0, _ in out_segments:
                            xformed_seg0 = transform_segment(out_seg, coordinate_transform)
                            quads_with_weights.append((xformed_seg0, weights0))

                        cps, weights, knots = MFEMPatchesData.quadratic_beziers_to_multispan(
                            quads_with_weights)
                        mfem_patches_data.add_nurbs_patch(
                            cps=cps,
                            degree=2,
                            weights=weights,
                            knots=knots,
                            attrib=attrib,
                            patch_comment=f"quadratic (order 2, {len(out_segments)} spans)",
                        )
                    else:
                        for out_seg, weights0, degree0 in out_segments:
                            xformed_seg0 = transform_segment(out_seg, coordinate_transform)
                            mfem_patches_data.add_bezier(xformed_seg0, degree0, weights0, attrib)
                else:
                    cubic, weights = segment_as_cubic(seg, reverse_paths)
                    xformed_cubic = transform_cubic(cubic, coordinate_transform)
                    mfem_data.add_cubic_bezier(xformed_cubic, weights, attrib)

    output_file = opts["outputfile"]
    if mfem_patches:
        mfem_patches_data.write_file(output_file)
        print(
            f"Wrote '{output_file}' with {mfem_patches_data.vert_cnt} vertices and NURBS {mfem_patches_data.elem_cnt} elements (patches format)"
        )
    else:
        mfem_data.write_file(output_file)
        print(
            f"Wrote '{output_file}' with {mfem_data.vert_cnt} vertices and NURBS {mfem_data.elem_cnt} elements"
        )


if __name__ == "__main__":
    exitcode = main()
    sys.exit(exitcode)
