#!/bin/sh
"exec" "python3" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: svg2contours.py

 description: 
  Reads in an SVG document and outputs an MFEM NURBS mesh.
  Depends on the svgpathtools module
"""

import sys
import os
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


def get_root_transform(doc: Document):
    """Create transform to convert 2D coordinate system from y pointing down to y pointing up"""

    attr = doc.tree.getroot().attrib
    width = attr.get("width", None)
    height = attr.get("height", None)
    viewBox = attr.get("viewBox", None)

    print(f"""SVG dimensions: {width=} {height=} {viewBox=}""")

    tf = np.identity(3)

    # Prefer viewBox since it defines the internal user-coordinate system used by path data.
    if viewBox:
        try:
            # SVG spec allows whitespace and/or comma separators.
            parts = [p for p in re.split(r"[,\s]+", viewBox.strip()) if p]
            if len(parts) == 4:
                _, min_y, _, vb_height = map(float, parts)
                tf[1, 1] = -1
                tf[1, 2] = 2.0 * min_y + vb_height
                return tf
        except Exception:
            pass

    # Fallback to root height attribute when viewBox is unavailable/invalid.
    if height and not height.endswith("%"):
        match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", height.strip())
        if match:
            height_number = float(match.group(0))

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


def dist_to_ellipse(center, radius, angle, pt):
    cx, cy = center.real, center.imag
    rx, ry = radius.real, radius.imag

    rot = np.exp(-1j * np.radians(angle))
    transformed_pt = rot * complex(pt.real - cx, pt.imag - cy)
    return transformed_pt.real**2 / rx**2 + transformed_pt.imag**2 / ry**2 - 1


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
        self.dof_ends.append(" ".join(map(str, [cubic.start.real, cubic.start.imag, cubic.end.real, cubic.end.imag])))
        self.dof_ints.append(
            " ".join(map(str, [cubic.control2.real, cubic.control2.imag, cubic.control1.real, cubic.control1.imag]))
        )

    def write_file(self, filename):
        mfem_file = []

        mfem_file.extend(
            [
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
            ]
        )

        mfem_file.extend(["dimension", "1", ""])

        mfem_file.extend(["elements", f"{self.elem_cnt}", "\n".join(self.elems), ""])

        mfem_file.extend(["boundary", "0", ""])

        mfem_file.extend(["edges", f"{self.elem_cnt}", "\n".join(self.edges), ""])

        mfem_file.extend(["vertices", f"{self.vert_cnt}", ""])

        mfem_file.extend(["knotvectors", f"{self.elem_cnt}", "\n".join(self.knots), ""])

        mfem_file.extend(["weights", "\n".join(self.wgts_ends), "\n".join(self.wgts_ints), ""])

        mfem_file.extend(
            [
                "FiniteElementSpace",
                "FiniteElementCollection: NURBS",
                "VDim: 2",
                "Ordering: 1",
                "",
                "\n".join(self.dof_ends),
                "\n".join(self.dof_ints),
                "",
            ]
        )

        with open(filename, mode="w") as f:
            f.write("\n".join(mfem_file))


def parse_args():

    parser = argparse.ArgumentParser(description="svg2contours: Convert the curves in an SVG to MFEM NURBS mesh")

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
        "-r", "--reverse", dest="reverse_paths", default=False, action="store_true", 
        help="reverses paths (can be helpful during coordinate system transformation from y-axis pointing down to up)"
    )

    parser.add_argument(
        "-v", "--verbose", dest="verbose", default=False, action="store_true", help="verbose output flag"
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

    ## Process SVG paths
    if verbose:
        print("SVG paths: \n", paths)

    mfem_data = MFEMData()

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

            if isinstance(seg, Arc) and seg.large_arc and is_d_path:
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
                cubic, weights = segment_as_cubic(seg, reverse_paths)
                xformed_cubic = transform_cubic(cubic, coordinate_transform)
                mfem_data.add_cubic_bezier(xformed_cubic, weights, attrib)

    output_file = opts["outputfile"]
    mfem_data.write_file(output_file)
    print(f"Wrote '{output_file}' with {mfem_data.vert_cnt} vertices and NURBS {mfem_data.elem_cnt} elements")


if __name__ == "__main__":
    exitcode = main()
    sys.exit(exitcode)
