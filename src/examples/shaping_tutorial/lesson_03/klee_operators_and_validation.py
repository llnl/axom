#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Loads and validates a Klee input file."""

import json
from argparse import ArgumentParser

import pyklee


def _materials_replaced_by(shape, materials):
    """Get the list of materials that the given shape will replace.

    Args:
        shape: The shape to check against.
        materials: The list of materials to check.

    Returns:
        The materials from the input list that the shape will replace.
    """

    return [material for material in materials if shape.replaces(material)]


def _dimensions_to_string(dim) -> str:
    """Format a Klee dimension enum for display."""

    if dim == pyklee.Dimensions.Two:
        return "2D"
    if dim == pyklee.Dimensions.Three:
        return "3D"
    return "<unknown>"


def _length_unit_to_string(unit) -> str:
    """Format a Klee length unit enum for display."""

    if unit == pyklee.LengthUnit.km:
        return "km"
    if unit == pyklee.LengthUnit.m:
        return "m"
    if unit == pyklee.LengthUnit.dm:
        return "dm"
    if unit == pyklee.LengthUnit.cm:
        return "cm"
    if unit == pyklee.LengthUnit.mm:
        return "mm"
    if unit == pyklee.LengthUnit.um:
        return "um"
    if unit == pyklee.LengthUnit.nm:
        return "nm"
    if unit == pyklee.LengthUnit.angstrom:
        return "A"
    if unit == pyklee.LengthUnit.miles:
        return "miles"
    if unit == pyklee.LengthUnit.feet:
        return "feet"
    if unit == pyklee.LengthUnit.inches:
        return "inches"
    if unit == pyklee.LengthUnit.mils:
        return "mils"
    return "<unspecified>"


def _formatted_units(start_units, end_units) -> str:
    """Format a pair of units"""

    start = _length_unit_to_string(start_units)
    end = _length_unit_to_string(end_units)
    return f"'{start}'" if start == end else f"'{start}' -> '{end}'"


def _formatted_dimensions(start_dims, end_dims) -> str:
    """Format a pair of dimensions"""

    start = _dimensions_to_string(start_dims)
    end = _dimensions_to_string(end_dims)
    return start if start == end else f"{start} -> {end}"


def _sorted_materials(shape_set):
    """Collect and sort the unique materials referenced by a ShapeSet."""

    return sorted({shape.getMaterial() for shape in shape_set.getShapes()})


def render_shape_set(shape_set) -> str:
    """Render a ShapeSet summary"""

    materials = _sorted_materials(shape_set)
    lines = [
        "[INFO] Klee ShapeSet Information:",
        f"  Overall dimensions: {_dimensions_to_string(shape_set.getDimensions())}",
        "",
        f"  Unique materials ({len(materials)}): {json.dumps(materials)}",
        "",
        f"  Details for the {len(shape_set.getShapes())} shapes:",
    ]

    # Print information about each shape.
    for shape in shape_set.getShapes():
        geom = shape.getGeometry()
        lines.extend([
            f"  - name: '{shape.getName()}'",
            f"    material: '{shape.getMaterial()}'",
            f"    format: '{geom.getFormat()}'",
            f"    units: {_formatted_units(geom.getStartProperties().units, geom.getEndProperties().units)}",
            f"    dimensions: {_formatted_dimensions(geom.getInputDimensions(), geom.getOutputDimensions())}",
            f"    replaces materials: {json.dumps(_materials_replaced_by(shape, materials))}",
            "",
        ])

    return "\n".join(lines) + "\n"


def main() -> int:
    """Load the input file and print the Klee ShapeSet summary."""

    parser = ArgumentParser(description="Klee Input Validator and Summary")
    parser.add_argument("input_file", help="Klee input file")
    args = parser.parse_args()

    # Parse and validate the Klee shape file, then extract summary information.
    try:
        shape_set = pyklee.readShapeSet(args.input_file)
    except RuntimeError as exc:
        # Surface parsing failures and exiting nonzero.
        print(exc)
        return 1

    print(render_shape_set(shape_set), end="")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
