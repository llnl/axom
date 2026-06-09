#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Parse, validate, and print 2D mesh metadata from structured input."""

from argparse import ArgumentParser
from dataclasses import dataclass

import pyinlet

MAX_INT = 2 ** 31 - 1


@dataclass
class MeshMetadata:
    """Mesh metadata matching the lesson's C++ structure."""

    @dataclass
    class BoundingBox:
        min_x: float
        min_y: float
        max_x: float
        max_y: float

    @dataclass
    class Resolution:
        x: int
        y: int

    bounding_box: BoundingBox
    resolution: Resolution


def define_schema(mesh_schema) -> None:
    """Define the mesh schema directly in Python."""

    bounding_box = mesh_schema.addStruct("bounding_box", "Mesh bounding box").required()

    minimum = bounding_box.addStruct("min", "Minimum coordinates").required()
    minimum.addDouble("x", "Minimum x coordinate").required()
    minimum.addDouble("y", "Minimum y coordinate").required()

    maximum = bounding_box.addStruct("max", "Maximum coordinates").required()
    maximum.addDouble("x", "Maximum x coordinate").required()
    maximum.addDouble("y", "Maximum y coordinate").required()

    resolution = mesh_schema.addStruct("resolution", "Mesh resolution").required()
    resolution.addInt("x", "Resolution in x direction").required().range(1, MAX_INT)
    resolution.addInt("y", "Resolution in y direction").required().range(1, MAX_INT)

    def verify_bounding_box(input_data) -> bool:
        min_x = float(input_data["min/x"])
        max_x = float(input_data["max/x"])
        min_y = float(input_data["min/y"])
        max_y = float(input_data["max/y"])

        valid = True
        if min_x >= max_x:
            print(f"Invalid bounding box range for x-coordinate: {min_x} >= {max_x}")
            valid = False
        if min_y >= max_y:
            print(f"Invalid bounding box range for y-coordinate: {min_y} >= {max_y}")
            valid = False
        return valid

    bounding_box.registerVerifier(verify_bounding_box)


def mesh_metadata_from_proxy(mesh) -> MeshMetadata:
    """Construct a MeshMetadata instance from a verified mesh proxy."""

    bounding_box = mesh["bounding_box"]
    resolution = mesh["resolution"]

    return MeshMetadata(
        bounding_box=MeshMetadata.BoundingBox(
            min_x=float(bounding_box["min/x"]),
            min_y=float(bounding_box["min/y"]),
            max_x=float(bounding_box["max/x"]),
            max_y=float(bounding_box["max/y"]),
        ),
        resolution=MeshMetadata.Resolution(
            x=int(resolution["x"]),
            y=int(resolution["y"]),
        ),
    )


def load_mesh_metadata(input_file: str) -> MeshMetadata:
    """Load validated 2D mesh metadata using the same flow as the C++ example."""

    reader = pyinlet.YAMLReader()
    reader.parseFile(input_file)
    inlet = pyinlet.Inlet(reader)
    mesh_schema = inlet.addStruct("mesh", "Mesh metadata").required()
    define_schema(mesh_schema)

    if not inlet.verify():
        print("Error: Input validation failed.")
        print("Missing required fields or invalid data.")
        raise ValueError("mesh metadata validation failed")

    return mesh_metadata_from_proxy(inlet["mesh"])


def print_metadata(metadata) -> None:
    """Print the parsed mesh metadata.

    Args:
        metadata: The parsed mesh metadata object.

    Returns:
        None.
    """

    print(f"Bounding Box Min: ({metadata.bounding_box.min_x}, {metadata.bounding_box.min_y})")
    print(f"Bounding Box Max: ({metadata.bounding_box.max_x}, {metadata.bounding_box.max_y})")
    print(f"Resolution: ({metadata.resolution.x}, {metadata.resolution.y})")


def main() -> int:
    """Run the validated inlet metadata example.

    Returns:
        Process exit status.
    """

    parser = ArgumentParser(description="Inlet Metadata Setup")
    parser.add_argument("input_file", help="YAML input file with inlet metadata")
    args = parser.parse_args()

    supported_extensions = [".yaml", ".yml"]

    if not any(args.input_file.endswith(ext) for ext in supported_extensions):
        supported = ", ".join(f'"{ext}"' for ext in supported_extensions)
        print(f"input_file: Invalid extension for file '{args.input_file}'; "
              f"supported extensions: [{supported}]")
        return 1

    try:
        metadata = load_mesh_metadata(args.input_file)
    except ValueError:
        return 1

    print_metadata(metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
