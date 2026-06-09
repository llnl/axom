#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Parse, validate, and print 2D or 3D mesh metadata from structured input."""

from argparse import ArgumentParser
from dataclasses import dataclass, field

import pyinlet

MAX_INT = 2 ** 31 - 1


@dataclass
class MeshMetadata:
    """Mesh metadata matching the lesson's C++ structure."""

    @dataclass
    class BoundingBox:
        min_corner: tuple[float, ...]
        max_corner: tuple[float, ...]

        def range(self) -> tuple[float, ...]:
            """Return the extents of the bounding box."""

            return tuple(max_value - min_value
                         for min_value, max_value in zip(self.min_corner, self.max_corner))

        def __str__(self) -> str:
            return f"{{ min:{self.min_corner}; max:{self.max_corner}; range:{self.range()} }}"

    dim: int = 2
    bb_min: list[float] = field(default_factory=lambda: [0.0, 0.0])
    bb_max: list[float] = field(default_factory=lambda: [1.0, 1.0])
    resolution: list[int] = field(default_factory=lambda: [10, 10])

    def bounding_box(self) -> BoundingBox:
        """Return the mesh bounding box."""

        return MeshMetadata.BoundingBox(tuple(self.bb_min), tuple(self.bb_max))


def define_schema(mesh_schema) -> None:
    """Define the mesh schema directly in Python."""

    mesh_schema.addInt("dim", "Dimension (2 or 3)").required().range(2, 3)

    bounding_box = mesh_schema.addStruct("bounding_box", "Mesh bounding box").required()

    minimum = bounding_box.addStruct("min", "Minimum coordinates").required()
    minimum.addDouble("x", "Minimum x coordinate").required()
    minimum.addDouble("y", "Minimum y coordinate").required()
    minimum.addDouble("z", "Minimum z coordinate (only specify when dim is 3)")

    maximum = bounding_box.addStruct("max", "Maximum coordinates").required()
    maximum.addDouble("x", "Maximum x coordinate").required()
    maximum.addDouble("y", "Maximum y coordinate").required()
    maximum.addDouble("z", "Maximum z coordinate (only specify when dim is 3)")

    resolution = mesh_schema.addStruct("resolution", "Mesh resolution").required()
    resolution.addInt("x", "Resolution in x direction").required().range(1, MAX_INT)
    resolution.addInt("y", "Resolution in y direction").required().range(1, MAX_INT)
    resolution.addInt("z",
                      "Resolution in z direction (only specify when dim is 3)").range(1, MAX_INT)

    def verify_bounding_box(input_data) -> bool:
        valid = True
        for axis in ("x", "y", "z"):
            min_path = f"min/{axis}"
            max_path = f"max/{axis}"
            if axis == "z" and (not input_data.contains(min_path)
                                or not input_data.contains(max_path)):
                continue

            min_value = float(input_data[min_path])
            max_value = float(input_data[max_path])
            if min_value >= max_value:
                print(f"Invalid bounding box range for {axis}-coordinate: "
                      f"{min_value} >= {max_value}")
                valid = False
        return valid

    def verify_dimension(input_data) -> bool:
        dim = int(input_data["dim"])
        valid = True

        for field in ("bounding_box/min/z", "bounding_box/max/z", "resolution/z"):
            if dim == 3 and not input_data.contains(field):
                print(f"Z-coordinate for '{field}' is required when dimension is 3")
                valid = False
            elif dim == 2 and input_data.contains(field):
                print(f"Z-coordinate for '{field}' should not be provided when dimension is 2")
                valid = False
        return valid

    bounding_box.registerVerifier(verify_bounding_box)
    mesh_schema.registerVerifier(verify_dimension)


def mesh_metadata_from_proxy(mesh) -> MeshMetadata:
    """Construct a MeshMetadata instance from a verified mesh proxy."""

    metadata = MeshMetadata(dim=int(mesh["dim"]))

    bounding_box = mesh["bounding_box"]
    resolution = mesh["resolution"]

    metadata.bb_min = [float(bounding_box["min/x"]), float(bounding_box["min/y"])]
    metadata.bb_max = [float(bounding_box["max/x"]), float(bounding_box["max/y"])]
    metadata.resolution = [int(resolution["x"]), int(resolution["y"])]

    if metadata.dim == 3:
        metadata.bb_min.append(float(bounding_box["min/z"]))
        metadata.bb_max.append(float(bounding_box["max/z"]))
        metadata.resolution.append(int(resolution["z"]))

    return metadata


def load_mesh_metadata(input_file: str) -> MeshMetadata:
    """Load validated 2D or 3D mesh metadata using the C++-style workflow."""

    if input_file.endswith((".yaml", ".yml")):
        reader = pyinlet.YAMLReader()
    else:
        reader = pyinlet.LuaReader()

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

    bounding_box = metadata.bounding_box()

    print("Parsed mesh metadata:")
    print(f"  - dimension: {metadata.dim}")
    print(f"  - bounding Box: {bounding_box}")
    print(f"  - resolution: {metadata.resolution}")


def main() -> int:
    """Run the improved inlet metadata example.

    Returns:
        Process exit status.
    """

    parser = ArgumentParser(description="Inlet Metadata Setup")
    parser.add_argument("input_file", help="YAML or Lua input file with inlet metadata")
    args = parser.parse_args()

    supported_extensions = [".yaml", ".yml"]
    if hasattr(pyinlet, "LuaReader"):
        supported_extensions.append(".lua")

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
