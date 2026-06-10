#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Shape a Klee input onto a computational mesh using Quest SamplingShaper."""

from argparse import ArgumentParser
from dataclasses import dataclass, field
from pathlib import Path

import pyklee
import pyinlet
import pyquest

MAX_INT = 2 ** 31 - 1


@dataclass
class MeshMetadata:
    dim: int = 2
    bb_min: list[float] = field(default_factory=lambda: [0.0, 0.0])
    bb_max: list[float] = field(default_factory=lambda: [1.0, 1.0])
    resolution: list[int] = field(default_factory=lambda: [10, 10])
    background_material: str = ""
    volume_fraction_order: int = 2
    mesh_order: int = 1
    sampling_resolution: int = 5
    sampling_method: pyquest.SamplingMethod = pyquest.SamplingMethod.InOut


def _make_reader(path: str):
    suffix = Path(path).suffix.lower()
    if suffix in {".yaml", ".yml"}:
        return pyinlet.YAMLReader()
    if suffix == ".lua" and hasattr(pyinlet, "LuaReader"):
        return pyinlet.LuaReader()

    supported = [".yaml", ".yml"]
    if hasattr(pyinlet, "LuaReader"):
        supported.append(".lua")

    raise RuntimeError(f"Unsupported mesh metadata extension for '{path}'. "
                       f"Expected {', '.join(supported[:-1])} or {supported[-1]}.")


def _define_mesh_schema(mesh_schema) -> None:
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

    mesh_schema.addString("background_material", "Optional background material")
    mesh_schema.addInt("volume_fraction_order",
                       "Order for volume fraction fields (>= 1)").range(1, MAX_INT)
    mesh_schema.addInt("mesh_order", "Order for mesh nodes (>= 1)").range(1, MAX_INT)
    mesh_schema.addInt("sampling_resolution", "Sampling resolution (>= 1)").range(1, MAX_INT)
    mesh_schema.addInt("quadrature_order",
                       "Legacy alias for sampling_resolution (>= 1)").range(1, MAX_INT)
    mesh_schema.addString("sampling_method", "Sampling method ('inout' or 'winding')").validValues(
        ["inout", "winding"])

    def verify_bounding_box(input_data) -> bool:
        valid = True
        for axis in ("x", "y", "z"):
            min_name = f"min/{axis}"
            max_name = f"max/{axis}"
            if axis == "z" and (not input_data.contains(min_name)
                                and not input_data.contains(max_name)):
                continue

            min_value = float(input_data[min_name])
            max_value = float(input_data[max_name])
            if min_value >= max_value:
                print(
                    f"Invalid bounding box range for {axis}-coordinate: {min_value} >= {max_value}")
                valid = False

        return valid

    def verify_dimension(input_data) -> bool:
        dim = int(input_data["dim"])
        valid = True

        for field_name in ("bounding_box/min/z", "bounding_box/max/z", "resolution/z"):
            if dim == 3:
                if not input_data.contains(field_name):
                    print(f"Z-coordinate for '{field_name}' is required when dimension is 3")
                    valid = False
            elif dim == 2 and input_data.contains(field_name):
                print(f"Z-coordinate for '{field_name}' should not be provided when dimension is 2")
                valid = False

        return valid

    bounding_box.registerVerifier(verify_bounding_box)
    mesh_schema.registerVerifier(verify_dimension)


def _mesh_metadata_from_inlet(mesh) -> MeshMetadata:
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

    if mesh.contains("background_material"):
        metadata.background_material = str(mesh["background_material"])
    if mesh.contains("volume_fraction_order"):
        metadata.volume_fraction_order = int(mesh["volume_fraction_order"])
    if mesh.contains("mesh_order"):
        metadata.mesh_order = int(mesh["mesh_order"])
    if mesh.contains("sampling_resolution"):
        metadata.sampling_resolution = int(mesh["sampling_resolution"])
    elif mesh.contains("quadrature_order"):
        metadata.sampling_resolution = int(mesh["quadrature_order"])

    if mesh.contains("sampling_method"):
        method = str(mesh["sampling_method"])
        metadata.sampling_method = (pyquest.SamplingMethod.WindingNumber
                                    if method == "winding" else pyquest.SamplingMethod.InOut)

    return metadata


def _read_mesh_metadata(path: str) -> MeshMetadata:
    reader = _make_reader(path)
    if not reader.parseFile(path):
        raise RuntimeError(f"Failed to parse '{path}'.")

    inlet = pyinlet.Inlet(reader)
    mesh_schema = inlet.addStruct("mesh", "Mesh metadata").required()
    _define_mesh_schema(mesh_schema)

    if not inlet.verify():
        raise RuntimeError("Mesh metadata validation failed.")

    return _mesh_metadata_from_inlet(inlet["mesh"])


def main() -> int:
    parser = ArgumentParser(
        description="Shaping pipeline using separate Inlet mesh metadata and Klee shapes.")
    parser.add_argument(
        "-m",
        "--mesh_file",
        required=True,
        help="Mesh metadata Inlet Lua or YAML file.",
    )
    parser.add_argument(
        "-k",
        "--klee_file",
        required=True,
        help="Klee shape-set YAML file.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging.")
    args = parser.parse_args()

    try:
        meta = _read_mesh_metadata(args.mesh_file)
        shape_set = pyklee.readShapeSet(args.klee_file)
    except RuntimeError as exc:
        print(exc)
        return 1

    shaper = pyquest.SamplingShaper(
        shape_set,
        meta.dim,
        meta.bb_min,
        meta.bb_max,
        meta.resolution,
        meta.mesh_order,
    )
    shaper.setVerbosity(args.verbose)
    shaper.setSamplingResolution(meta.sampling_resolution)
    shaper.setVolumeFractionOrder(meta.volume_fraction_order)
    shaper.setSamplingMethod(meta.sampling_method)

    if meta.sampling_method == pyquest.SamplingMethod.InOut:
        shaper.setSamplesPerKnotSpan(50)

    if meta.background_material:
        shaper.importBackgroundMaterial(meta.background_material, meta.volume_fraction_order)

    for shape in shape_set.getShapes():
        shape_dim = shape.getGeometry().getInputDimensions()
        shaper.loadShape(shape)
        shaper.prepareShapeQuery(shape_dim, shape)
        shaper.runShapeQuery(shape)
        shaper.applyReplacementRules(shape)
        shaper.finalizeShapeQuery()

    shaper.adjustVolumeFractions()
    shaper.save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
