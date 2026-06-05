#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""Parse, validate, and print 2D or 3D mesh metadata from structured input."""

from argparse import ArgumentParser

import pyinlet


def print_metadata(metadata) -> None:
    """Print the parsed mesh metadata.

    Args:
        metadata: The parsed mesh metadata object.

    Returns:
        None.
    """

    print(f"Dimension: {metadata.dim}")
    print(f"Bounding Box Min: {tuple(metadata.bb_min)}")
    print(f"Bounding Box Max: {tuple(metadata.bb_max)}")
    print(f"Resolution: {tuple(metadata.resolution)}")


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
        print(
            f"input_file: Invalid extension for file '{args.input_file}'; "
            f"supported extensions: [{supported}]"
        )
        return 1

    metadata = pyinlet.load_mesh_metadata(args.input_file)
    print_metadata(metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
