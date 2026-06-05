#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""Parse, validate, and print 2D mesh metadata from structured input."""

from argparse import ArgumentParser

import pyinlet


def print_metadata(metadata) -> None:
    """Print the parsed mesh metadata.

    Args:
        metadata: The parsed mesh metadata object.

    Returns:
        None.
    """

    print(f"Bounding Box Min: ({metadata.bb_min[0]}, {metadata.bb_min[1]})")
    print(f"Bounding Box Max: ({metadata.bb_max[0]}, {metadata.bb_max[1]})")
    print(f"Resolution: ({metadata.resolution[0]}, {metadata.resolution[1]})")


def main() -> int:
    """Run the validated inlet metadata example.

    Returns:
        Process exit status.
    """

    parser = ArgumentParser(description="Inlet Metadata Setup")
    parser.add_argument("input_file", help="YAML input file with inlet metadata")
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

    errors = pyinlet.validate_mesh_metadata(args.input_file)
    if errors:
        print("Error: Input validation failed.")
        print("Missing required fields or invalid data.")
        return 1

    metadata = pyinlet.load_mesh_metadata(args.input_file)
    print_metadata(metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
