#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path

import pyquest


def main() -> int:
    lesson_dir = Path(__file__).resolve().parent

    parser = ArgumentParser(description="Run the shaping tutorial pipeline from Python.")
    parser.add_argument(
        "-m",
        "--mesh_file",
        default=str(lesson_dir / "circle_input.lua"),
        help="Mesh metadata file in YAML or Lua format.",
    )
    parser.add_argument(
        "-k",
        "--klee_file",
        default=str(lesson_dir / "circles.yaml"),
        help="Klee shape-set file.",
    )
    parser.add_argument(
        "-o",
        "--output_name",
        default="py_shaping",
        help="Output collection name for saved Sidre data.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    result = pyquest.runShaping(args.mesh_file, args.klee_file, args.output_name, args.verbose)
    result.save()
    print("Saved shaping output as:", result.getCollectionName())
    print("Blueprint group:", result.getBlueprintGroup().getPathName())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
