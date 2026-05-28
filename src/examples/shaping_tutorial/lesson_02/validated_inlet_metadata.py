#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path

import pyinlet


def main() -> int:
    parser = ArgumentParser(description="Validate and print 2D mesh metadata.")
    parser.add_argument("input_file",
                        nargs="?",
                        default=str(Path(__file__).with_name("input2D.yaml")))
    args = parser.parse_args()

    errors = pyinlet.validate_mesh_metadata(args.input_file)
    if errors:
        for err in errors:
            print(err)
        return 1

    meta = pyinlet.load_mesh_metadata(args.input_file)
    bbox = meta.getBoundingBox2D()
    print("Bounding Box Min:", tuple(bbox.getMin()))
    print("Bounding Box Max:", tuple(bbox.getMax()))
    print("Resolution:", tuple(meta.resolution))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
