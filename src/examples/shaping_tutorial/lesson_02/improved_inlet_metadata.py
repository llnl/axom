#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path

import pyinlet


def main() -> int:
    parser = ArgumentParser(description="Validate and print 2D or 3D mesh metadata.")
    parser.add_argument("input_file", nargs="?", default=str(Path(__file__).with_name("input3D.yaml")))
    args = parser.parse_args()

    meta = pyinlet.load_mesh_metadata(args.input_file)
    if meta.dim == 2:
        bbox = meta.getBoundingBox2D()
    else:
        bbox = meta.getBoundingBox3D()

    print("Dimension:", meta.dim)
    print("Bounding Box Min:", tuple(bbox.getMin()))
    print("Bounding Box Max:", tuple(bbox.getMax()))
    print("Resolution:", tuple(meta.resolution))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
