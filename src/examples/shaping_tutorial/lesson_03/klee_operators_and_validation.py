#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path

import pyklee


def main() -> int:
    parser = ArgumentParser(description="Load and summarize a Klee shape-set file.")
    parser.add_argument("input_file", nargs="?", default=str(Path(__file__).with_name("ice_cream.yaml")))
    args = parser.parse_args()

    shape_set = pyklee.readShapeSet(args.input_file)
    print("Overall dimensions:", shape_set.getDimensions())
    print("Number of shapes:", len(shape_set.getShapes()))
    for shape in shape_set.getShapes():
        geom = shape.getGeometry()
        print(
            f"- {shape.getName()}: material={shape.getMaterial()}, "
            f"format={geom.getFormat()}, dims={geom.getInputDimensions()}->{geom.getOutputDimensions()}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
