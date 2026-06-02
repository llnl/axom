#!/usr/bin/env python3

import sys

import pycore


def main() -> int:
    print("Python executable:", sys.executable)
    # Use the pycore nanobind binding to show Axom configuration
    pycore.about()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
