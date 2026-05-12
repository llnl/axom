#!/usr/bin/env python3

from importlib import import_module
from pathlib import Path
import sys


def describe_module(name: str) -> str:
    module = import_module(name)
    module_path = Path(getattr(module, "__file__", "<builtin>"))
    return f"{name}: {module_path}"


def main() -> int:
    print("Python executable:", sys.executable)
    for module_name in ["pysidre", "pyprimal", "pyinlet", "pyklee", "pyquest"]:
        print(describe_module(module_name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
