# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""Helpers for locating Axom Python wheel configuration files."""

from __future__ import annotations

import argparse
from importlib import resources
from pathlib import Path


def share_dir() -> Path:
    """Return the installed Axom Python package share directory."""
    return Path(resources.files("axom").joinpath("share"))


def host_config_path() -> Path:
    """Return the CMake host-config generated for this Axom Python wheel."""
    return share_dir() / "axom-python-host-config.cmake"


def env_script_path() -> Path:
    """Return the shell environment helper generated for this Axom Python wheel."""
    return share_dir() / "axom-python-env.sh"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Print Axom Python wheel configuration paths.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--host-config", action="store_true", help="print the CMake host-config path")
    group.add_argument("--env-script", action="store_true", help="print the shell environment script path")
    group.add_argument("--share-dir", action="store_true", help="print the Axom Python share directory")
    group.add_argument("--cmake-args", action="store_true", help="print CMake arguments using the host-config")
    args = parser.parse_args(argv)

    if args.env_script:
        print(env_script_path())
    elif args.share_dir:
        print(share_dir())
    elif args.cmake_args:
        print(f"-C {host_config_path()}")
    else:
        print(host_config_path())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
