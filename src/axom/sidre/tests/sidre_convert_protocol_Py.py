# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import pysidre


def _build_dir() -> Path:
    return Path(pysidre.__file__).resolve().parents[1]


def _repo_dir() -> Path:
    return _build_dir().parent


def _run_tool(command, cwd):
    result = subprocess.run(
        command,
        cwd=cwd,
        env=os.environ.copy(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout


def _compare_conversion(strip_value=None):
    build_dir = _build_dir()
    repo_dir = _repo_dir()
    input_root = repo_dir / "data" / "quest" / "box_2D_r3.root"
    cpp_tool = build_dir / "bin" / "convert_sidre_protocol"
    py_tool = repo_dir / "src" / "tools" / "convert_sidre_protocol.py"

    with tempfile.TemporaryDirectory() as cpp_tmp, tempfile.TemporaryDirectory() as py_tmp:
        cpp_args = [
            str(cpp_tool),
            "--input",
            str(input_root),
            "--output",
            "converted",
        ]
        py_args = [
            sys.executable,
            str(py_tool),
            "--input",
            str(input_root),
            "--output",
            "converted",
        ]

        if strip_value is not None:
            cpp_args.extend(["--strip", str(strip_value)])
            py_args.extend(["--strip", str(strip_value)])

        _run_tool(cpp_args, cpp_tmp)
        _run_tool(py_args, py_tmp)

        assert (Path(cpp_tmp) / "converted.root").read_bytes() == (
            Path(py_tmp) / "converted.root"
        ).read_bytes()
        assert (Path(cpp_tmp) / "converted_0000000.json").read_bytes() == (
            Path(py_tmp) / "converted_0000000.json"
        ).read_bytes()


def test_convert_sidre_protocol_matches_cpp():
    _compare_conversion()


def test_convert_sidre_protocol_strip_matches_cpp():
    _compare_conversion(strip_value=3)
