#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

from pathlib import Path

import pyklee


def _read_shape_set(tmp_path: Path, yaml_contents: str):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(yaml_contents, encoding="utf-8")
    return pyklee.readShapeSet(str(shape_file))


def _assert_runtime_error_contains(func, *expected_substrings):
    try:
        func()
    except RuntimeError as exc:
        message = str(exc)
        for substring in expected_substrings:
            assert substring in message
    else:
        assert False, "Expected RuntimeError"


def test_replaces_no_lists_given(tmp_path: Path):
    shape_set = _read_shape_set(
        tmp_path,
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    geometry:
      format: none
""",
    )

    shape = shape_set.getShapes()[0]
    assert shape.replaces("some_material")


def test_replaces_replacement_list_given(tmp_path: Path):
    shape_set = _read_shape_set(
        tmp_path,
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    replaces: [mat1, mat2]
    geometry:
      format: none
""",
    )

    shape = shape_set.getShapes()[0]
    assert shape.replaces("mat1")
    assert shape.replaces("mat2")
    assert not shape.replaces("mat3")


def test_replaces_non_replacement_list_given(tmp_path: Path):
    shape_set = _read_shape_set(
        tmp_path,
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    does_not_replace: [mat1, mat2]
    geometry:
      format: none
""",
    )

    shape = shape_set.getShapes()[0]
    assert not shape.replaces("mat1")
    assert not shape.replaces("mat2")
    assert shape.replaces("mat3")


def test_both_replacement_lists_given(tmp_path: Path):
    def read_invalid():
        _read_shape_set(
            tmp_path,
            """dimensions: 2
shapes:
  - name: wheel
    material: steel
    replaces: [replaced]
    does_not_replace: [not_replaced]
    geometry:
      format: none
""",
        )

    _assert_runtime_error_contains(read_invalid, "replaces", "does_not_replace")
