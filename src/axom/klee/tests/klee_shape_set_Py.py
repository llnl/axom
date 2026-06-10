#!/usr/bin/env python3

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

from pathlib import Path

import pyklee


def _assert_runtime_error_contains(func, *expected_substrings):
    try:
        func()
    except RuntimeError as exc:
        message = str(exc)
        for substring in expected_substrings:
            assert substring in message
    else:
        assert False, "Expected RuntimeError"


def test_read_shapeset_no_shapes(tmp_path: Path):
    shape_file = tmp_path / "no_shapes.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes: []
""",
        encoding="utf-8",
    )

    shape_set = pyklee.readShapeSet(str(shape_file))

    assert shape_set.getShapes() == []
    assert shape_set.getDimensions() == pyklee.Dimensions.Two
    assert shape_set.getDimensions() != pyklee.Dimensions.Three
    assert shape_set.getDimensions() != pyklee.Dimensions.Unspecified


def test_read_shapeset_invalid_dimensions(tmp_path: Path):
    shape_file = tmp_path / "invalid_dimensions.yaml"
    shape_file.write_text(
        """dimensions: 5
shapes: []
""",
        encoding="utf-8",
    )

    _assert_runtime_error_contains(lambda: pyklee.readShapeSet(str(shape_file)), "dimensions")
