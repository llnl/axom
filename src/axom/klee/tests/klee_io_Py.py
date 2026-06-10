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


def _write_mesh(path: Path) -> None:
    path.write_text(
        "MFEM mesh v1.0\n\n"
        "dimension\n2\n\n"
        "elements\n0\n\n"
        "boundary\n0\n\n"
        "vertices\n0\n",
        encoding="utf-8",
    )


def test_read_shapeset_file(tmp_path: Path):
    mesh_file = tmp_path / "shape.mesh"
    shape_file = tmp_path / "shape.yaml"
    _write_mesh(mesh_file)
    shape_file.write_text(
        f"""dimensions: 2
shapes:
  - name: wheel
    material: steel
    geometry:
      format: mfem
      path: {mesh_file.name}
      units: cm
""",
        encoding="utf-8",
    )

    shape_set = pyklee.readShapeSet(str(shape_file))

    assert shape_set.getPath() == str(shape_file)
    assert len(shape_set.getShapes()) == 1


def test_read_shapeset_shape_with_no_replacement_lists(tmp_path: Path):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    geometry:
      format: none
""",
        encoding="utf-8",
    )

    shape_set = pyklee.readShapeSet(str(shape_file))
    shape = shape_set.getShapes()[0]

    assert shape.replaces("mat1")
    assert shape.replaces("mat2")
    assert shape.getName() == "wheel"
    assert shape.getMaterial() == "steel"
    assert shape.getGeometry().getFormat() == "none"
    assert not shape.getGeometry().getPath()
    assert not shape.getGeometry().hasGeometry()
    assert shape.getGeometry().getInputDimensions() == pyklee.Dimensions.Two
    assert shape.getGeometry().getOutputDimensions() == pyklee.Dimensions.Two


def test_read_shapeset_shape_with_replaces_list(tmp_path: Path):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    replaces: [mat1, mat2]
    geometry:
      format: none
""",
        encoding="utf-8",
    )

    shape_set = pyklee.readShapeSet(str(shape_file))
    shape = shape_set.getShapes()[0]

    assert shape.replaces("mat1")
    assert shape.replaces("mat2")
    assert not shape.replaces("material_not_in_list")


def test_read_shapeset_shape_with_does_not_replace_list(tmp_path: Path):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    does_not_replace: [mat1, mat2]
    geometry:
      format: none
""",
        encoding="utf-8",
    )

    shape_set = pyklee.readShapeSet(str(shape_file))
    shape = shape_set.getShapes()[0]

    assert not shape.replaces("mat1")
    assert not shape.replaces("mat2")
    assert shape.replaces("material_not_in_list")


def test_read_shapeset_missing_name(tmp_path: Path):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes:
  - material: steel
    geometry:
      format: none
""",
        encoding="utf-8",
    )

    _assert_runtime_error_contains(lambda: pyklee.readShapeSet(str(shape_file)), "name")


def test_read_shapeset_missing_material(tmp_path: Path):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes:
  - name: wheel
    geometry:
      format: none
""",
        encoding="utf-8",
    )

    _assert_runtime_error_contains(lambda: pyklee.readShapeSet(str(shape_file)), "material")


def test_read_shapeset_missing_geometry_path(tmp_path: Path):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    geometry:
      format: test_format
""",
        encoding="utf-8",
    )

    _assert_runtime_error_contains(lambda: pyklee.readShapeSet(str(shape_file)), "Provided format")


def test_read_shapeset_format_geometry_format(tmp_path: Path):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    geometry:
      path: my/file.format
""",
        encoding="utf-8",
    )

    _assert_runtime_error_contains(lambda: pyklee.readShapeSet(str(shape_file)), "format")


def test_read_shapeset_shape_with_replaces_and_does_not_replace_lists(tmp_path: Path):
    shape_file = tmp_path / "shape.yaml"
    shape_file.write_text(
        """dimensions: 2
shapes:
  - name: wheel
    material: steel
    replaces: [mat1, mat2]
    does_not_replace: [mat1, mat2]
    geometry:
      format: none
""",
        encoding="utf-8",
    )

    _assert_runtime_error_contains(lambda: pyklee.readShapeSet(str(shape_file)), "replaces",
                                   "does_not_replace")
