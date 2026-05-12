from pathlib import Path

import pyinlet


def test_load_shaping_metadata_yaml(tmp_path: Path):
    input_file = tmp_path / "mesh.yaml"
    input_file.write_text(
        """mesh:
  dim: 2
  bounding_box:
    min:
      x: 0.0
      y: 1.0
    max:
      x: 2.0
      y: 3.0
  resolution:
    x: 8
    y: 16
  sampling_method: inout
""",
        encoding="utf-8",
    )

    meta = pyinlet.load_mesh_metadata(str(input_file))
    bbox = meta.getBoundingBox2D()

    assert meta.dim == 2
    assert meta.resolution == [8, 16]
    assert bbox.getMin() == [0.0, 1.0]
    assert bbox.getMax() == [2.0, 3.0]


def test_validate_shaping_metadata_reports_errors(tmp_path: Path):
    input_file = tmp_path / "bad_mesh.yaml"
    input_file.write_text(
        """mesh:
  dim: 2
  bounding_box:
    min:
      x: 3.0
      y: 1.0
    max:
      x: 2.0
      y: 4.0
  resolution:
    x: 8
    y: 16
    z: 32
""",
        encoding="utf-8",
    )

    errors = pyinlet.validate_mesh_metadata(str(input_file))
    assert errors
