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


def test_proxy_get_mesh_metadata(tmp_path: Path):
    input_file = tmp_path / "mesh.yaml"
    input_file.write_text(
        """mesh:
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
""",
        encoding="utf-8",
    )

    reader = pyinlet.YAMLReader()
    reader.parseFile(str(input_file))
    inlet = pyinlet.Inlet(reader)
    mesh_schema = inlet.addStruct("mesh", "Mesh metadata").required()
    pyinlet.define_validated_mesh_schema(mesh_schema)

    assert inlet.verify()

    meta = inlet["mesh"].getMeshMetadata()
    assert meta.dim == 2
    assert meta.bb_min == [0.0, 1.0]
    assert meta.bb_max == [2.0, 3.0]
    assert meta.resolution == [8, 16]


def test_python_defined_container_verifier(tmp_path: Path):
    input_file = tmp_path / "mesh.yaml"
    input_file.write_text(
        """mesh:
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
""",
        encoding="utf-8",
    )

    reader = pyinlet.YAMLReader()
    reader.parseFile(str(input_file))
    inlet = pyinlet.Inlet(reader)
    mesh_schema = inlet.addStruct("mesh", "Mesh metadata").required()

    bounding_box = mesh_schema.addStruct("bounding_box", "Mesh bounding box").required()

    minimum = bounding_box.addStruct("min", "Minimum coordinates").required()
    minimum.addDouble("x", "Minimum x coordinate").required()
    minimum.addDouble("y", "Minimum y coordinate").required()

    maximum = bounding_box.addStruct("max", "Maximum coordinates").required()
    maximum.addDouble("x", "Maximum x coordinate").required()
    maximum.addDouble("y", "Maximum y coordinate").required()

    resolution = mesh_schema.addStruct("resolution", "Mesh resolution").required()
    resolution.addInt("x", "Resolution in x direction").required().range(1, 2 ** 31 - 1)
    resolution.addInt("y", "Resolution in y direction").required().range(1, 2 ** 31 - 1)

    bounding_box.registerVerifier(lambda container: float(container["min/x"]) < float(container[
        "max/x"]) and float(container["min/y"]) < float(container["max/y"]))

    assert inlet.verify()


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
