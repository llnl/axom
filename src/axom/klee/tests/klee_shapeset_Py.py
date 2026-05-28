from pathlib import Path

import pyklee


def test_read_shapeset(tmp_path: Path):
    geometry_file = tmp_path / "shape.mfem"
    geometry_file.write_text(
        "MFEM mesh v1.0\n\ndimension\n2\n\nelements\n0\n\nboundary\n0\n\nvertices\n0\n",
        encoding="utf-8")

    shape_file = tmp_path / "shapes.yaml"
    shape_file.write_text(
        f"""dimensions: 2
shapes:
  - name: test_shape
    material: steel
    geometry:
      format: mfem
      path: {geometry_file.name}
      units: cm
""",
        encoding="utf-8",
    )

    shape_set = pyklee.readShapeSet(str(shape_file))

    assert shape_set.getDimensions() == pyklee.Dimensions.Two
    shape = shape_set.getShapes()[0]
    assert shape.getName() == "test_shape"
    assert shape.getMaterial() == "steel"
    assert shape.getGeometry().getFormat() == "mfem"
