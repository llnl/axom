import pyprimal


def test_bbox_2d_basic_properties():
    bbox = pyprimal.BoundingBox2D.fromPoints([0.0, 1.0], [2.0, 3.0])

    assert bbox.isValid()
    assert bbox.getMin() == [0.0, 1.0]
    assert bbox.getMax() == [2.0, 3.0]
    assert bbox.contains([1.0, 2.0])
    assert not bbox.contains([3.0, 2.0])


def test_bbox_3d_intersection_and_containment():
    outer = pyprimal.BoundingBox3D.fromPoints([0.0, 0.0, 0.0], [4.0, 4.0, 4.0])
    inner = pyprimal.BoundingBox3D.fromPoints([1.0, 1.0, 1.0], [2.0, 2.0, 2.0])
    disjoint = pyprimal.BoundingBox3D.fromPoints([5.0, 5.0, 5.0], [6.0, 6.0, 6.0])

    assert outer.containsBox(inner)
    assert outer.intersectsWith(inner)
    assert not outer.intersectsWith(disjoint)
