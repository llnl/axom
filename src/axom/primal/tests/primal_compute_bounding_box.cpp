// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <limits>
#include <algorithm>

#include "gtest/gtest.h"

#include "axom/core/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/operators/compute_bounding_box.hpp"

namespace primal = axom::primal;

TEST(primal_compute_bounding_box, compute_oriented_box_test)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  QPoint pt1;      // origin
  QVector u[DIM];  // make standard axes
  for(int i = 0; i < DIM; i++)
  {
    u[i] = QVector();
    u[i][i] = 1.;
  }
  QVector e(1.);

  // test basic cube
  QOBBox obbox1(pt1, u, e);

  std::vector<QPoint> v = obbox1.vertices();
  QPoint verts[8];
  for(int i = 0; i < 8; i++)
  {
    verts[i] = v[i];
  }

  QOBBox obbox2 = primal::compute_oriented_bounding_box<CoordType, DIM>(verts, 8);

  EXPECT_TRUE(obbox2.contains(obbox1));

  // now test line of points
  QPoint pt2(1.);
  QPoint pt3(2.);
  verts[0] = pt1;
  verts[1] = pt2;
  verts[2] = pt3;

  QOBBox obbox3 = primal::compute_oriented_bounding_box<CoordType, DIM>(verts, 3);

  EXPECT_TRUE(obbox3.contains(pt1));
  EXPECT_TRUE(obbox3.contains(pt2));
  EXPECT_TRUE(obbox3.contains(pt3));

  EXPECT_TRUE(obbox3.getExtents()[1] < 1E-8);
  EXPECT_TRUE(obbox3.getExtents()[2] < 1E-8);
}

TEST(primal_compute_bounding_box, merge_oriented_box_test)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  const double ONE_OVER_SQRT_TWO = 0.7071;
  QPoint pt1;      // origin
  QVector u[DIM];  // make axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);
  QOBBox obbox1(pt1, u, e);

  // test when they're identical
  EXPECT_TRUE((obbox1 == primal::merge_boxes<CoordType, DIM>(obbox1, obbox1)));

  // test when they're completely separate
  QPoint pt2(10.);

  QOBBox obbox2(pt2, u, e);
  QOBBox obbox3 = primal::merge_boxes<CoordType, DIM>(obbox1, obbox2);
  EXPECT_TRUE(obbox3.contains(obbox1));
  EXPECT_TRUE(obbox3.contains(obbox2));

  // test when they only partially intersect
  QPoint pt3(0.5);
  QOBBox obbox4(pt3, u, e);
  QOBBox obbox5 = primal::merge_boxes<CoordType, DIM>(obbox1, obbox4);
  EXPECT_TRUE(obbox5.contains(obbox1));
  EXPECT_TRUE(obbox5.contains(obbox4));
}

TEST(primal_compute_bounding_box, merge_aligned_box_test)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  QPoint pt1(0.);
  QPoint pt2(1.);
  QBBox bbox1(pt1);
  QBBox bbox2(pt2);
  QBBox bbox3(pt1, pt2);

  EXPECT_TRUE((bbox3 == primal::merge_boxes<CoordType, DIM>(bbox1, bbox2)));

  QBBox bbox4(bbox3);
  QVector s(10.);
  bbox3.shift(s);

  QBBox bbox5 = primal::merge_boxes<CoordType, DIM>(bbox3, bbox4);

  EXPECT_TRUE(bbox5.contains(bbox3));
  EXPECT_TRUE(bbox5.contains(bbox4));
}

TEST(primal_compute_bounding_box, compute_quad_2d_box_test)
{
  constexpr int DIM = 2;
  using CoordinateType = double;
  using PointType = primal::Point<CoordinateType, DIM>;
  using BoundingBoxType = primal::BoundingBox<CoordinateType, DIM>;
  using QuadrilateralType = primal::Quadrilateral<CoordinateType, DIM>;

  PointType A {-1.0, 0.1};
  PointType B {-0.1, 0.5};
  PointType C {2.0, 0.0};
  PointType D {0.0, 1.0};

  QuadrilateralType quad {A, B, C, D};
  BoundingBoxType box = primal::compute_bounding_box<CoordinateType, DIM>(quad);

  EXPECT_TRUE(box.contains(A));
  EXPECT_TRUE(box.contains(B));
  EXPECT_TRUE(box.contains(C));
  EXPECT_TRUE(box.contains(D));
  EXPECT_EQ(box.getMin(), (PointType {-1.0, 0.0}));
  EXPECT_EQ(box.getMax(), (PointType {2.0, 1.0}));
}

TEST(primal_compute_bounding_box, compute_quad_3d_box_test)
{
  constexpr int DIM = 3;
  using CoordinateType = double;
  using PointType = primal::Point<CoordinateType, DIM>;
  using BoundingBoxType = primal::BoundingBox<CoordinateType, DIM>;
  using QuadrilateralType = primal::Quadrilateral<CoordinateType, DIM>;

  PointType A {-1.0, 0.1, 0.0};
  PointType B {-0.1, 0.5, 0.0};
  PointType C {2.0, 0.0, 1.0};
  PointType D {0.0, 1.0, 1.0};

  QuadrilateralType quad {A, B, C, D};
  BoundingBoxType box = primal::compute_bounding_box<CoordinateType, DIM>(quad);

  EXPECT_TRUE(box.contains(A));
  EXPECT_TRUE(box.contains(B));
  EXPECT_TRUE(box.contains(C));
  EXPECT_TRUE(box.contains(D));
  EXPECT_EQ(box.getMin(), (PointType {-1.0, 0.0, 0.0}));
  EXPECT_EQ(box.getMax(), (PointType {2.0, 1.0, 1.0}));
}

TEST(primal_compute_bounding_box, compute_oct_box_test)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;
  using QOct = primal::Octahedron<CoordType, DIM>;

  QPoint p1({1, 0, 0});
  QPoint p2({1, 1, 0});
  QPoint p3({0, 1, 0});
  QPoint p4({0, 1, 1});
  QPoint p5({0, 0, 1});
  QPoint p6({1, 0, 1});

  QOct oct(p1, p2, p3, p4, p5, p6);

  QBBox box = primal::compute_bounding_box<CoordType, DIM>(oct);
  EXPECT_TRUE(box.contains(p1));
  EXPECT_TRUE(box.contains(p2));
  EXPECT_TRUE(box.contains(p3));
  EXPECT_TRUE(box.contains(p4));
  EXPECT_TRUE(box.contains(p5));
  EXPECT_TRUE(box.contains(p6));
  EXPECT_EQ(box.getMin(), QPoint::zero());
  EXPECT_EQ(box.getMax(), QPoint::ones());
}

TEST(primal_compute_bounding_box, compute_poly_box_test)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;
  using QPoly = primal::Polyhedron<CoordType, DIM>;

  QPoly poly;
  QPoint p1({0, 0, 0});
  QPoint p2({1, 0, 0});
  QPoint p3({1, 1, 0});
  QPoint p4({0, 1, 0});
  QPoint p5({0, 0, 1});
  QPoint p6({1, 0, 1});
  QPoint p7({1, 1, 1});
  QPoint p8({0, 1, 1});

  poly.addVertex(p1);
  poly.addVertex(p2);
  poly.addVertex(p3);
  poly.addVertex(p4);
  poly.addVertex(p5);
  poly.addVertex(p6);
  poly.addVertex(p7);
  poly.addVertex(p8);

  QBBox box = primal::compute_bounding_box<CoordType, DIM>(poly);
  EXPECT_TRUE(box.contains(p1));
  EXPECT_TRUE(box.contains(p2));
  EXPECT_TRUE(box.contains(p3));
  EXPECT_TRUE(box.contains(p4));
  EXPECT_TRUE(box.contains(p5));
  EXPECT_TRUE(box.contains(p6));
  EXPECT_TRUE(box.contains(p7));
  EXPECT_TRUE(box.contains(p8));
  EXPECT_EQ(box.getMin(), QPoint::zero());
  EXPECT_EQ(box.getMax(), QPoint::ones());
}

TEST(primal_compute_bounding_box, compute_polygon_2)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;
  using QPoly = primal::Polygon<CoordType, DIM>;

  // empty polygon generates an empty (invalid) bounding box
  {
    QPoly empty_poly;
    QBBox box = primal::compute_bounding_box(empty_poly);
    EXPECT_FALSE(box.isValid());
  }

  // non-empty polygon generates expected axis-aligned bounding box
  {
    QPoint p1({1, 0});
    QPoint p2({5, 0});
    QPoint p3({7, 2});
    QPoint p4({9, 4});
    QPoint p5({-1, 4});

    QPoly polygon(5);
    polygon.addVertex(p1);
    polygon.addVertex(p2);
    polygon.addVertex(p3);
    polygon.addVertex(p4);
    polygon.addVertex(p5);

    QBBox box = primal::compute_bounding_box(polygon);
    EXPECT_TRUE(box.contains(p1));
    EXPECT_TRUE(box.contains(p2));
    EXPECT_TRUE(box.contains(p3));
    EXPECT_TRUE(box.contains(p4));
    EXPECT_TRUE(box.contains(p5));

    EXPECT_EQ(box.getMin(), (QPoint {-1., 0.}));
    EXPECT_EQ(box.getMax(), (QPoint {9., 4.}));
  }
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

#ifdef COMPUTE_BOUNDING_BOX_TESTER_SHOULD_SEED
  std::srand(std::time(0));
#else
  std::srand(105);
#endif

  int result = RUN_ALL_TESTS();
  return result;
}
