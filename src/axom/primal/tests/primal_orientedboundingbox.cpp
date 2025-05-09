// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <limits>
#include <algorithm>

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

namespace primal = axom::primal;

TEST(primal_OBBox, obb_default_constructor)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  QOBBox bbox;
  EXPECT_FALSE(bbox.isValid()) << "Default constructed bounding box is invalid";
  EXPECT_FALSE(bbox.contains(QPoint()))
    << "Default constructed bounding box should not contain any points";
  EXPECT_FALSE(bbox.contains(QPoint(1000)))
    << "Default constructed bounding box should not contain any points";
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_ctor_from_singlePt)
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
  QVector e;      // zero extents
  QPoint pt2(2);  // (2,2,2)

  QOBBox obbox1(pt1);

  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox1.contains(pt1));
  EXPECT_FALSE(obbox1.contains(pt2));
  EXPECT_TRUE(obbox1.getCentroid() == pt1)
    << "OBBox only has a single point, so bbox1.getCentroid()==pt1";
  const QVector* u_in = obbox1.getAxes();

  for(int i = 0; i < DIM; i++)
  {
    EXPECT_TRUE(u_in[i] == u[i]);
  }

  EXPECT_TRUE(obbox1.getExtents() == e) << "Extents should be 0";

  QOBBox obbox2(pt2);

  EXPECT_TRUE(obbox2.isValid());
  EXPECT_TRUE(obbox2.contains(pt2));
  EXPECT_FALSE(obbox2.contains(pt1));
  EXPECT_TRUE(obbox2.getCentroid() == pt2)
    << "OBBox only has a single point, so obbox2.getCentroid()==pt2";
  u_in = obbox2.getAxes();
  for(int i = 0; i < DIM; i++)
  {
    EXPECT_TRUE(u_in[i] == u[i]);
  }
  EXPECT_TRUE(e == obbox2.getExtents()) << "Extents should be 0";
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_ctor_from_data)
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
  QVector e(1.);    //extents
  QPoint pt2(0.5);  // (.5,.5,.5)
  QPoint pt3(2);    // (2,2,2)

  QOBBox obbox1(pt1, u, e);

  // check containments
  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox1.contains(pt1));
  EXPECT_TRUE(obbox1.contains(pt2));
  EXPECT_FALSE(obbox1.contains(pt3));

  // check settings
  const QVector* u_in = obbox1.getAxes();
  for(int i = 0; i < DIM; i++)
  {
    EXPECT_TRUE(u_in[i] == u[i]);
  }
  EXPECT_TRUE(obbox1.getCentroid() == pt1);
  EXPECT_TRUE(obbox1.getExtents() == e);
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_ctor_from_point_array)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  QPoint pt1;  // origin
  QPoint pt2({1.0, 0.0, 0.0});
  QPoint pt3({0.0, 1.0, 0.0});
  QPoint pt4({0.0, 0.0, 1.0});
  QPoint pt5({1.0, 1.0, 0.0});
  QPoint pt6({1.0, 0.0, 1.0});
  QPoint pt7({0.0, 1.0, 1.0});
  QPoint pt8({1.0, 1.0, 1.0});

  /* -1D OBB */
  QPoint* pts_00 = nullptr;
  QOBBox obbox00(pts_00, 0);

  EXPECT_FALSE(obbox00.isValid());

  /* 0D OBB */
  QPoint pts_0d[] = {pt1};
  QOBBox obbox0(pts_0d, 1);

  EXPECT_TRUE(obbox0.isValid());
  EXPECT_TRUE(obbox0.contains(pt1));

  EXPECT_NEAR(obbox0.getCentroid()[0], 0.0, 1e-6);
  EXPECT_NEAR(obbox0.getCentroid()[1], 0.0, 1e-6);
  EXPECT_NEAR(obbox0.getCentroid()[2], 0.0, 1e-6);

  /* 1D OBB */
  QPoint pts_1d[] = {pt1, pt2};
  QOBBox obbox1(pts_1d, 2);

  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox1.contains(pt1));
  EXPECT_TRUE(obbox1.contains(pt2));

  EXPECT_NEAR(obbox1.getCentroid()[0], 0.5, 1e-6);
  EXPECT_NEAR(obbox1.getCentroid()[1], 0.0, 1e-6);
  EXPECT_NEAR(obbox1.getCentroid()[2], 0.0, 1e-6);

  /* 2D OBB */
  QPoint pts_2d[] = {pt1, pt2, pt3, pt5};
  QOBBox obbox2(pts_2d, 4);

  EXPECT_TRUE(obbox2.isValid());
  EXPECT_TRUE(obbox2.contains(pt1));
  EXPECT_TRUE(obbox2.contains(pt2));
  EXPECT_TRUE(obbox2.contains(pt3));
  EXPECT_TRUE(obbox2.contains(pt5));

  EXPECT_NEAR(obbox2.getCentroid()[0], 0.5, 1e-6);
  EXPECT_NEAR(obbox2.getCentroid()[1], 0.5, 1e-6);
  EXPECT_NEAR(obbox2.getCentroid()[2], 0.0, 1e-6);

  /* 3D OBB */
  QPoint pts_3d[] = {pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8};
  QOBBox obbox3(pts_3d, 8);

  // check containments
  EXPECT_TRUE(obbox3.isValid());
  for(int i = 0; i < 8; i++)
  {
    EXPECT_TRUE(obbox3.contains(pts_3d[i]));
  }

  // check settings
  EXPECT_TRUE(obbox3.getCentroid() == QPoint(0.5));
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_test_clear)
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
  QVector e(1.);    //extents
  QPoint pt2(0.5);  // (.5,.5,.5)
  QPoint pt3(2);    // (2,2,2)

  QOBBox obbox1(pt1, u, e);

  // clear
  obbox1.clear();

  EXPECT_FALSE(obbox1.isValid());
  EXPECT_FALSE(obbox1.contains(pt1));
  EXPECT_FALSE(obbox1.contains(pt2));
  EXPECT_FALSE(obbox1.contains(pt3));
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_test_vertices)
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
  QVector e(1.);  //extents

  QOBBox obbox1(pt1, u, e);
  std::vector<QPoint> l = obbox1.vertices();

  axom::NumericArray<CoordType, DIM> v;

  for(int i = 0; i < 2; i++)
  {
    v[0] = 1. - 2. * i;
    for(int j = 0; j < 2; j++)
    {
      v[1] = 1. - 2. * j;
      for(int k = 0; k < 2; k++)
      {
        v[2] = 1. - 2. * k;
        EXPECT_TRUE(std::find(l.begin(), l.end(), QPoint(v)) != l.end());
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_test_add_point)
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
  QVector e(1.);  //extents

  QOBBox obbox1(pt1, u, e);
  QOBBox obbox2 = obbox1;
  // add point already inside, vertiy nothing changes
  obbox2.addPoint(pt1);

  EXPECT_EQ(obbox1, obbox2);

  // add point
  QPoint pt2(10.);
  obbox2.addPoint(pt2);
  EXPECT_TRUE(obbox2.contains(pt2));

  QVector e2 = obbox2.getExtents();

  for(int i = 0; i < 3; i++)
  {
    EXPECT_EQ(e2[i], 10.);
  }

  // test when it's shifted away from the origin
  obbox2.shift(QVector(10.));
  obbox2.addPoint(QPoint(-10.));

  QOBBox obbox4 = QOBBox(pt2, u, QVector(20.));
  EXPECT_EQ(obbox4, obbox2);

  // test adding a vertex does nothing
  QOBBox obbox3 = obbox2;
  obbox2.addPoint(obbox2.vertices()[0]);
  EXPECT_EQ(obbox2, obbox3);
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_test_add_box)
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
  QVector e(1.);  //extents

  QOBBox obbox1(pt1, u, e);
  QOBBox obbox2 = obbox1;

  // test addding obbox1 to itself does nothing
  obbox1.addBox(obbox1);
  EXPECT_EQ(obbox1, obbox2);

  // bisect and then add back together; should contain original
  QOBBox obbox3;
  QOBBox obbox4;
  obbox1.bisect(obbox3, obbox4);
  obbox3.addBox(obbox4);
  EXPECT_TRUE(obbox3.contains(obbox1));

  // adding empty/invalid box is a no-op
  QOBBox obbox5;
  obbox3.addBox(obbox5);
  EXPECT_FALSE(obbox5.isValid());
  EXPECT_TRUE(obbox3.isValid());
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_test_expand)
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
  QVector e(1.);    //extents
  QPoint pt2(0.5);  // (.5,.5,.5)
  QPoint pt3(2);    // (2,2,2)

  QOBBox obbox1(pt1, u, e);

  (obbox1.expand(0.5)).expand(-0.5);
  EXPECT_TRUE(pt1 == obbox1.getCentroid());
  EXPECT_TRUE(e == obbox1.getExtents());
  const QVector* u_in = obbox1.getAxes();

  for(int i = 0; i < DIM; i++)
  {
    EXPECT_TRUE(u_in[i] == u[i]);
  }
  obbox1.expand(10.);
  QVector newE = QVector(11.);
  EXPECT_TRUE(pt1 == obbox1.getCentroid());
  EXPECT_TRUE(newE == obbox1.getExtents());

  for(int i = 0; i < DIM; i++)
  {
    EXPECT_TRUE(u_in[i] == u[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_test_scale)
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
  QVector e(1.);    //extents
  QPoint pt2(0.5);  // (.5,.5,.5)
  QPoint pt3(2);    // (2,2,2)

  QOBBox obbox1(pt1, u, e);
  QOBBox obbox2(pt1, u, e);
  QOBBox obbox3(pt1);

  obbox1.scale(1.);
  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox1 == obbox2);
  EXPECT_TRUE(obbox1.contains(pt2));

  obbox1.scale(-1.);
  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox1 == obbox2);
  EXPECT_TRUE(obbox1.contains(pt2));

  obbox1.scale(0.);
  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox1 == obbox3);
  EXPECT_TRUE(obbox1.contains(pt1));
  EXPECT_FALSE(obbox1.contains(pt2));

  obbox1.scale(10.);
  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox1 == obbox3);
  EXPECT_TRUE(obbox1.contains(pt1));

  obbox2.scale(5.);
  EXPECT_TRUE(obbox2.isValid());
  EXPECT_TRUE(obbox2.contains(pt3));
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_test_shift)
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
  QVector e(1.);    //extents
  QPoint pt2(0.5);  // (.5,.5,.5)
  QPoint pt3(2);    // (2,2,2)
  QVector disp = QVector(2.);

  QOBBox obbox1(pt1, u, e);
  QOBBox obbox2(pt3, u, e);

  obbox1.shift(disp);
  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox1 == obbox2);
  EXPECT_FALSE(obbox1.contains(pt1));
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_copy_and_assignment)
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

  QVector e(1.);    //extents
  QPoint pt2(0.5);  // (.5,.5,.5)
  QPoint pt3(2);    // (2,2,2)

  QOBBox obbox1(pt1, u, e);
  QOBBox obbox2 = obbox1;
  QOBBox obbox3(obbox1);

  EXPECT_TRUE(obbox2.isValid());
  EXPECT_TRUE(obbox3.isValid());

  EXPECT_TRUE(obbox1 == obbox2);
  EXPECT_TRUE(obbox1 == obbox3);
  EXPECT_TRUE(obbox2 == obbox3);

  const QVector* u_in = obbox2.getAxes();
  for(int i = 0; i < DIM; i++)
  {
    EXPECT_TRUE(u_in[i] == u[i]);
  }
  EXPECT_TRUE(e == obbox2.getExtents());
  EXPECT_TRUE(pt1 == obbox2.getCentroid());

  // test out infinitesimal box
  obbox2 = QOBBox(pt1);
  EXPECT_TRUE(obbox2.contains(pt1));
  EXPECT_FALSE(obbox2.contains(pt2));
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_contains_obb)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  QPoint pt1;      // origin
  QVector u[DIM];  // make standard axes
  QVector u_o[DIM];
  for(int i = 0; i < DIM; i++)
  {
    u[i] = QVector();
    u[i][i] = 1.;
    u_o[i] = QVector();
  }
  u_o[0][0] = 1.;
  u_o[0][2] = -1.;
  u_o[1][0] = 1.;
  u_o[1][2] = 1.;
  u_o[2][1] = 1.;

  QVector e(1.);
  QVector e_o(.0001);

  QOBBox obbox1(pt1, u, e);
  QOBBox obbox2(pt1, u_o, e_o);

  EXPECT_TRUE(obbox1.isValid());
  EXPECT_TRUE(obbox2.isValid());
  EXPECT_TRUE(obbox1.contains(obbox2));
  EXPECT_FALSE(obbox2.contains(obbox1));

  QPoint pt3(2.);  // (2,2,2)
  QOBBox obbox3(pt3);
  EXPECT_FALSE(obbox2.contains(obbox3));
  EXPECT_FALSE(obbox1.contains(obbox3));
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_to_local)
{
  constexpr int DIM = 3;
  constexpr double EPS = 1e-8;

  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  // standard axes
  QVector u[DIM] = {QVector {1, 0, 0}, QVector {0, 1, 0}, QVector {0, 0, 1}};

  // another othogonal basis
  QVector u_o[DIM] = {QVector {1, 0, -1}, QVector {1, 0, 1}, QVector {0, 1, 0}};

  QPoint pt1;
  QVector vec0(0.);

  QPoint pt2(10.);
  QVector vec1(10.);

  QVector e(1.);

  QOBBox obbox1(pt1, u, e);

  // if box is standard one centered at origin, nothing should change
  EXPECT_EQ(vec0, QVector(obbox1.toLocal(pt1)));
  EXPECT_EQ(vec1, QVector(obbox1.toLocal(pt2)));

  obbox1.shift(-vec1);

  // now pt2's local coords are (20., 20., 20.)
  EXPECT_EQ(2. * vec1, QVector(obbox1.toLocal(pt2)));

  QOBBox obbox2(pt1, u_o, e);

  // local coordinates of getCentroid should be 0
  EXPECT_EQ(vec0, QVector(obbox2.toLocal(obbox2.getCentroid())));

  // can roughly compute local coords of pt2
  QVector vec2(obbox2.toLocal(pt2));
  EXPECT_NEAR(vec2[0], 0., EPS);
  EXPECT_TRUE(((vec2[1] - 14.142) < 0.1) && ((14.142 - vec2[1]) < 0.1));
  EXPECT_TRUE(((vec2[2] - 10.) < 0.1) && ((10. - vec2[2]) < 0.1));
}

//------------------------------------------------------------------------------
TEST(primal_OBBox, obb_bisect)
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

  QVector e = QVector();
  e[0] = 2.;
  e[1] = 1.;
  e[2] = 1.;

  QOBBox obbox1(pt1, u, e);
  QOBBox obbox2;
  QOBBox obbox3;

  obbox1.bisect(obbox2, obbox3);

  QVector shift;  // zeros
  shift[0] = 2.;
  obbox3.shift(shift);
  EXPECT_TRUE(obbox2 == obbox3);
}

TEST(primal_OBBox, obb_test_furthest_point)
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
  QPoint pt2(-1.);
  QPoint pt3(2.);

  QOBBox obbox1(pt1, u, e);

  EXPECT_TRUE(obbox1.furthestPoint(pt3) == pt2);

  QPoint pt4(-0.1);
  QPoint pt5(1.);
  EXPECT_TRUE(obbox1.furthestPoint(pt4) == pt5);

  pt4.array()[1] = 0.1;
  pt5.array()[1] = -1.;
  EXPECT_TRUE(obbox1.furthestPoint(pt4) == pt5);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

#ifdef ORIENTEDBOUNDINGBOX_TESTER_SHOULD_SEED
  std::srand(std::time(0));
#else
  std::srand(105);
#endif

  int result = RUN_ALL_TESTS();
  return result;
}
