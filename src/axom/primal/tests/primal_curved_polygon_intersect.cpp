// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file primal_curved_polygon_intersect.cpp
 * /brief This file tests intersections of CurvedPolygon instances
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/primal.hpp"

namespace primal = axom::primal;

/**
 * Helper function to compute the area and centroid of a curved polygon and to check that they match expectations, stored in \a expArea and \a expCentroid. Areas and Moments are computed within tolerance \a eps and checks use \a test_eps.
 */
template <typename CoordType, int DIM>
void checkMoments(const primal::CurvedPolygon<CoordType, DIM>& bPolygon,
                  const CoordType expArea,
                  const primal::Point<CoordType, DIM>& expMoment,
                  double eps,
                  double test_eps)
{
  EXPECT_NEAR(expArea, primal::area(bPolygon, eps), test_eps);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(expMoment[i], primal::centroid(bPolygon, eps)[i], test_eps);
  }
}

template <typename CoordType, int DIM>
primal::CurvedPolygon<CoordType, DIM> createPolygon(
  std::vector<primal::Point<CoordType, DIM>> ControlPoints,
  std::vector<int> orders)
{
  using PointType = primal::Point<CoordType, DIM>;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int num_edges = orders.size();
  const int num_unique_control_points = ControlPoints.size();

  //std::cout << num_edges << ", " << num_unique_control_points << std::endl;
  //std::cout << ControlPoints << std::endl;
  for(int i = 0; i < num_edges; ++i)
  {
    std::cout << orders[i] << std::endl;
  }

  //checks if the orders and control points given will give a valid polygon
  EXPECT_EQ(accumulate(orders.begin(), orders.end(), 0) + 1,
            num_unique_control_points);

  CurvedPolygonType bPolygon;
  int iter = 0;
  for(int j = 0; j < num_edges; ++j)
  {
    std::vector<PointType> subCP;
    subCP.assign(ControlPoints.begin() + iter,
                 ControlPoints.begin() + iter + orders[j] + 1);
    BezierCurveType addCurve(subCP, orders[j]);
    bPolygon.addEdge(addCurve);
    iter += (orders[j]);
  }
  std::cout << bPolygon << std::endl;
  return bPolygon;
}

TEST(primal_curvedpolygon, intersection_triangle_linear)
{
  const int DIM = 2;
  const int order = 1;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two linear triangular CurvedPolygons (single region).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);
  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] += .11;
      }
    }
  }

  PointType expcontrolPoints[order + 1] = {
    PointType::make_point(0.3091666666666666666666, 1.9755555555555555555),
    PointType::make_point(0.11, 1.71)};

  PointType expcontrolPoints2[order + 1] = {
    PointType::make_point(0.11, 1.71),
    PointType::make_point(0.5083333333333333333, 1.44444444444444444444)};

  PointType expcontrolPoints3[order + 1] = {
    PointType::make_point(0.5083333333333333333, 1.44444444444444444444),
    PointType::make_point(0.3091666666666666666666, 1.9755555555555555555)};

  CurvedPolygonType expbPolygon;

  BezierCurveType expbCurve(expcontrolPoints, order);
  expbPolygon.addEdge(expbCurve);
  BezierCurveType expbCurve2(expcontrolPoints2, order);
  expbPolygon.addEdge(expbCurve2);
  BezierCurveType expbCurve3(expcontrolPoints3, order);
  expbPolygon.addEdge(expbCurve3);

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);
  std::cout << bPolygons3.size() << std::endl;
  // int nEd= bPolygons3.size();
  for(int i = 0; i < static_cast<int>(bPolygons3.size()); ++i)
  {
    std::cout << bPolygons3[i] << std::endl;
  }
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int idxcurve = 0; idxcurve < static_cast<int>(bPolygons3.size());
          ++idxcurve)
      {
        for(int k = 0; k < bPolygons3[idxcurve].numEdges(); ++k)
        {
          EXPECT_NEAR(expbPolygon[k][j][i], bPolygons3[idxcurve][k][j][i], 1e-10);
          /*EXPECT_TRUE(axom::utilities::isNearlyEqual(expbPolygon[k][j][i],bPolygons3[idxcurve][k][j][i],2.));*/
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, intersection_triangle_quadratic)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two quadratic triangular CurvedPolygons (single "
    "region).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.4, 1.3),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.27, 1.5),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.1, 1.5),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);
  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] += .11;
      }
    }
  }

  PointType expcontrolPoints[order + 1] = {
    PointType::make_point(0.335956890729522, 1.784126953773395),
    PointType::make_point(0.297344765794753, 1.718171485335525),
    PointType::make_point(0.239567753301698, 1.700128235793372)};

  PointType expcontrolPoints2[order + 1] = {
    PointType::make_point(0.2395677533016981, 1.700128235793371),
    PointType::make_point(0.221884203146682, 1.662410644580941),
    PointType::make_point(0.199328465398189, 1.636873522352205)};

  PointType expcontrolPoints3[order + 1] = {
    PointType::make_point(0.199328465398188, 1.636873522352206),
    PointType::make_point(0.277429214338182, 1.579562422716502),
    PointType::make_point(0.408882616650578, 1.495574996394597)};

  PointType expcontrolPoints4[order + 1] = {
    PointType::make_point(0.408882616650588, 1.495574996394586),
    PointType::make_point(0.368520120719339, 1.616453177259694),
    PointType::make_point(0.335956890729522, 1.784126953773394)};

  CurvedPolygonType expbPolygon;

  BezierCurveType expbCurve(expcontrolPoints, order);
  expbPolygon.addEdge(expbCurve);
  BezierCurveType expbCurve2(expcontrolPoints2, order);
  expbPolygon.addEdge(expbCurve2);
  BezierCurveType expbCurve3(expcontrolPoints3, order);
  expbPolygon.addEdge(expbCurve3);
  BezierCurveType expbCurve4(expcontrolPoints4, order);
  expbPolygon.addEdge(expbCurve4);

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);

  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int idxcurve = 0; idxcurve < static_cast<int>(bPolygons3.size());
          ++idxcurve)
      {
        for(int k = 0; k < bPolygons3[idxcurve].numEdges(); ++k)
        {
          EXPECT_TRUE(axom::utilities::isNearlyEqual(expbPolygon[k][j][i],
                                                     bPolygons3[idxcurve][k][j][i],
                                                     2.));
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, area_intersection_triangle_linear)
{
  const int DIM = 2;
  const int order = 1;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test finding area of intersection two linear triangular CurvedPolygons "
    "(single region).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);
  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] += .11;
      }
    }
  }

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);

  CoordType A = 0.0;
  for(int i = 0; i < static_cast<int>(bPolygons3.size()); ++i)
  {
    A += primal::area(bPolygons3[i], 1e-14);
  }
  CoordType expA = -0.0793347222222222222;
  EXPECT_NEAR(A, expA, 1e-10);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, area_intersection_triangle_quadratic)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two quadratic triangular CurvedPolygons (single "
    "region).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.4, 1.3),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.27, 1.5),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.1, 1.5),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);
  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] += .11;
      }
    }
  }

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3, 1e-14);
  EXPECT_TRUE(didIntersect);

  CoordType A = 0.0;
  for(int i = 0; i < static_cast<int>(bPolygons3.size()); ++i)
  {
    A += primal::area(bPolygons3[i], 1e-8);
  }
  CoordType expA = -0.024649833203616;
  EXPECT_NEAR(A, expA, 1e-10);
}

TEST(primal_curvedpolygon, area_intersection_triangle_quadratic_two_regions)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two quadratic triangular CurvedPolygons (two regions).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.4, 1.3),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.27, 1.5),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.1, 1.5),
                                         PointType::make_point(0.6, 1.2)};

  PointType controlPoints4[order + 1] = {PointType::make_point(1.0205, 1.6699),
                                         PointType::make_point(0.8339, 1.5467),
                                         PointType::make_point(0.1777, 1.8101)};

  PointType controlPoints5[order + 1] = {PointType::make_point(0.1777, 1.8101),
                                         PointType::make_point(0.5957, 1.5341),
                                         PointType::make_point(0.3741, 1.3503)};

  PointType controlPoints6[order + 1] = {PointType::make_point(0.3741, 1.3503),
                                         PointType::make_point(0.5107, 1.3869),
                                         PointType::make_point(1.0205, 1.6699)};
  CurvedPolygonType bPolygon2;

  BezierCurveType bCurve(controlPoints, order);
  BezierCurveType bCurve4(controlPoints4, order);
  bPolygon2.addEdge(bCurve4);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  BezierCurveType bCurve5(controlPoints5, order);
  bPolygon2.addEdge(bCurve5);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  BezierCurveType bCurve6(controlPoints6, order);
  bPolygon2.addEdge(bCurve6);
  bPolygon.addEdge(bCurve3);

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);
  EXPECT_EQ(bPolygons3.size(), 2);
}

TEST(primal_curvedpolygon, area_intersection_triangle_inclusion)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two quadratic triangular CurvedPolygons (inclusion).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.0, 0.0),
                                        PointType::make_point(0.5, 0.0),
                                        PointType::make_point(1.0, 0.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(1.0, 0.0),
                                         PointType::make_point(0.5, 0.5),
                                         PointType::make_point(0.0, 1.0)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.0),
                                         PointType::make_point(0.0, 0.5),
                                         PointType::make_point(0.0, 0.0)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);

  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] = bPolygon2[k][j][i] * .5 + .05;
      }
    }
  }

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  bPolygons3.clear();
  bool didIntersect2 = intersect(bPolygon2, bPolygon, bPolygons3);
  EXPECT_TRUE(didIntersect);
  EXPECT_TRUE(didIntersect2);
  EXPECT_EQ(bPolygons3.size(), 1);
}
//----------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
