// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file primal_curved_polygon_intersect.cpp
 * \brief This file tests intersections of CurvedPolygon instances
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include <numeric>

namespace primal = axom::primal;

/*!
 * Helper function to compute the set of intersection polygons given two input polygons 
 * and to check that they match expectations, stored in \a expbPolygon. 
 * Intersection polygon is computed to within tolerance \a eps and checks use \a test_eps.
 */
template <typename CoordType, int DIM>
void checkIntersection(
  const primal::CurvedPolygon<CoordType, DIM>& bPolygon1,
  const primal::CurvedPolygon<CoordType, DIM>& bPolygon2,
  const std::vector<primal::CurvedPolygon<CoordType, DIM>> expbPolygon,
  const double eps = 1e-15,
  const double test_eps = 1e-13)
{
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;

  std::vector<CurvedPolygonType> intersectionPolys;

  //Compute intersection using algorithm with tolerance of eps
  intersect(bPolygon1, bPolygon2, intersectionPolys, eps);
  //Check that expected number of intersection regions are found
  EXPECT_EQ(expbPolygon.size(), intersectionPolys.size());

  //Check that expected intersection curves are found to within test_eps
  for(int i = 0; i < DIM; ++i)
  {
    int sz = intersectionPolys.size();
    for(int idxcurve = 0; idxcurve < sz; ++idxcurve)
    {
      int nEd = intersectionPolys[idxcurve].numEdges();
      for(int k = 0; k < nEd; ++k)
      {
        int ord = intersectionPolys[idxcurve][k].getOrder();
        for(int j = 0; j <= ord; ++j)
        {
          EXPECT_NEAR(expbPolygon[idxcurve][k][j][i],
                      intersectionPolys[idxcurve][k][j][i],
                      test_eps);
        }
      }
    }
  }
}

/*!
 * Helper function to create a CurvedPolygon from a list of control points and a list
 * of orders of component curves. Control points should be given as a list of Points
 * in order of orientation with no duplicates except that the first control point 
 * should also be the last control point (if the polygon is closed).  Orders should
 * be given as a list of ints in order of orientation, representing the orders of the component curves.
 */
template <typename CoordType, int DIM>
primal::CurvedPolygon<CoordType, DIM> createPolygon(
  const std::vector<primal::Point<CoordType, DIM>> ControlPoints,
  const std::vector<int> orders)
{
  using PointType = primal::Point<CoordType, DIM>;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int num_edges = orders.size();
  const int num_unique_control_points = ControlPoints.size();

  //checks if the orders and control points given will yield a valid curved polygon
  {
    const int sum_of_orders = std::accumulate(orders.begin(), orders.end(), 0);
    EXPECT_EQ(sum_of_orders + 1, num_unique_control_points);
  }

  //Converts the control points to BezierCurves of specified orders and stores them in a CurvedPolygon object.
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
  return bPolygon;
}

//----------------------------------------------------------------------------------

TEST(primal_curvedpolygon, intersection_triangle_linear)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of two linear triangles (single region)");

  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.3, 2.0},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};
  std::vector<int> orders = {1, 1, 1};

  CurvedPolygonType bPolygon1 = createPolygon(CP, orders);

  std::vector<PointType> CP2 = {PointType {0.71, 1.31},
                                PointType {0.41, 2.11},
                                PointType {0.11, 1.71},
                                PointType {0.71, 1.31}};

  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);

  std::vector<PointType> expCP = {
    PointType {0.3091666666666666666666, 1.9755555555555555555},
    PointType {0.11, 1.71},
    PointType {0.5083333333333333333, 1.44444444444444444444},
    PointType {0.3091666666666666666666, 1.9755555555555555555}};
  std::vector<int> exporders = {1, 1, 1};
  CurvedPolygonType expbPolygon = createPolygon(expCP, exporders);

  std::vector<CurvedPolygonType> expbPolygons {expbPolygon};
  checkIntersection(bPolygon1, bPolygon2, expbPolygons);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, intersection_triangle_quadratic)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersecting two quadratic triangles (single region)");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.4, 1.3},
                               PointType {0.3, 2.0},
                               PointType {0.27, 1.5},
                               PointType {0.0, 1.6},
                               PointType {0.1, 1.5},
                               PointType {0.6, 1.2}};
  std::vector<int> orders = {2, 2, 2};
  CurvedPolygonType bPolygon1 = createPolygon(CP, orders);

  std::vector<PointType> CP2 = {PointType {0.71, 1.31},
                                PointType {0.51, 1.41},
                                PointType {0.41, 2.11},
                                PointType {0.38, 1.61},
                                PointType {0.11, 1.71},
                                PointType {0.21, 1.61},
                                PointType {0.71, 1.31}};
  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);

  std::vector<PointType> expCP = {
    PointType {0.335956890729522, 1.784126953773395},
    PointType {0.297344765794753, 1.718171485335525},
    PointType {0.2395677533016981, 1.700128235793371},
    PointType {0.221884203146682, 1.662410644580941},
    PointType {0.199328465398189, 1.636873522352205},
    PointType {0.277429214338182, 1.579562422716502},
    PointType {0.408882616650578, 1.495574996394597},
    PointType {0.368520120719339, 1.616453177259694},
    PointType {0.335956890729522, 1.784126953773394}};
  std::vector<int> exporders = {2, 2, 2, 2};
  CurvedPolygonType expbPolygon = createPolygon(expCP, exporders);
  std::vector<CurvedPolygonType> expbPolygons = {expbPolygon};

  checkIntersection(bPolygon1, bPolygon2, expbPolygons, 1e-15, 1e-13);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, intersection_triangle_quadratic_two_regions)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of two quadratic triangles (two regions)");

  std::vector<PointType> CP1 = {PointType {0.6, 1.2},
                                PointType {0.4, 1.3},
                                PointType {0.3, 2.0},
                                PointType {0.27, 1.5},
                                PointType {0.0, 1.6},
                                PointType {0.1, 1.5},
                                PointType {0.6, 1.2}};

  std::vector<PointType> CP2 = {PointType {1.0205, 1.6699},
                                PointType {0.8339, 1.5467},
                                PointType {0.1777, 1.8101},
                                PointType {0.5957, 1.5341},
                                PointType {0.3741, 1.3503},
                                PointType {0.5107, 1.3869},
                                PointType {1.0205, 1.6699}};
  std::vector<int> orders = {2, 2, 2};

  std::vector<PointType> expCP1 = {
    PointType {0.343364196589264, 1.747080669655736},
    PointType {0.305984025190458, 1.760433098612141},
    PointType {0.266743999290327, 1.775316659915674},
    PointType {0.263419346128088, 1.763343410502168},
    PointType {0.259796003065908, 1.752116885838515},
    PointType {0.320641367919239, 1.705796408318085},
    PointType {0.362111919147859, 1.662268860466508},
    PointType {0.352450139541348, 1.702947255097842},
    PointType {0.343364196589264, 1.747080669655736}};

  std::vector<PointType> expCP2 = {
    PointType {0.454478985809487, 1.379250566393211},
    PointType {0.444689566319939, 1.400290430035245},
    PointType {0.435276730907216, 1.423589798138227},
    PointType {0.416268597450954, 1.385275578571685},
    PointType {0.374100000000000, 1.350300000000000},
    PointType {0.404839872482010, 1.358536305511285},
    PointType {0.454478985809487, 1.379250566393211}};

  std::vector<int> exporder1 = {2, 2, 2, 2};
  std::vector<int> exporder2 = {2, 2, 2};
  CurvedPolygonType bPolygon1 = createPolygon(CP1, orders);
  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);
  CurvedPolygonType expbPolygon1 = createPolygon(expCP1, exporder1);
  CurvedPolygonType expbPolygon2 = createPolygon(expCP2, exporder2);
  std::vector<CurvedPolygonType> expIntersections = {expbPolygon1, expbPolygon2};

  checkIntersection(bPolygon1, bPolygon2, expIntersections);
}

TEST(primal_curvedpolygon, area_intersection_triangle_inclusion)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of two quadratic triangles (inclusion)");

  std::vector<PointType> CP1 = {PointType {0.0, 0.0},
                                PointType {0.5, 0.0},
                                PointType {1.0, 0.0},
                                PointType {0.5, 0.5},
                                PointType {0.0, 1.0},
                                PointType {0.0, 0.5},
                                PointType {0.0, 0.0}};

  std::vector<PointType> CP2 = {PointType {0.05, 0.05},
                                PointType {0.30, 0.05},
                                PointType {0.55, 0.05},
                                PointType {0.30, 0.30},
                                PointType {0.05, 0.55},
                                PointType {0.05, 0.30},
                                PointType {0.05, 0.05}};
  std::vector<int> orders = {2, 2, 2};
  CurvedPolygonType bPolygon1 = createPolygon(CP1, orders);
  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);
  std::vector<CurvedPolygonType> expIntersections = {bPolygon2};

  checkIntersection(bPolygon1, bPolygon2, expIntersections);
  checkIntersection(bPolygon2, bPolygon1, expIntersections);
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
