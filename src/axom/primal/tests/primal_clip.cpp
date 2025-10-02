// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/core/Types.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/numerics/transforms.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/construct.hpp"

#include "axom/primal/operators/clip.hpp"
#include "axom/primal/operators/intersection_volume.hpp"
#include "axom/primal/operators/split.hpp"
#include "axom/primal/operators/compute_bounding_box.hpp"

#include <limits>
#include <cmath>

#include "gtest/gtest.h"

namespace Primal3D
{
using PointType = axom::primal::Point<double, 3>;
using VectorType = axom::primal::Vector<double, 3>;
using BoundingBoxType = axom::primal::BoundingBox<double, 3>;
using HexahedronType = axom::primal::Hexahedron<double, 3>;
using TriangleType = axom::primal::Triangle<double, 3>;
using TetrahedronType = axom::primal::Tetrahedron<double, 3>;
using OctahedronType = axom::primal::Octahedron<double, 3>;
using PolyhedronType = axom::primal::Polyhedron<double, 3>;
using PolygonType = axom::primal::Polygon<double, 3>;
using PlaneType = axom::primal::Plane<double, 3>;
using PolyhedronType = axom::primal::Polyhedron<double, 3>;
}  // namespace Primal3D

TEST(primal_clip, simple_clip)
{
  using namespace Primal3D;

  // all checks in this test are against the cube [-1,-1,-1] to [1,1,1]
  BoundingBoxType bbox;
  bbox.addPoint(PointType {-1, -1, -1});
  bbox.addPoint(PointType {1, 1, 1});

  // define a slightly inflated bbox for containment tests
  constexpr double EPS = 1e-8;
  BoundingBoxType expanded_bbox(bbox);
  expanded_bbox.expand(EPS);

  std::vector<std::pair<TriangleType, int>> test_cases;

  // add a bunch of test cases
  {
    // this triangle is clearly outside the reference cube
    test_cases.push_back(
      std::make_pair(TriangleType {PointType {2, 2, 2}, PointType {2, 2, 4}, PointType {2, 4, 2}}, 0));

    // triangle at z=0 that spans entire cube
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {-100, -100, 0.}, PointType {-100, 100, 0.}, PointType {100, 0, 0.}},
      4));

    // triangle at z=0 w/ two vertices inside unit cube
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {0.25, 0.25, 0.}, PointType {0.75, 0.25, 0.}, PointType {1.5, 0.5, 0.}},
      4));

    // triangle at z=0 w/ vertices aligned w/ bounding box planes
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {2, 1, 0.}, PointType {2, 2, 0.}, PointType {1, 2, 0.}},
      0));

    // triangle at z=0 w/ one vertex coincident w/ cube
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {1, 2, 0.}, PointType {1, 1, 0.}, PointType {2, 1, 0.}},
      0));

    // triangle at z=0 w/ one edge midpoint coincident w/ cube
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {0, 2, 0.}, PointType {2, 0, 0.}, PointType {2, 2, 0.}},
      0));

    // triangle at z=0 w/ one edge coincident w/ unit cube and the other outside
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {-10, 1, 0.}, PointType {10, 1, 0.}, PointType {0, 10, 0.}},
      0));

    // triangle at z=0 w/ one edge coincident w/ unit cube and the last vertex inside
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {-10, 1, 0.}, PointType {10, 1, 0.}, PointType {0, 0, 0.}},
      5));

    // triangle at z=0 w/ one edge coincident w/ unit cube and the last vertex on the other side
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {-10, 1, 0.}, PointType {10, 1, 0.}, PointType {0, -10, 0.}},
      4));

    // triangle at z=1 w/ one edge coincident w/ unit cube and the last vertex on the other side
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {-10, 1, 1.}, PointType {10, 1, 1.}, PointType {0, -10, 1.}},
      4));

    // triangle at z=.5 w/ two vertices inside cube and one to the left
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {-2, .25, .5}, PointType {.25, .25, .5}, PointType {.25, .75, .5}},
      4));

    // triangle at z=.5 w/ two vertices inside cube and one to the left
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {-2, .25, .5}, PointType {.75, .25, .5}, PointType {.75, .75, .5}},
      4));

    // all vertices outside cube, one vertex above, one vertex to the left and one to the top left
    test_cases.push_back(std::make_pair(
      TriangleType {PointType {.9, 100, .5}, PointType {100, -1, .5}, PointType {100, 100, .5}},
      0));
  }

  // check each test case
  for(const auto& pair : test_cases)
  {
    const auto& tri = pair.first;
    const int expected_verts = pair.second;

    PolygonType poly = axom::primal::clip(tri, bbox);

    SLIC_INFO(
      axom::fmt::format("Intersection of triangle {} and bbox {} is polygon {}", tri, bbox, poly));

    // The clipped polygon
    // ... should have the expected number of vertices
    EXPECT_EQ(expected_verts, poly.numVertices());
    // ... should lie inside the bounding box (inflated by EPS to deal w/ boundaries)
    EXPECT_TRUE(expanded_bbox.contains(axom::primal::compute_bounding_box(poly)));
    // ... and its area should be at most that of the original triangle
    EXPECT_GE(tri.area(), poly.area());
  }
}

TEST(primal_clip, unit_simplex)
{
  using namespace Primal3D;
  constexpr double delta = 1e-5;

  // Test the "unit simplex", and a jittered version
  PointType points[] = {PointType {1, 0, 0},
                        PointType {0, 1, 0},
                        PointType {0, 0, 1},
                        PointType {1 + delta, delta, delta},
                        PointType {delta, 1 + delta, delta},
                        PointType {delta, delta, 1 + delta}};

  BoundingBoxType bbox;
  bbox.addPoint(PointType::zero());
  bbox.addPoint(PointType(.75));

  // intersection of this triangle and cube is a hexagon
  {
    TriangleType tri(points[0], points[1], points[2]);

    PolygonType poly = axom::primal::clip(tri, bbox);
    EXPECT_EQ(6, poly.numVertices());

    SLIC_INFO("Intersection of triangle " << tri << " and bounding box " << bbox << " is polygon"
                                          << poly);
  }
  {
    TriangleType tri(points[3], points[4], points[5]);

    PolygonType poly = axom::primal::clip(tri, bbox);
    EXPECT_EQ(6, poly.numVertices());

    SLIC_INFO("Intersection of triangle " << tri << " and bounding box " << bbox << " is polygon"
                                          << poly);
  }
}

TEST(primal_clip, boundingBoxOptimization)
{
  using namespace Primal3D;

  SLIC_INFO("Checking correctness of optimization for skipping clipping "
            << " of planes that the triangle's bounding box doesn't cover");

  constexpr double VAL1 = 3.;
  constexpr double VAL2 = 2.;

  BoundingBoxType bbox;
  bbox.addPoint(PointType(-1.));
  bbox.addPoint(PointType(1.));

  PointType midpoint = PointType::zero();

  PointType points[] = {
    PointType {VAL1, VAL2, 0},
    PointType {-VAL1, VAL2, 0},
    PointType {VAL1, -VAL2, 0},
    PointType {-VAL1, -VAL2, 0},

    PointType {VAL1, 0, VAL2},
    PointType {-VAL1, 0, VAL2},
    PointType {VAL1, 0, -VAL2},
    PointType {-VAL1, 0, -VAL2},

    PointType {0, VAL2, VAL1},
    PointType {0, VAL2, -VAL1},
    PointType {0, -VAL2, VAL1},
    PointType {0, -VAL2, -VAL1},

    PointType {0, VAL1, VAL2},
    PointType {0, -VAL1, VAL2},
    PointType {0, VAL1, -VAL2},
    PointType {0, -VAL1, -VAL2},
  };

  for(int i = 0; i < 16; i += 2)
  {
    TriangleType tri(midpoint, points[i], points[i + 1]);
    PolygonType poly = axom::primal::clip(tri, bbox);
    SLIC_INFO(poly);
    EXPECT_EQ(5, poly.numVertices());
  }
}

TEST(primal_clip, experimentalData)
{
  using namespace Primal3D;

  constexpr double EPS = 1e-8;

  // Triangle 248 from sphere mesh
  TriangleType tri(PointType {0.405431, 3.91921, 3.07821},
                   PointType {1.06511, 3.96325, 2.85626},
                   PointType {0.656002, 4.32465, 2.42221});

  // Block index {grid pt: (19,29,24); level: 5} from InOutOctree
  BoundingBoxType box12(PointType {0.937594, 4.06291, 2.50025},
                        PointType {1.25012, 4.37544, 2.81278});

  PolygonType poly = axom::primal::clip(tri, box12);
  EXPECT_EQ(3, poly.numVertices());

  SLIC_INFO("Intersection of triangle " << tri << " \n\t and bounding box " << box12
                                        << " \n\t is polygon" << poly << " with centroid "
                                        << poly.vertexMean());

  // Check that the polygon vertices are on the triangle
  for(int i = 0; i < poly.numVertices(); ++i)
  {
    PointType bary = tri.physToBarycentric(poly[i]);
    PointType reconstructed = tri.baryToPhysical(bary);

    SLIC_INFO("Testing clipped polygon point " << poly[i] << "-- w/ barycentric coords " << bary
                                               << "\n\t-- reconstructed point is: " << reconstructed
                                               << "...\n");

    double barySum = bary[0] + bary[1] + bary[2];
    EXPECT_NEAR(1., barySum, EPS);

    for(int dim = 0; dim < 3; ++dim)
    {
      EXPECT_GE(bary[dim], -EPS);
      EXPECT_NEAR(poly[i][dim], reconstructed[dim], EPS);
    }

    EXPECT_TRUE(box12.contains(poly[i]));
  }

  // Check that the polygon centroid is on the triangle
  {
    PointType centroid = poly.vertexMean();
    PointType bary = tri.physToBarycentric(centroid);
    PointType reconstructed = tri.baryToPhysical(bary);

    SLIC_INFO("Testing clipped polygon centroid "
              << centroid << "-- w/ barycentric coords " << bary
              << "\n\t-- reconstructed point is: " << reconstructed << "...\n");

    double barySum = bary[0] + bary[1] + bary[2];
    EXPECT_NEAR(1., barySum, EPS);

    for(int dim = 0; dim < 3; ++dim)
    {
      EXPECT_GE(bary[dim], -EPS);
      EXPECT_NEAR(centroid[dim], reconstructed[dim], EPS);
    }
    EXPECT_TRUE(box12.contains(centroid));
  }
}

template <typename ExecPolicy>
void unit_check_poly_clip()
{
  using namespace Primal3D;

  constexpr double EPS = 1.e-24;

  // Create a square in 3D space as our "polyhedron", and a plane splitting it
  // between the middle. The polyhedron after a call to poly_clip_vertices
  // should have two new vertices where the plane crosses the 0-1 and 2-3 edges:
  //      |
  // 3 ---|--- 2        3 -- 2' -- 2
  // |    |    |        |          |
  // |    p->  |   -->  |          |
  // |    |    |        |          |
  // 0 ---|--- 1        0 -- 1' -- 1
  //      |
  // In addition, vertices 0 and 3 should be marked as clipped.

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  PolyhedronType* out_square = axom::allocate<PolyhedronType>(1);
  unsigned int* out_clipped = axom::allocate<unsigned int>(1);

  PolyhedronType square;
  square.addVertex({0.0, 0.0, 0.0});
  square.addVertex({1.0, 0.0, 0.0});
  square.addVertex({1.0, 1.0, 0.0});
  square.addVertex({0.0, 1.0, 0.0});

  square.addNeighbors(0, {1, 3});
  square.addNeighbors(1, {0, 2});
  square.addNeighbors(2, {1, 3});
  square.addNeighbors(3, {0, 2});

  axom::copy(out_square, &square, sizeof(PolyhedronType));

  axom::for_all<ExecPolicy>(
    0,
    1,
    AXOM_LAMBDA(int /* idx */) {
      PlaneType plane(VectorType {1, 0, 0}, PointType {0.5, 0.0, 0.0});

      out_clipped[0] = 0;

      axom::primal::detail::poly_clip_vertices(out_square[0], plane, EPS, out_clipped[0]);
    });

  PolyhedronType clippedSquare;
  axom::copy(&clippedSquare, out_square, sizeof(PolyhedronType));

  unsigned int out_clipped_host;
  axom::copy(&out_clipped_host, out_clipped, sizeof(unsigned int));

  EXPECT_EQ(clippedSquare.numVertices(), 6);
  // Check that existing vertices were not modified
  EXPECT_EQ(clippedSquare[0], square[0]);
  EXPECT_EQ(clippedSquare[1], square[1]);
  EXPECT_EQ(clippedSquare[2], square[2]);
  EXPECT_EQ(clippedSquare[3], square[3]);

  int lo_idx = -1, hi_idx = -1;
  for(int i = 4; i < 6; i++)
  {
    if(clippedSquare[i] == PointType {0.5, 0.0, 0.0})
    {
      lo_idx = i;
    }
    else if(clippedSquare[i] == PointType {0.5, 1.0, 0.0})
    {
      hi_idx = i;
    }
  }

  // Check that plane-crossing vertices were generated
  EXPECT_NE(lo_idx, -1);
  EXPECT_NE(hi_idx, -1);

  // Check that vertices outside the plane were marked as clipped
  EXPECT_EQ(out_clipped_host, (1 | (1 << 3)));

  // Generate sets of expected neighbors
  std::vector<std::set<int>> expectedNbrs(6);
  expectedNbrs[0] = {3, lo_idx};
  expectedNbrs[1] = {2, lo_idx};
  expectedNbrs[2] = {1, hi_idx};
  expectedNbrs[3] = {0, hi_idx};
  expectedNbrs[lo_idx] = {0, 1};
  expectedNbrs[hi_idx] = {2, 3};

  // Check connectivity
  for(int iv = 0; iv < 6; iv++)
  {
    // Each vertex should only have two neighbors
    EXPECT_EQ(clippedSquare.getNumNeighbors(iv), 2);
    std::set<int> resultNbrs {clippedSquare.getNeighbors(iv)[0], clippedSquare.getNeighbors(iv)[1]};
    // Check that neighbors match expected connectivity
    EXPECT_EQ(expectedNbrs[iv], resultNbrs);
  }

  axom::deallocate(out_square);
  axom::deallocate(out_clipped);
  axom::setDefaultAllocator(current_allocator);
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
template <typename ExecPolicy>
void check_hex_tet_clip(double EPS)
{
  using namespace Primal3D;

  constexpr bool CHECK_ORIENTATION = true;

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Set new default to device if available
  axom::setDefaultAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  // Allocate memory for shapes
  TetrahedronType* tet = axom::allocate<TetrahedronType>(1);
  HexahedronType* hex = axom::allocate<HexahedronType>(1);
  PolyhedronType* res = axom::allocate<PolyhedronType>(1);

  // Shapes on host
  TetrahedronType tet_host;
  HexahedronType hex_host;
  PolyhedronType res_host;

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      tet[0] = TetrahedronType(PointType {1, 0, 0},
                               PointType {1, 1, 0},
                               PointType {0, 1, 0},
                               PointType {1, 0, 1});

      hex[0] = HexahedronType(PointType {0, 0, 0},
                              PointType {1, 0, 0},
                              PointType {1, 1, 0},
                              PointType {0, 1, 0},
                              PointType {0, 0, 1},
                              PointType {1, 0, 1},
                              PointType {1, 1, 1},
                              PointType {0, 1, 1});

      res[i] = axom::primal::clip(hex[i], tet[i]);
    });

  axom::copy(&tet_host, tet, sizeof(TetrahedronType));
  axom::copy(&hex_host, hex, sizeof(HexahedronType));
  axom::copy(&res_host, res, sizeof(PolyhedronType));

  EXPECT_NEAR(0.1666, res_host.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(hex_host, tet_host), EPS);

  // Test tryFixOrientation optional parameter using shapes with negative volumes
  axom::utilities::swap<PointType>(tet_host[1], tet_host[2]);
  axom::utilities::swap<PointType>(hex_host[1], hex_host[3]);
  axom::utilities::swap<PointType>(hex_host[5], hex_host[7]);

  EXPECT_LT(tet_host.signedVolume(), 0.0);
  EXPECT_LT(hex_host.signedVolume(), 0.0);

  axom::copy(tet, &tet_host, sizeof(TetrahedronType));
  axom::copy(hex, &hex_host, sizeof(HexahedronType));

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) { res[i] = axom::primal::clip(hex[i], tet[i], EPS, CHECK_ORIENTATION); });

  axom::copy(&res_host, res, sizeof(PolyhedronType));

  EXPECT_NEAR(0.1666, res_host.volume(), EPS);
  EXPECT_NEAR(0.1666,
              axom::primal::intersection_volume<double>(hex_host, tet_host, EPS, CHECK_ORIENTATION),
              EPS);

  axom::deallocate(tet);
  axom::deallocate(hex);
  axom::deallocate(res);

  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecPolicy>
void check_oct_tet_clip(double EPS)
{
  using namespace Primal3D;

  constexpr bool CHECK_ORIENTATION = true;

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Set new default to device if available
  axom::setDefaultAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  // Allocate memory for shapes
  TetrahedronType* tet = axom::allocate<TetrahedronType>(1);
  OctahedronType* oct = axom::allocate<OctahedronType>(1);
  PolyhedronType* res = axom::allocate<PolyhedronType>(1);

  // Shapes on host
  TetrahedronType tet_host;
  OctahedronType oct_host;
  PolyhedronType res_host;

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      tet[0] = TetrahedronType(PointType {1, 0, 0},
                               PointType {1, 1, 0},
                               PointType {0, 1, 0},
                               PointType {1, 0, 1});

      oct[0] = OctahedronType(PointType {1, 0, 0},
                              PointType {1, 1, 0},
                              PointType {0, 1, 0},
                              PointType {0, 1, 1},
                              PointType {0, 0, 1},
                              PointType {1, 0, 1});

      res[i] = axom::primal::clip(oct[i], tet[i]);
    });

  axom::copy(&tet_host, tet, sizeof(TetrahedronType));
  axom::copy(&oct_host, oct, sizeof(OctahedronType));
  axom::copy(&res_host, res, sizeof(PolyhedronType));

  EXPECT_NEAR(0.1666, res_host.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(oct_host, tet_host), EPS);

  // Test tryFixOrientation optional parameter using shapes with negative volumes
  axom::utilities::swap<PointType>(tet_host[1], tet_host[2]);
  axom::utilities::swap<PointType>(oct_host[1], oct_host[2]);
  axom::utilities::swap<PointType>(oct_host[4], oct_host[5]);

  EXPECT_LT(tet_host.signedVolume(), 0.0);

  axom::copy(tet, &tet_host, sizeof(TetrahedronType));
  axom::copy(oct, &oct_host, sizeof(OctahedronType));

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) { res[i] = axom::primal::clip(oct[i], tet[i], EPS, CHECK_ORIENTATION); });

  axom::copy(&res_host, res, sizeof(PolyhedronType));

  EXPECT_NEAR(0.1666, res_host.volume(), EPS);
  EXPECT_NEAR(0.1666,
              axom::primal::intersection_volume<double>(oct_host, tet_host, EPS, CHECK_ORIENTATION),
              EPS);

  axom::deallocate(tet);
  axom::deallocate(oct);
  axom::deallocate(res);

  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecPolicy>
void check_tet_tet_clip(double EPS)
{
  using namespace Primal3D;

  constexpr bool CHECK_ORIENTATION = true;

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Set new default to device if available
  axom::setDefaultAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  // Allocate memory for shapes
  TetrahedronType* tet1 = axom::allocate<TetrahedronType>(1);
  TetrahedronType* tet2 = axom::allocate<TetrahedronType>(1);
  PolyhedronType* res = axom::allocate<PolyhedronType>(1);

  // Shapes on host
  TetrahedronType tet1_host;
  TetrahedronType tet2_host;
  PolyhedronType res_host;

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      tet1[0] = TetrahedronType(PointType {1, 0, 0},
                                PointType {1, 1, 0},
                                PointType {0, 1, 0},
                                PointType {1, 1, 1});

      tet2[0] = TetrahedronType(PointType {0, 0, 0},
                                PointType {1, 0, 0},
                                PointType {1, 1, 0},
                                PointType {1, 1, 1});

      res[i] = axom::primal::clip(tet1[i], tet2[i]);
    });

  axom::copy(&tet1_host, tet1, sizeof(TetrahedronType));
  axom::copy(&tet2_host, tet2, sizeof(TetrahedronType));
  axom::copy(&res_host, res, sizeof(PolyhedronType));

  EXPECT_NEAR(0.0833, res_host.volume(), EPS);
  EXPECT_NEAR(0.0833, axom::primal::intersection_volume<double>(tet1_host, tet2_host), EPS);

  // Test tryFixOrientation optional parameter using shapes with negative volumes
  axom::utilities::swap<PointType>(tet1_host[1], tet1_host[2]);
  axom::utilities::swap<PointType>(tet2_host[1], tet2_host[2]);

  EXPECT_LT(tet1_host.signedVolume(), 0.0);
  EXPECT_LT(tet2_host.signedVolume(), 0.0);

  axom::copy(tet1, &tet1_host, sizeof(TetrahedronType));
  axom::copy(tet2, &tet2_host, sizeof(TetrahedronType));

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) { res[i] = axom::primal::clip(tet1[i], tet2[i], EPS, CHECK_ORIENTATION); });

  axom::copy(&res_host, res, sizeof(PolyhedronType));

  EXPECT_NEAR(0.0833, res_host.volume(), EPS);
  EXPECT_NEAR(0.0833,
              axom::primal::intersection_volume<double>(tet1_host, tet2_host, EPS, CHECK_ORIENTATION),
              EPS);

  axom::deallocate(tet1);
  axom::deallocate(tet2);
  axom::deallocate(res);

  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecPolicy>
void check_polygon_polygon_clip(double EPS)
{
  // Will be clipping two triangles, max number of vertices for output polygon
  // is 3 + 3 = 6
  const int MAX_NUM_VERTS = 6;

  constexpr bool CHECK_SIGN = true;

  using PolygonStatic2D =
    axom::primal::Polygon<double, 2, axom::primal::PolygonArray::Static, MAX_NUM_VERTS>;
  using Point2D = axom::primal::Point<double, 2>;

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int kernel_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  axom::Array<PolygonStatic2D> subject_polygon_device(1, 1, kernel_allocator);
  auto subject_polygon_view = subject_polygon_device.view();

  axom::Array<PolygonStatic2D> clip_polygon_device(1, 1, kernel_allocator);
  auto clip_polygon_view = clip_polygon_device.view();

  axom::Array<PolygonStatic2D> output_polygon_device(1, 1, kernel_allocator);
  auto output_polygon_view = output_polygon_device.view();

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      subject_polygon_view[i] =
        PolygonStatic2D({Point2D({0.0, 0.0}), Point2D({1.0, 0.0}), Point2D({0.5, 1})});

      clip_polygon_view[i] = PolygonStatic2D(
        {Point2D({0.0, 2.0 / 3.0}), Point2D({0.5, -1.0 / 3.0}), Point2D({1.0, 2.0 / 3.0})});

      output_polygon_view[i] = axom::primal::clip(subject_polygon_view[i], clip_polygon_view[i], EPS);
    });

  // Copy output polygon back to host
  axom::Array<PolygonStatic2D> output_polygon_host =
    axom::Array<PolygonStatic2D>(output_polygon_device, host_allocator);

  EXPECT_NEAR(0.3333, output_polygon_host[0].signedArea(), EPS);
  EXPECT_EQ(output_polygon_host[0].numVertices(), 6);

  // Test tryFixOrientation optional parameter using same polygons with
  // negative volumes (clockwise ordering)
  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      subject_polygon_view[i].clear();
      clip_polygon_view[i].clear();
      output_polygon_view[i].clear();

      subject_polygon_view[i].addVertex(Point2D({0.5, 1}));
      subject_polygon_view[i].addVertex(Point2D({1.0, 0.0}));
      subject_polygon_view[i].addVertex(Point2D({0.0, 0.0}));

      clip_polygon_view[i].addVertex(Point2D({0.5, -1.0 / 3.0}));
      clip_polygon_view[i].addVertex(Point2D({0.0, 2.0 / 3.0}));
      clip_polygon_view[i].addVertex(Point2D({1.0, 2.0 / 3.0}));

      output_polygon_view[i] =
        axom::primal::clip(subject_polygon_view[i], clip_polygon_view[i], EPS, CHECK_SIGN);
    });

  // Copy output polygon back to host
  output_polygon_host = axom::Array<PolygonStatic2D>(output_polygon_device, host_allocator);

  EXPECT_NEAR(0.3333, output_polygon_host[0].signedArea(), EPS);
  EXPECT_EQ(output_polygon_host[0].numVertices(), 6);
}

TEST(primal_clip, unit_poly_clip_vertices_sequential) { unit_check_poly_clip<axom::SEQ_EXEC>(); }

TEST(primal_clip, clip_hex_tet_sequential)
{
  constexpr double EPS = 1e-4;
  check_hex_tet_clip<axom::SEQ_EXEC>(EPS);
}

TEST(primal_clip, clip_oct_tet_sequential)
{
  constexpr double EPS = 1e-4;
  check_oct_tet_clip<axom::SEQ_EXEC>(EPS);
}

TEST(primal_clip, clip_tet_tet_sequential)
{
  constexpr double EPS = 1e-4;
  check_tet_tet_clip<axom::SEQ_EXEC>(EPS);
}

TEST(primal_clip, clip_polygon_polygon_sequential)
{
  constexpr double EPS = 1e-4;
  check_polygon_polygon_clip<axom::SEQ_EXEC>(EPS);
}

  #ifdef AXOM_USE_OPENMP
TEST(primal_clip, unit_poly_clip_vertices_omp) { unit_check_poly_clip<axom::OMP_EXEC>(); }

TEST(primal_clip, clip_hex_tet_omp)
{
  constexpr double EPS = 1e-4;
  check_hex_tet_clip<axom::OMP_EXEC>(EPS);
}

TEST(primal_clip, clip_oct_tet_omp)
{
  constexpr double EPS = 1e-4;
  check_oct_tet_clip<axom::OMP_EXEC>(EPS);
}

TEST(primal_clip, clip_tet_tet_omp)
{
  constexpr double EPS = 1e-4;
  check_tet_tet_clip<axom::OMP_EXEC>(EPS);
}

TEST(primal_clip, clip_polygon_polygon_omp)
{
  constexpr double EPS = 1e-4;
  check_polygon_polygon_clip<axom::OMP_EXEC>(EPS);
}
  #endif /* AXOM_USE_OPENMP */

  #if defined(AXOM_USE_CUDA)
AXOM_CUDA_TEST(primal_clip, unit_poly_clip_vertices_cuda)
{
  unit_check_poly_clip<axom::CUDA_EXEC<256>>();
}

AXOM_CUDA_TEST(primal_clip, clip_hex_tet_cuda)
{
  constexpr double EPS = 1e-4;
  check_hex_tet_clip<axom::CUDA_EXEC<256>>(EPS);
}

AXOM_CUDA_TEST(primal_clip, clip_oct_tet_cuda)
{
  constexpr double EPS = 1e-4;
  check_oct_tet_clip<axom::CUDA_EXEC<256>>(EPS);
}

TEST(primal_clip, clip_tet_tet_cuda)
{
  constexpr double EPS = 1e-4;
  check_tet_tet_clip<axom::CUDA_EXEC<256>>(EPS);
}

TEST(primal_clip, clip_polygon_polygon_cuda)
{
  constexpr double EPS = 1e-4;
  check_polygon_polygon_clip<axom::CUDA_EXEC<256>>(EPS);
}
  #endif /* AXOM_USE_CUDA */

  #if defined(AXOM_USE_HIP)
TEST(primal_clip, unit_poly_clip_vertices_hip) { unit_check_poly_clip<axom::HIP_EXEC<256>>(); }

TEST(primal_clip, clip_hex_tet_hip)
{
  constexpr double EPS = 1e-4;
  check_hex_tet_clip<axom::HIP_EXEC<256>>(EPS);
}

TEST(primal_clip, clip_oct_tet_hip)
{
  constexpr double EPS = 1e-4;
  check_oct_tet_clip<axom::HIP_EXEC<256>>(EPS);
}

TEST(primal_clip, clip_tet_tet_hip)
{
  constexpr double EPS = 1e-4;
  check_tet_tet_clip<axom::HIP_EXEC<256>>(EPS);
}

TEST(primal_clip, clip_polygon_polygon_hip)
{
  constexpr double EPS = 1e-4;
  check_polygon_polygon_clip<axom::HIP_EXEC<256>>(EPS);
}
  #endif /* AXOM_USE_HIP */

#endif /* AXOM_USE_RAJA && AXOM_USE_UMPIRE */

// Tetrahedron does not clip hexahedron.
TEST(primal_clip, hex_tet_clip_nonintersect)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {-1, -1, -1},
                      PointType {-1, 0, 0},
                      PointType {-1, -1, 0},
                      PointType {0, 0, 0});
  HexahedronType hex(PointType {1, 0, 0},
                     PointType {2, 0, 0},
                     PointType {2, 1, 0},
                     PointType {1, 1, 0},
                     PointType {1, 0, 1},
                     PointType {2, 0, 1},
                     PointType {2, 1, 1},
                     PointType {1, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);
  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(hex, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, hex));
  EXPECT_EQ(axom::primal::clip(tet, hex).volume(), poly.volume());
}

// Tetrahedron is encapsulated by the hexahedron
TEST(primal_clip, hex_tet_clip_encapsulate)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {1, 0, 0},
                      PointType {1, 1, 0},
                      PointType {0, 1, 0},
                      PointType {1, 0, 1});
  HexahedronType hex(PointType {0, 0, 0},
                     PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1},
                     PointType {1, 1, 1},
                     PointType {0, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(hex, tet), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(tet, hex), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, hex).volume(), poly.volume(), EPS);
}

// Hexahedron is encapsulated inside the tetrahedron
TEST(primal_clip, hex_tet_clip_encapsulate_inv)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {0, 0, 0},
                      PointType {0, 3, 0},
                      PointType {0, 0, 3},
                      PointType {3, 0, 0});
  HexahedronType hex(PointType {0, 0, 0},
                     PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1},
                     PointType {1, 1, 1},
                     PointType {0, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  // Expected result should be 1.0, volume of hex.
  EXPECT_NEAR(1.0, poly.volume(), EPS);
  EXPECT_NEAR(1.0, axom::primal::intersection_volume<double>(hex, tet), EPS);
  EXPECT_NEAR(1.0, axom::primal::intersection_volume<double>(tet, hex), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, hex).volume(), poly.volume(), EPS);
}

// Half of the hexahedron is clipped by the tetrahedron
TEST(primal_clip, hex_tet_clip_half)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {0, 0, 0},
                      PointType {0, 3, 0},
                      PointType {0, 0, 3},
                      PointType {3, 0, 0});
  HexahedronType hex(PointType {-1, 0, 0},
                     PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {-1, 1, 0},
                     PointType {-1, 0, 1},
                     PointType {1, 0, 1},
                     PointType {1, 1, 1},
                     PointType {-1, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  // Expected result should be 1.0, half the volume of hex.
  EXPECT_NEAR(1.0, poly.volume(), EPS);
  EXPECT_NEAR(1.0, axom::primal::intersection_volume<double>(hex, tet), EPS);
  EXPECT_NEAR(1.0, axom::primal::intersection_volume<double>(tet, hex), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, hex).volume(), poly.volume(), EPS);
}

// Hexahedron is adjacent to tetrahedron
TEST(primal_clip, hex_tet_clip_adjacent)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {0, -1, 0},
                      PointType {0, 0, 1},
                      PointType {1, 0, 1},
                      PointType {1, 0, 0});
  HexahedronType hex(PointType {0, 0, 0},
                     PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1},
                     PointType {1, 1, 1},
                     PointType {0, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(hex, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, hex));
  EXPECT_EQ(axom::primal::clip(tet, hex).volume(), poly.volume());
}

// Tetrahedron clips hexahedron at a single vertex
TEST(primal_clip, hex_tet_clip_point)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {0, -1, 0},
                      PointType {0, 0, 1},
                      PointType {1, 0, 1},
                      PointType {1, 0, 0});
  HexahedronType hex(PointType {-1, -2, -1},
                     PointType {1, -2, -1},
                     PointType {1, -1, -1},
                     PointType {-1, -1, -1},
                     PointType {-1, -2, 1},
                     PointType {1, -2, 1},
                     PointType {1, -1, 1},
                     PointType {-1, -1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(hex, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, hex));
  EXPECT_EQ(axom::primal::clip(tet, hex).volume(), poly.volume());
}

// Tetrahedron with a negative volume clips hexahedron
TEST(primal_clip, hex_tet_negative_tet_vol)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;
  constexpr bool CHECK_ORIENTATION = true;

  TetrahedronType tet(PointType {0.0774795, -11.3217, -33.5237},
                      PointType {0.482086, -12.283, -33.5237},
                      PointType {-0.503954, -12.216, -33.5237},
                      PointType {0.0368068, -12.1403, -33.4339});
  HexahedronType hex(PointType {0, -12.1295, -33.6323},
                     PointType {0, -12.2774, -33.6323},
                     PointType {0.185429, -12.276, -33.6323},
                     PointType {0.183195, -12.1281, -33.6323},
                     PointType {0, -12.1295, -33.5186},
                     PointType {0, -12.2774, -33.5186},
                     PointType {0.185429, -12.276, -33.5186},
                     PointType {0.183195, -12.1281, -33.5186});

  // Clipped polyhedron has volume of zero due to tetrahedron orientation
  PolyhedronType zeroPoly = axom::primal::clip(hex, tet);

  EXPECT_EQ(0.0, zeroPoly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(hex, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, hex));
  EXPECT_EQ(axom::primal::clip(tet, hex).volume(), zeroPoly.volume());

  // Run clip operation with orientation flag enabled to fix tetrahedron
  PolyhedronType fixedPoly = axom::primal::clip(hex, tet, EPS, CHECK_ORIENTATION);

  EXPECT_NEAR(0.00011, fixedPoly.volume(), EPS);
  EXPECT_NEAR(0.00011,
              axom::primal::intersection_volume<double>(hex, tet, EPS, CHECK_ORIENTATION),
              EPS);
  EXPECT_NEAR(0.00011,
              axom::primal::intersection_volume<double>(tet, hex, EPS, CHECK_ORIENTATION),
              EPS);
  EXPECT_EQ(axom::primal::clip(tet, hex, EPS, CHECK_ORIENTATION).volume(), fixedPoly.volume());
}

// Tetrahedron does not clip octahedron.
TEST(primal_clip, oct_tet_clip_nonintersect)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {-1, -1, -1},
                      PointType {-1, 0, 0},
                      PointType {-1, -1, 0},
                      PointType {0, 0, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);
  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(oct, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, oct));
  EXPECT_EQ(axom::primal::clip(tet, oct).volume(), poly.volume());
}

// Tetrahedron is encapsulated by the octahedron
TEST(primal_clip, oct_tet_clip_encapsulate)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {1, 0, 0},
                      PointType {1, 1, 0},
                      PointType {0, 1, 0},
                      PointType {1, 0, 1});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(oct, tet), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(tet, oct), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct).volume(), poly.volume(), EPS);
}

// Octahedron is encapsulated inside the tetrahedron
TEST(primal_clip, oct_tet_clip_encapsulate_inv)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {0, 0, 0},
                      PointType {0, 2, 0},
                      PointType {0, 0, 2},
                      PointType {2, 0, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.6666, volume of oct.
  EXPECT_NEAR(0.6666, poly.volume(), EPS);
  EXPECT_NEAR(0.6666, axom::primal::intersection_volume<double>(oct, tet), EPS);
  EXPECT_NEAR(0.6666, axom::primal::intersection_volume<double>(tet, oct), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct).volume(), poly.volume(), EPS);
}

// Half of the octahedron is clipped by the tetrahedron
TEST(primal_clip, oct_tet_clip_half)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {0.5, 0.5, 2},
                      PointType {2, -1, 0},
                      PointType {-1, -1, 0},
                      PointType {-1, 2, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.3333, half the volume of oct.
  EXPECT_NEAR(0.3333, poly.volume(), EPS);
  EXPECT_NEAR(0.3333, axom::primal::intersection_volume<double>(oct, tet), EPS);
  EXPECT_NEAR(0.3333, axom::primal::intersection_volume<double>(tet, oct), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct).volume(), poly.volume(), EPS);
}

// Octahedron is adjacent to tetrahedron
TEST(primal_clip, oct_tet_clip_adjacent)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {0, -1, 0},
                      PointType {0, 0, 1},
                      PointType {1, 0, 1},
                      PointType {1, 0, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(oct, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, oct));
  EXPECT_EQ(axom::primal::clip(tet, oct).volume(), poly.volume());
}

// Tetrahedron clips octahedron at a single vertex
TEST(primal_clip, oct_tet_clip_point)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {-1, -1, 0},
                      PointType {-0.5, 0.5, 0},
                      PointType {0, 0, 2},
                      PointType {0.5, -0.5, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(oct, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, oct));
  EXPECT_EQ(axom::primal::clip(tet, oct).volume(), poly.volume());
}

TEST(primal_clip, oct_tet_clip_special_case_1)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;
  constexpr bool CHECK_ORIENTATION = true;

  TetrahedronType tet(PointType {0.5, 0.5, 0.5},
                      PointType {1, 1, 0},
                      PointType {1, 0, 0},
                      PointType {0.5, 0.5, 0});

  OctahedronType oct(PointType {0.5, 0.853553, 0.146447},
                     PointType {0.853553, 0.853553, 0.5},
                     PointType {0.853553, 0.5, 0.146447},
                     PointType {1, 0.5, 0.5},
                     PointType {0.5, 0.5, 0},
                     PointType {0.5, 1, 0.5});

  // NOTE: Order of vertices 1,2 and 4,5 are flipped due
  // to winding being opposite of what volume() expects
  PolyhedronType octPoly;
  octPoly.addVertex(PointType {0.5, 0.853553, 0.146447});
  octPoly.addVertex(PointType {0.853553, 0.5, 0.146447});
  octPoly.addVertex(PointType {0.853553, 0.853553, 0.5});
  octPoly.addVertex(PointType {1, 0.5, 0.5});
  octPoly.addVertex(PointType {0.5, 1, 0.5});
  octPoly.addVertex(PointType {0.5, 0.5, 0});

  octPoly.addNeighbors(0, {1, 5, 4, 2});
  octPoly.addNeighbors(1, {0, 2, 3, 5});
  octPoly.addNeighbors(2, {0, 4, 3, 1});
  octPoly.addNeighbors(3, {1, 2, 4, 5});
  octPoly.addNeighbors(4, {0, 5, 3, 2});
  octPoly.addNeighbors(5, {0, 1, 3, 4});

  EXPECT_NEAR(0.0251, octPoly.volume(), EPS);

  PolyhedronType poly = axom::primal::clip(oct, tet, EPS, CHECK_ORIENTATION);

  EXPECT_NEAR(0.0041, poly.volume(), EPS);
  EXPECT_NEAR(0.0041, axom::primal::intersection_volume<double>(oct, tet, EPS, CHECK_ORIENTATION), EPS);
  EXPECT_NEAR(0.0041, axom::primal::intersection_volume<double>(tet, oct, EPS, CHECK_ORIENTATION), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct, EPS, CHECK_ORIENTATION).volume(), poly.volume(), EPS);
}

TEST(primal_clip, oct_tet_clip_special_case_2)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;
  constexpr bool CHECK_ORIENTATION = true;

  TetrahedronType tet(PointType {0.5, 0.5, 0.5},
                      PointType {0, 1, 0},
                      PointType {1, 1, 0},
                      PointType {0.5, 0.5, 0});

  OctahedronType oct(PointType {0.5, 0.853553, 0.146447},
                     PointType {0.853553, 0.853553, 0.5},
                     PointType {0.853553, 0.5, 0.146447},
                     PointType {1, 0.5, 0.5},
                     PointType {0.5, 0.5, 0},
                     PointType {0.5, 1, 0.5});

  // NOTE: Order of vertices 1,2 and 4,5 are flipped due
  // to winding being opposite of what volume() expects
  PolyhedronType octPoly;
  octPoly.addVertex(PointType {0.5, 0.853553, 0.146447});
  octPoly.addVertex(PointType {0.853553, 0.5, 0.146447});
  octPoly.addVertex(PointType {0.853553, 0.853553, 0.5});
  octPoly.addVertex(PointType {1, 0.5, 0.5});
  octPoly.addVertex(PointType {0.5, 1, 0.5});
  octPoly.addVertex(PointType {0.5, 0.5, 0});

  octPoly.addNeighbors(0, {1, 5, 4, 2});
  octPoly.addNeighbors(1, {0, 2, 3, 5});
  octPoly.addNeighbors(2, {0, 4, 3, 1});
  octPoly.addNeighbors(3, {1, 2, 4, 5});
  octPoly.addNeighbors(4, {0, 5, 3, 2});
  octPoly.addNeighbors(5, {0, 1, 3, 4});

  EXPECT_NEAR(0.0251, octPoly.volume(), EPS);

  PolyhedronType poly = axom::primal::clip(oct, tet, EPS, CHECK_ORIENTATION);

  EXPECT_NEAR(0.0041, poly.volume(), EPS);
  EXPECT_NEAR(0.0041, axom::primal::intersection_volume<double>(oct, tet, EPS, CHECK_ORIENTATION), EPS);
  EXPECT_NEAR(0.0041, axom::primal::intersection_volume<double>(tet, oct, EPS, CHECK_ORIENTATION), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct, EPS, CHECK_ORIENTATION).volume(), poly.volume(), EPS);
}

// Tetrahedron does not clip tetrahedron.
TEST(primal_clip, tet_tet_clip_nonintersect)
{
  using namespace Primal3D;

  {
    TetrahedronType tet1(PointType {-1, -1, -1},
                         PointType {-1, 0, 0},
                         PointType {-1, -1, 0},
                         PointType {0, 0, 0});
    EXPECT_TRUE(tet1.signedVolume() > 0.);

    TetrahedronType tet2(PointType {1, 0, 0},
                         PointType {1, 1, 0},
                         PointType {0, 1, 0},
                         PointType {1, 0, 1});
    EXPECT_TRUE(tet2.signedVolume() > 0.);

    PolyhedronType poly = axom::primal::clip(tet1, tet2);
    EXPECT_EQ(0.0, poly.volume());
    EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet1, tet2));
  }

  // User-provided test case; tetrahedra do not overlap
  {
    TetrahedronType tet1(PointType {0, 1, 0},
                         PointType {-1, 0, 0},
                         PointType {0, -1, 0},
                         PointType {-0.25, 0, 0.25});
    EXPECT_TRUE(tet1.signedVolume() > 0.);

    TetrahedronType tet2(PointType {0.74899999999999922, 0, 1},
                         PointType {1.7489999999999992, 0, 0},
                         PointType {0.74899999999999922, -1, 0},
                         PointType {0.99899999999999922, 0, 0.25});
    EXPECT_TRUE(tet2.signedVolume() > 0.);

    PolyhedronType poly = axom::primal::clip(tet1, tet2);
    EXPECT_EQ(0.0, poly.volume());
    EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet1, tet2));
  }
}

// Tetrahedron is adjacent to tetrahedron
TEST(primal_clip, tet_tet_clip_adjacent)
{
  using namespace Primal3D;

  TetrahedronType tet1(PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {0, 1, 0},
                       PointType {1, 0, 1});

  TetrahedronType tet2(PointType {1, 0, 1},
                       PointType {0, 1, 0},
                       PointType {1, 0, 0},
                       PointType {0, 0, 0});

  PolyhedronType poly = axom::primal::clip(tet1, tet2);
  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet1, tet2));
}

// Tetrahedron clips tetrahedron at a single vertex
TEST(primal_clip, tet_tet_clip_point)
{
  using namespace Primal3D;

  TetrahedronType tet1(PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {0, 1, 0},
                       PointType {1, 0, 1});

  TetrahedronType tet2(PointType {0, 1, 0},
                       PointType {0, 0, 0},
                       PointType {-1, 0, 0},
                       PointType {0, 0, 1});

  PolyhedronType poly = axom::primal::clip(tet1, tet2);
  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet1, tet2));
}

// Tetrahedrons are the same
TEST(primal_clip, tet_tet_equal)
{
  using namespace Primal3D;
  const double EPS = 1e-4;

  TetrahedronType tet(PointType {1, 0, 0},
                      PointType {1, 1, 0},
                      PointType {0, 1, 0},
                      PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(tet, tet);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(tet, tet), EPS);
}

// Tetrahedron is encapsulated inside the other tetrahedron
TEST(primal_clip, tet_tet_encapsulate)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet1(PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {0, 1, 0},
                       PointType {1, 0, 1});

  TetrahedronType tet2(PointType {3, 0, 0},
                       PointType {0, 3, 0},
                       PointType {-3, 0, 0},
                       PointType {0, 0, 3});

  PolyhedronType poly = axom::primal::clip(tet1, tet2);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(tet1, tet2), EPS);
}

// Half of the tetrahedron is clipped by the other tetrahedron
TEST(primal_clip, tet_tet_half)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet1(PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {0, 1, 0},
                       PointType {1, 1, 1});

  TetrahedronType tet2(PointType {0, 0, 0},
                       PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {1, 1, 1});

  PolyhedronType poly = axom::primal::clip(tet1, tet2);

  // Expected result should be 0.666 / 4 / 2 = 0.0833, volume of tet.
  EXPECT_NEAR(0.0833, poly.volume(), EPS);
  EXPECT_NEAR(0.0833, axom::primal::intersection_volume<double>(tet1, tet2), EPS);
}

// Half of the octahedron is clipped by the tetrahedron.
// Then, we split the octahedron into 8 tetrahedrons,
// and verify the total volume from clipping the 8 split tetrahedrons
// by the starting tetrahedron are the same as with the octahedron.
TEST(primal_clip, tet_tet_clip_split)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;
  constexpr bool CHECK_ORIENTATION = true;

  TetrahedronType tet(PointType {0.5, 0.5, 2},
                      PointType {2, -1, 0},
                      PointType {-1, -1, 0},
                      PointType {-1, 2, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.3333, half the volume of oct.
  EXPECT_NEAR(0.3333, poly.volume(), EPS);
  EXPECT_NEAR(0.3333, axom::primal::intersection_volume<double>(oct, tet), EPS);
  EXPECT_NEAR(0.3333, axom::primal::intersection_volume<double>(tet, oct), EPS);

  // Split the octahedron into 8 tetrahedrons
  double tet_volumes = 0.0;

  axom::Array<TetrahedronType> split_tets(8);

  axom::primal::split(oct, split_tets);

  for(int i = 0; i < split_tets.size(); i++)
  {
    tet_volumes += (axom::primal::clip(split_tets[i], tet, EPS, CHECK_ORIENTATION)).volume();
  }

  // Expected result should still be 0.3333
  EXPECT_NEAR(0.3333, tet_volumes, EPS);
}

TEST(primal_clip, tet_tet_clip_special_case_1)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_ORIENTATION = true;

  // Tets do not intersect, but share a face
  TetrahedronType tet1(PointType {0.5, 0.5, -0.125},
                       PointType {0, 0, -0.25},
                       PointType {1, 0, -0.25},
                       PointType {0.5, 0, -0.125});

  TetrahedronType tet2(PointType {0.1875, 0.0625, -0.234375},
                       PointType {0.125, 0.125, -0.25},
                       PointType {0.125, 0, -0.25},
                       PointType {0.125, 0.0625, -0.234375});

  PolyhedronType poly = axom::primal::clip(tet1, tet2, EPS, CHECK_ORIENTATION);

  EXPECT_NEAR(0.00, poly.volume(), EPS);
  EXPECT_NEAR(0.00, axom::primal::intersection_volume<double>(tet2, tet1, EPS, CHECK_ORIENTATION), EPS);
  EXPECT_NEAR(0.00, axom::primal::intersection_volume<double>(tet1, tet2, EPS, CHECK_ORIENTATION), EPS);
}

TEST(primal_clip, tet_plane_intersect_none_below)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Plane intersects one vertex of tet
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, 0.0, -1.0}, PointType {0.0, 0.0, -0.5});

  PolyhedronType poly = axom::primal::clip(tet, plane, EPS, CHECK_SIGN);

  EXPECT_NEAR(0.0, poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_intersect_none_above)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Plane intersects one vertex of tet
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, 0.0, 1.0}, PointType {0.0, 0.0, -0.5});

  PolyhedronType poly = axom::primal::clip(tet, plane, EPS, CHECK_SIGN);

  EXPECT_NEAR(tet.signedVolume(), poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_border_vertex_below)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Plane intersects one vertex of tet
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, 0.0, 1.0}, PointType {0.0, 0.0, 1.0});

  PolyhedronType poly = axom::primal::clip(tet, plane, EPS, CHECK_SIGN);

  EXPECT_NEAR(0.0, poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_border_vertex_above)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Plane intersects one vertex of tet
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, 0.0, -1.0}, PointType {0.0, 0.0, 1.0});

  PolyhedronType poly = axom::primal::clip(tet, plane, EPS, CHECK_SIGN);

  EXPECT_NEAR(tet.signedVolume(), poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_border_edge_below)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Plane intersects one edge of tet
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, -1.0, -1.0}, 0.0);

  PolyhedronType poly = axom::primal::clip(tet, plane, EPS, CHECK_SIGN);

  EXPECT_NEAR(0.0, poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_border_edge_above)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Plane intersects one edge of tet
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, 1.0, 1.0}, 0.0);

  PolyhedronType poly = axom::primal::clip(tet, plane, EPS, CHECK_SIGN);

  EXPECT_NEAR(tet.signedVolume(), poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_border_face_below)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Tet and plane do not intersect, but border each other
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, 0.0, -1.0}, 0.0);

  PolyhedronType poly = axom::primal::clip(plane, tet, EPS, CHECK_SIGN);

  EXPECT_NEAR(0.0, poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_border_face_above)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Tet and plane do not intersect, but border each other
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, 0.0, 1.0}, 0.0);

  PolyhedronType poly = axom::primal::clip(plane, tet, EPS, CHECK_SIGN);

  EXPECT_NEAR(tet.signedVolume(), poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_intersect_three_edges)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Plane intersects three edges of tet
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 0.0, 0.0},
                      PointType {0.0, 1.0, 0.0},
                      PointType {0.0, 0.0, 1.0});

  PlaneType plane(VectorType {-1.0, 0.0, 1.0}, 0.0);

  PolyhedronType poly = axom::primal::clip(plane, tet, EPS, CHECK_SIGN);

  EXPECT_NEAR(tet.signedVolume() / 2.0, poly.signedVolume(), EPS);
}

TEST(primal_clip, tet_plane_intersect_four_edges)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Plane intersects four edges of tet
  TetrahedronType tet(PointType {0.0, 0.0, 0.0},
                      PointType {1.0, 1.0, 0.0},
                      PointType {0.0, 1.0, 1.0},
                      PointType {1.0, 0.0, 1.0});

  PlaneType plane(VectorType {0.0, 0.0, 1.0}, 0.5);

  PolyhedronType poly = axom::primal::clip(plane, tet, EPS, CHECK_SIGN);

  EXPECT_NEAR(tet.signedVolume() / 2.0, poly.signedVolume(), EPS);
}

TEST(primal_clip, polyhedron_plane)
{
  // Use float (float uncovered compilation problems in Polyhedron)
  using Precision = float;
  using PlaneType = axom::primal::Plane<Precision, 3>;
  using VectorType = typename PlaneType::VectorType;

  constexpr double EPS = 2e-7;

  constexpr int nSides = 8;
  const auto poly = axom::primal::regular_prism<Precision>(nSides);
  PlaneType plane(VectorType {1.f, 0.f, 0.f}, 0.f);

  // Clip away half of the polyhedron.
  const auto clipped = axom::primal::clip(plane, poly, EPS);

  const double V = poly.signedVolume();
  const double V_2 = V / 2.0;
  const double Vc = clipped.signedVolume();

  // Entire polyhedron volume should be 2*sqrt(2)
  EXPECT_NEAR(V, 2 * sqrt(2.), EPS);

  // We clipped away half of the volume.
  EXPECT_NEAR(V_2, Vc, EPS);

  //----------------------------------------------------------------------------
  // Make the same shape tipped over (to test transforms)
  const auto poly2 = axom::primal::regular_prism<Precision>(
    nSides,
    1.f,
    1.f,
    axom::numerics::transforms::xRotation<Precision>(-M_PI / 2, 4));
  PlaneType plane2(VectorType {0.f, 1.f, 0.f}, 0.5f);

  // Clip away half of the polyhedron.
  const auto clipped2 = axom::primal::clip(plane2, poly2, EPS);
  const double Vc2 = clipped2.signedVolume();

  // We clipped away half of the volume.
  EXPECT_NEAR(Vc2, sqrt(2.), EPS);
}

TEST(primal_clip, empty_polygons)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;

  Polygon2D subjectPolygon;

  Polygon2D clipPolygon;

  Polygon2D poly = axom::primal::clip(subjectPolygon, clipPolygon);

  EXPECT_EQ(poly.isValid(), false);
}

TEST(primal_clip, polygon_intersects_polygon)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;
  using PolygonStatic2D = axom::primal::Polygon<double, 2, axom::primal::PolygonArray::Static>;
  using Point2D = axom::primal::Point<double, 2>;
  constexpr double EPS = 1e-10;
  constexpr bool CHECK_SIGN = true;

  // Expected counter-clockwise vertex ordering
  {
    Polygon2D subjectPolygon;
    subjectPolygon.addVertex(Point2D {0, 0});
    subjectPolygon.addVertex(Point2D {1, 0});
    subjectPolygon.addVertex(Point2D {1, 1});
    subjectPolygon.addVertex(Point2D {0, 1});

    Polygon2D clipPolygon;
    clipPolygon.addVertex(Point2D {0.5, -1});
    clipPolygon.addVertex(Point2D {1.5, -1});
    clipPolygon.addVertex(Point2D {1.5, 2});
    clipPolygon.addVertex(Point2D {0.5, 2});

    Polygon2D poly = axom::primal::clip(subjectPolygon, clipPolygon, EPS);

    EXPECT_NEAR(subjectPolygon.signedArea(), 1.0, EPS);
    EXPECT_NEAR(clipPolygon.signedArea(), 3.0, EPS);
    EXPECT_NEAR(poly.signedArea(), 0.5, EPS);
  }

  // Same polygons with static arrays
  {
    PolygonStatic2D subjectPolygon;
    subjectPolygon.addVertex(Point2D {0, 0});
    subjectPolygon.addVertex(Point2D {1, 0});
    subjectPolygon.addVertex(Point2D {1, 1});
    subjectPolygon.addVertex(Point2D {0, 1});

    PolygonStatic2D clipPolygon;
    clipPolygon.addVertex(Point2D {0.5, -1});
    clipPolygon.addVertex(Point2D {1.5, -1});
    clipPolygon.addVertex(Point2D {1.5, 2});
    clipPolygon.addVertex(Point2D {0.5, 2});

    PolygonStatic2D poly = axom::primal::clip(subjectPolygon, clipPolygon, EPS);

    EXPECT_NEAR(subjectPolygon.signedArea(), 1.0, EPS);
    EXPECT_NEAR(clipPolygon.signedArea(), 3.0, EPS);
    EXPECT_NEAR(poly.signedArea(), 0.5, EPS);
  }

  // Negative clockwise vertex ordering, use tryFixOrientation optional
  // parameter to fix.
  {
    Polygon2D subjectPolygon;
    subjectPolygon.addVertex(Point2D {0, 0});
    subjectPolygon.addVertex(Point2D {0, 1});
    subjectPolygon.addVertex(Point2D {1, 1});
    subjectPolygon.addVertex(Point2D {1, 0});

    Polygon2D clipPolygon;
    clipPolygon.addVertex(Point2D {0.5, -1});
    clipPolygon.addVertex(Point2D {0.5, 2});
    clipPolygon.addVertex(Point2D {1.5, 2});
    clipPolygon.addVertex(Point2D {1.5, -1});

    Polygon2D poly_fix_orientation = axom::primal::clip(subjectPolygon, clipPolygon, EPS, CHECK_SIGN);

    // Negative areas
    EXPECT_NEAR(subjectPolygon.signedArea(), -1.0, EPS);
    EXPECT_NEAR(clipPolygon.signedArea(), -3.0, EPS);

    // Positive area with tryFixOrientation flag enabled
    EXPECT_NEAR(poly_fix_orientation.signedArea(), 0.5, EPS);
  }

  // Edge case with 3 consecutive vertices having the orientations
  // ON_POSITIVE_SIDE, ON_BOUNDARY, ON_POSITIVE in respect to the
  // last edge of the clip polygon. Prior to fix, resulted in the
  // vertex on the boundary being added twice.
  {
    Polygon2D subjectPolygon;
    subjectPolygon.addVertex(Point2D {0.0, 0.0});
    subjectPolygon.addVertex(Point2D {1.0, 0.0});
    subjectPolygon.addVertex(Point2D {1.0, 1.0});

    Polygon2D clipPolygon;
    clipPolygon.addVertex(Point2D {0.0, 0.0});
    clipPolygon.addVertex(Point2D {1.0, 0.0});
    clipPolygon.addVertex(Point2D {0.0, 1.0});

    Polygon2D poly = axom::primal::clip(subjectPolygon, clipPolygon, EPS);

    EXPECT_NEAR(poly.signedArea(), 0.25, EPS);
    EXPECT_EQ(poly.numVertices(), 3);
  }

  // Non-clipping - polygons intersect at a line
  {
    Polygon2D subjectPolygon;
    subjectPolygon.addVertex(Point2D {0, 0});
    subjectPolygon.addVertex(Point2D {1, 0});
    subjectPolygon.addVertex(Point2D {1, 1});
    subjectPolygon.addVertex(Point2D {0, 1});

    Polygon2D clipPolygon;
    clipPolygon.addVertex(Point2D {-1, 0});
    clipPolygon.addVertex(Point2D {0, 0});
    clipPolygon.addVertex(Point2D {0, 1});
    clipPolygon.addVertex(Point2D {-1, 1});

    Polygon2D poly = axom::primal::clip(subjectPolygon, clipPolygon, EPS, CHECK_SIGN);

    EXPECT_EQ(poly.isValid(), false);
    EXPECT_EQ(poly.numVertices(), 0);
  }

  // Non-clipping - polygons intersect at a point
  {
    Polygon2D subjectPolygon;
    subjectPolygon.addVertex(Point2D {0, 0});
    subjectPolygon.addVertex(Point2D {1, 0});
    subjectPolygon.addVertex(Point2D {1, 1});
    subjectPolygon.addVertex(Point2D {0, 1});

    Polygon2D clipPolygon;
    clipPolygon.addVertex(Point2D {-1, -1});
    clipPolygon.addVertex(Point2D {0, -1});
    clipPolygon.addVertex(Point2D {0, 0});
    clipPolygon.addVertex(Point2D {-1, 0});

    Polygon2D poly = axom::primal::clip(subjectPolygon, clipPolygon, EPS, CHECK_SIGN);

    EXPECT_EQ(poly.isValid(), false);
    EXPECT_EQ(poly.numVertices(), 0);
  }

  // Rosetta Code example
  {
    Polygon2D subjectPolygon;
    subjectPolygon.addVertex(Point2D {50, 150});
    subjectPolygon.addVertex(Point2D {200, 50});
    subjectPolygon.addVertex(Point2D {350, 150});
    subjectPolygon.addVertex(Point2D {350, 300});
    subjectPolygon.addVertex(Point2D {250, 300});
    subjectPolygon.addVertex(Point2D {200, 250});
    subjectPolygon.addVertex(Point2D {150, 350});
    subjectPolygon.addVertex(Point2D {100, 250});
    subjectPolygon.addVertex(Point2D {100, 200});

    Polygon2D clipPolygon;
    clipPolygon.addVertex(Point2D {100, 100});
    clipPolygon.addVertex(Point2D {300, 100});
    clipPolygon.addVertex(Point2D {300, 300});
    clipPolygon.addVertex(Point2D {100, 300});

    Polygon2D poly = axom::primal::clip(subjectPolygon, clipPolygon, EPS);
    EXPECT_NEAR(poly.signedArea(), 37083.3333333333, EPS);
    EXPECT_EQ(poly.numVertices(), 10);

    // Check vertices
    EXPECT_NEAR(poly[0][0], 100, EPS);
    EXPECT_NEAR(poly[0][1], 116.6666666666, EPS);

    EXPECT_NEAR(poly[1][0], 125, EPS);
    EXPECT_NEAR(poly[1][1], 100, EPS);

    EXPECT_NEAR(poly[2][0], 275, EPS);
    EXPECT_NEAR(poly[2][1], 100, EPS);

    EXPECT_NEAR(poly[3][0], 300, EPS);
    EXPECT_NEAR(poly[3][1], 116.6666666666, EPS);

    EXPECT_NEAR(poly[4][0], 300, EPS);
    EXPECT_NEAR(poly[4][1], 300, EPS);

    EXPECT_NEAR(poly[5][0], 250, EPS);
    EXPECT_NEAR(poly[5][1], 300, EPS);

    EXPECT_NEAR(poly[6][0], 200, EPS);
    EXPECT_NEAR(poly[6][1], 250, EPS);

    EXPECT_NEAR(poly[7][0], 175, EPS);
    EXPECT_NEAR(poly[7][1], 300, EPS);

    EXPECT_NEAR(poly[8][0], 125, EPS);
    EXPECT_NEAR(poly[8][1], 300, EPS);

    EXPECT_NEAR(poly[9][0], 100, EPS);
    EXPECT_NEAR(poly[9][1], 250, EPS);
  }

  /*
   * Edge case where polygons overlap along diagonal of one polygon,
   * and the edge of the other:
   *    +-----------+
   *   /| \         |
   *  / |  \        |
   * +  |   +       |
   *  \ |  /        |
   *   \| /         |
   *    +-----------+
   */
  {
    Polygon2D subjectPolygon;
    subjectPolygon.addVertex(Point2D {0., 0.});
    subjectPolygon.addVertex(Point2D {1., 0.});
    subjectPolygon.addVertex(Point2D {1., 1.});
    subjectPolygon.addVertex(Point2D {0., 1.});

    Polygon2D clipPolygon;
    clipPolygon.addVertex(Point2D {0., 0.});
    clipPolygon.addVertex(Point2D {0.3, 0.3});
    clipPolygon.addVertex(Point2D {0., 1.});
    clipPolygon.addVertex(Point2D {-0.3, 0.7});

    Polygon2D poly = axom::primal::clip(subjectPolygon, clipPolygon, EPS);
    EXPECT_NEAR(poly.signedArea(), 0.15, EPS);
    EXPECT_EQ(poly.numVertices(), 3);

    // Check vertices
    EXPECT_NEAR(poly[0][0], 0, EPS);
    EXPECT_NEAR(poly[0][1], 1, EPS);

    EXPECT_NEAR(poly[1][0], 0, EPS);
    EXPECT_NEAR(poly[1][1], 0, EPS);

    EXPECT_NEAR(poly[2][0], 0.3, EPS);
    EXPECT_NEAR(poly[2][1], 0.3, EPS);
  }
}

/*
 * 3x3 grid of quads that occupies the unit square,
 * that clip against a lower or upper triangle:
 *
 * +---+---+---+
 * | \ |   |   |
 * |  \|   |   |
 * +---+---+---+
 * |   | \ |   |
 * |   |  \|   |
 * +---+---+---+
 * |   |   | \ |
 * |   |   |  \|
 * +---+---+---+
 */
TEST(primal_clip, polygon_intersect_robustness)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;
  using Point2D = axom::primal::Point<double, 2>;
  constexpr double EPS = 1e-8;

  const double x[] =
    {0., 1. / 3., 2. / 3, 1., 0., 1. / 3., 2. / 3, 1., 0., 1. / 3., 2. / 3, 1., 0., 1. / 3., 2. / 3, 1.};
  const double y[] =
    {0., 0., 0., 0., 1. / 3., 1. / 3, 1. / 3., 1. / 3., 2. / 3., 2. / 3, 2. / 3., 2. / 3., 1., 1., 1., 1.};
  const int conn[][4] = {{0, 1, 5, 4},
                         {1, 2, 6, 5},
                         {2, 3, 7, 6},
                         {4, 5, 9, 8},
                         {5, 6, 10, 9},
                         {6, 7, 11, 10},
                         {8, 9, 13, 12},
                         {9, 10, 14, 13},
                         {10, 11, 15, 14}};
  Polygon2D fine[9];
  for(int p = 0; p < 9; p++)
  {
    fine[p].addVertex(Point2D {x[conn[p][0]], y[conn[p][0]]});
    fine[p].addVertex(Point2D {x[conn[p][1]], y[conn[p][1]]});
    fine[p].addVertex(Point2D {x[conn[p][2]], y[conn[p][2]]});
    fine[p].addVertex(Point2D {x[conn[p][3]], y[conn[p][3]]});
  }

  // Overlap polygons
  Polygon2D clipLower;
  clipLower.addVertex(Point2D {0., 0.});
  clipLower.addVertex(Point2D {1., 0.});
  clipLower.addVertex(Point2D {0., 1.});

  Polygon2D clipUpper;
  clipUpper.addVertex(Point2D {1., 0.});
  clipUpper.addVertex(Point2D {1., 1.});
  clipUpper.addVertex(Point2D {0., 1.});

  // Expected VFs (volume fractions, areas) for overlaps.
  const double lower_vf[] = {1., 1., 0.5, 1., 0.5, 0., 0.5, 0., 0.};
  const double upper_vf[] = {0., 0., 0.5, 0., 0.5, 1., 0.5, 1., 1.};

  // Each quad has an area of 1/9 (3x3 quads in the unit square).
  constexpr double fine_area = (1. / 3.) * (1. / 3.);

  for(int p : std::vector<int> {0, 1, 2, 3, 4, 6})
  {
    // Clip clipLower with overlapping fine quad
    const auto overlapPoly = axom::primal::clip(clipLower, fine[p], EPS);
    EXPECT_TRUE(overlapPoly.isValid());

    if(overlapPoly.isValid())
    {
      // Area of overlapped polygon is some fraction of 1/9.
      // Divide by 1/9 to get a volume fraction.
      const double vf = overlapPoly.area() / fine_area;
      EXPECT_NEAR(vf, lower_vf[p], EPS);
    }
  }

  for(int p : std::vector<int> {5, 7, 8})
  {
    // Clip clipLower with non-intersecting fine quad
    const auto overlapPoly = axom::primal::clip(clipLower, fine[p], EPS);
    EXPECT_FALSE(overlapPoly.isValid());
    EXPECT_EQ(overlapPoly.numVertices(), 0);
  }

  for(int p : std::vector<int> {2, 4, 5, 6, 7, 8})
  {
    // Clip clipUpper with overlapping fine quad
    const auto overlapPoly = axom::primal::clip(clipUpper, fine[p], EPS);
    EXPECT_TRUE(overlapPoly.isValid());

    if(overlapPoly.isValid())
    {
      // Area of overlapped polygon is some fraction of 1/9.
      // Divide by 1/9 to get a volume fraction.
      const double vf = overlapPoly.area() / fine_area;
      EXPECT_NEAR(vf, upper_vf[p], EPS);
    }
  }

  for(int p : std::vector<int> {0, 1, 3})
  {
    // Clip clipUpper with non-intersecting fine quad
    const auto overlapPoly = axom::primal::clip(clipUpper, fine[p], EPS);
    EXPECT_FALSE(overlapPoly.isValid());
    EXPECT_EQ(overlapPoly.numVertices(), 0);
  }
}

TEST(primal_clip, polygon_clip_regression)
{
  // Use float
  using Polygon2D = axom::primal::Polygon<float, 2, axom::primal::PolygonArray::Dynamic, 5>;
  using Point2D = axom::primal::Point<float, 2>;
  using Vector2D = axom::primal::Vector<float, 2>;
  using Plane2D = axom::primal::Plane<float, 2>;

  // Shape to clip
  Polygon2D poly;
  poly.addVertex(Point2D {148.f, 141.f});
  poly.addVertex(Point2D {149.f, 141.f});
  poly.addVertex(Point2D {149.f, 142.f});
  poly.addVertex(Point2D {148.f, 142.f});

  // Clipping plane
  const float normal[] = {-0.7599482f, -0.6499836f};
  const float pt[] = {148.5418f, 141.5357f};
  constexpr bool normalize = false;
  Plane2D plane(Vector2D(normal, 2), Point2D(pt, 2), normalize);

  // Clip the shape
  auto clippedShape = axom::primal::clip(poly, plane);

  // Expected polygon.
  Polygon2D expectedPoly;
  expectedPoly.addVertex(Point2D {148.f, 141.f});
  expectedPoly.addVertex(Point2D {149.f, 141.f});
  expectedPoly.addVertex(Point2D {148.144683f, 142.f});
  expectedPoly.addVertex(Point2D {148.f, 142.f});

  // Comparisons
  constexpr float EPS = 1.6e-5f;
  EXPECT_EQ(clippedShape.numVertices(), 4);
  for(int i = 0; i < 4; i++)
  {
    EXPECT_NEAR(clippedShape[i][0], expectedPoly[i][0], EPS);
    EXPECT_NEAR(clippedShape[i][1], expectedPoly[i][1], EPS);
  }
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
