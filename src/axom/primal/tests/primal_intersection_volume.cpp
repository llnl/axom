// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal.hpp"

#include "gtest/gtest.h"

#include <algorithm>

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

namespace
{
using Primal3D::HexahedronType;
using Primal3D::TetrahedronType;

double intersection_volume_via_hex_triangulation(const HexahedronType& hex,
                                                 const TetrahedronType& tet,
                                                 double eps)
{
  axom::StackArray<TetrahedronType, 24> tets;
  hex.triangulate(tets);

  double vol = 0.;
  constexpr bool TRY_FIX_ORIENTATION = true;
  for(int i = 0; i < 24; ++i)
  {
    vol += axom::primal::intersection_volume<double>(tets[i], tet, eps, TRY_FIX_ORIENTATION);
  }
  return vol;
}

}  // namespace

TEST(primal_intersection_volume, hex_tet_user_regression_cases)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-10;
  constexpr double sixth = 1. / 6.;

  HexahedronType hex(PointType {42. + sixth, -66, -178.5},
                     PointType {52. - sixth, -66, -178.5},
                     PointType {52. - sixth, -55, -178.5},
                     PointType {42. + sixth, -55, -178.5},
                     PointType {42. + sixth, -66, -170},
                     PointType {52. - sixth, -66, -170},
                     PointType {52. - sixth, -55, -170},
                     PointType {42. + sixth, -55, -170});

  const TetrahedronType cases[] = {
    TetrahedronType {PointType {27.5859, -19.5363, -148.01},
                     PointType {44.2539, -58.1624, -171.152},
                     PointType {43.7957, -57.9539, -146.494},
                     PointType {15.8564, -32.9522, -147.246}},
    TetrahedronType {PointType {76.8265, -45.3561, -146.396},
                     PointType {79.5055, -43.6324, -171.152},
                     PointType {77.7366, -63.3987, -171.152},
                     PointType {44.2539, -58.1624, -171.152}},
    TetrahedronType {PointType {77.4836, -63.4271, -145.519},
                     PointType {76.8265, -45.3561, -146.396},
                     PointType {44.2539, -58.1624, -171.152},
                     PointType {43.7957, -57.9539, -146.494}},
    TetrahedronType {PointType {77.4836, -63.4271, -145.519},
                     PointType {76.8265, -45.3561, -146.396},
                     PointType {77.7366, -63.3987, -171.152},
                     PointType {44.2539, -58.1624, -171.152}},
  };

  const double expected_volumes[] = {0.30774371225316693,
                                     15.302276033131164,
                                     0.2431953788553831,
                                     3.1668745904354338};

  for(int i = 0; i < 4; ++i)
  {
    const auto& tet = cases[i];
    const TetrahedronType tet_flipped(tet[1], tet[0], tet[2], tet[3]);

    constexpr bool fix_orient = true;
    const double direct_hex_tet =
      axom::primal::intersection_volume<double>(hex, tet, EPS, fix_orient);
    const double ref_subdiv_hex_tet = intersection_volume_via_hex_triangulation(hex, tet, EPS);

    // check that intersection volumes are non-negative
    EXPECT_GT(direct_hex_tet, 0.) << "case " << i;
    EXPECT_GT(ref_subdiv_hex_tet, 0.) << "case " << i;

    // check that primal::intersection_volume results matches expectations and other ways to compute it
    const double expected = expected_volumes[i];
    EXPECT_NEAR(direct_hex_tet, expected, EPS) << "case " << i;
    EXPECT_NEAR(direct_hex_tet, ref_subdiv_hex_tet, EPS) << "case " << i;
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
