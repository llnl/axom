// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/primal.hpp"

#include "axom/fmt.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <random>

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
using SphereType = axom::primal::Sphere<double, 3>;
}  // namespace Primal3D

namespace
{
using Primal3D::HexahedronType;
using Primal3D::PointType;
using Primal3D::SphereType;
using Primal3D::TetrahedronType;

double intersection_volume_via_hex_triangulation(const HexahedronType& hex,
                                                 const TetrahedronType& tet,
                                                 double eps)
{
  axom::StackArray<TetrahedronType, HexahedronType::NUM_TRIANGULATE> tets;
  hex.triangulate(tets);

  double vol = 0.;
  constexpr bool TRY_FIX_ORIENTATION = true;
  for(int i = 0; i < HexahedronType::NUM_TRIANGULATE; ++i)
  {
    vol += axom::primal::intersection_volume<double>(tets[i], tet, eps, TRY_FIX_ORIENTATION);
  }
  return vol;
}

double intersection_volume_via_hex_triangulation(const SphereType& sphere, const HexahedronType& hex)
{
  axom::StackArray<TetrahedronType, HexahedronType::NUM_TRIANGULATE> tets;
  hex.triangulate(tets);

  double volume = 0.;
  for(const TetrahedronType& tet : tets)
  {
    volume += axom::primal::intersection_volume(sphere, tet);
  }
  return volume;
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
template <typename ExecPolicy>
void check_sphere_tet_intersection_volume_on_device()
{
  constexpr int num_cases = 5;
  const SphereType unit_sphere(PointType::zero(), 1.);
  const std::array<SphereType, num_cases> spheres {
    {unit_sphere, unit_sphere, unit_sphere, unit_sphere, unit_sphere}};
  const std::array<TetrahedronType, num_cases> tets {{TetrahedronType(PointType {0., 0., 0.},
                                                                      PointType {0.25, 0., 0.},
                                                                      PointType {0., 0.25, 0.},
                                                                      PointType {0., 0., 0.25}),
                                                      TetrahedronType(PointType {3., 3., 3.},
                                                                      PointType {3., -3., -3.},
                                                                      PointType {-3., 3., -3.},
                                                                      PointType {-3., -3., 3.}),
                                                      TetrahedronType(PointType {2., 2., 2.},
                                                                      PointType {3., 2., 2.},
                                                                      PointType {2., 3., 2.},
                                                                      PointType {2., 2., 3.}),
                                                      TetrahedronType(PointType {0., 0., 0.},
                                                                      PointType {3., 0., 0.},
                                                                      PointType {0., 3., 0.},
                                                                      PointType {0., 0., 3.}),
                                                      TetrahedronType(PointType {1., 0., 0.},
                                                                      PointType {2., 1., 0.},
                                                                      PointType {2., 0., 1.},
                                                                      PointType {2., -1., -1.})}};
  const std::array<double, num_cases> expected {{tets[0].volume(), 4. * M_PI / 3., 0., M_PI / 6., 0.}};

  const int allocator_id = axom::execution_space<ExecPolicy>::allocatorID();
  SphereType* device_spheres = axom::allocate<SphereType>(num_cases, allocator_id);
  TetrahedronType* device_tets = axom::allocate<TetrahedronType>(num_cases, allocator_id);
  double* device_results = axom::allocate<double>(2 * num_cases, allocator_id);
  axom::copy(device_spheres, spheres.data(), sizeof(SphereType) * num_cases);
  axom::copy(device_tets, tets.data(), sizeof(TetrahedronType) * num_cases);

  axom::for_all<ExecPolicy>(num_cases, [=] AXOM_HOST_DEVICE(axom::IndexType i) {
    device_results[i] = axom::primal::intersection_volume(device_spheres[i], device_tets[i]);
    device_results[num_cases + i] =
      axom::primal::intersection_volume(device_tets[i], device_spheres[i]);
  });

  std::array<double, 2 * num_cases> results;
  axom::copy(results.data(), device_results, sizeof(double) * results.size());
  for(int i = 0; i < num_cases; ++i)
  {
    EXPECT_NEAR(results[i], expected[i], 1.e-12) << "case " << i;
    EXPECT_NEAR(results[num_cases + i], expected[i], 1.e-12) << "reverse case " << i;
  }

  axom::deallocate(device_results);
  axom::deallocate(device_tets);
  axom::deallocate(device_spheres);
}
#endif

}  // namespace

TEST(primal_intersection_volume, hex_tet_user_regression_cases)
{
  using namespace Primal3D;
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

  // Use intersection volumes computed using exact arithmetic in cgal
  const double expected_volumes[] = {0.30774371225334463,
                                     15.302276033130795,
                                     0.24319537885473608,
                                     3.1668745904352349};

  constexpr bool DO_NOT_FIX_ORIENTATION = false;
  constexpr bool TRY_FIX_ORIENTATION = true;

  // clang-format off
  constexpr bool plot_all_values = false;
  const auto eps_values = plot_all_values
    ? std::vector<double> {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8,
                           1e-9, 1e-10,1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16}
    : std::vector<double> {1e-10};
  // clang-format on

  // print the header
  axom::fmt::print("case, eps, expected, direct, direct_fixed, direct_fixed_flipped, ref_subdiv\n");

  for(int i = 0; i < 4; ++i)
  {
    const auto& tet = cases[i];
    const TetrahedronType tet_flipped(tet[1], tet[0], tet[2], tet[3]);

    for(const double eps : eps_values)
    {
      const double direct_hex_tet =
        axom::primal::intersection_volume<double>(hex, tet, eps, DO_NOT_FIX_ORIENTATION);
      const double direct_hex_tet_fixed =
        axom::primal::intersection_volume<double>(hex, tet, eps, TRY_FIX_ORIENTATION);
      const double direct_hex_tet_fixed_flipped =
        axom::primal::intersection_volume<double>(hex, tet_flipped, eps, TRY_FIX_ORIENTATION);
      const double ref_subdiv_hex_tet = intersection_volume_via_hex_triangulation(hex, tet, eps);

      axom::fmt::print("{},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}\n",
                       i,
                       eps,
                       expected_volumes[i],
                       direct_hex_tet,
                       direct_hex_tet_fixed,
                       direct_hex_tet_fixed_flipped,
                       ref_subdiv_hex_tet);

      // check that intersection volumes are non-negative
      EXPECT_GT(direct_hex_tet, 0.) << "case " << i;
      EXPECT_GT(ref_subdiv_hex_tet, 0.) << "case " << i;
      EXPECT_GT(direct_hex_tet_fixed_flipped, 0.) << "case " << i;
      EXPECT_GT(ref_subdiv_hex_tet, 0.) << "case " << i;

      // check that primal::intersection_volume results matches expectations and other ways to compute it
      // we expect about 12 digits of precision to baseline values
      if(eps > 1e-13)
      {
        EXPECT_NEAR(direct_hex_tet_fixed, expected_volumes[i], eps) << "case " << i;
        EXPECT_NEAR(direct_hex_tet_fixed, ref_subdiv_hex_tet, eps) << "case " << i;
        EXPECT_NEAR(direct_hex_tet_fixed, direct_hex_tet_fixed_flipped, eps) << "case " << i;
        EXPECT_NEAR(direct_hex_tet_fixed, ref_subdiv_hex_tet, eps) << "case " << i;
      }
    }
  }
}

TEST(primal_intersection_volume, sphere_tet_trivial_and_analytic_cases)
{
  using namespace Primal3D;

  const SphereType unit_sphere(PointType::zero(), 1.0);
  const TetrahedronType contained_tet(PointType {0., 0., 0.},
                                      PointType {0.25, 0., 0.},
                                      PointType {0., 0.25, 0.},
                                      PointType {0., 0., 0.25});
  EXPECT_DOUBLE_EQ(axom::primal::intersection_volume(unit_sphere, contained_tet),
                   contained_tet.volume());

  const TetrahedronType containing_tet(PointType {3., 3., 3.},
                                       PointType {3., -3., -3.},
                                       PointType {-3., 3., -3.},
                                       PointType {-3., -3., 3.});
  EXPECT_NEAR(axom::primal::intersection_volume(unit_sphere, containing_tet), 4. * M_PI / 3., 1.e-14);

  const TetrahedronType disjoint_tet(PointType {2., 2., 2.},
                                     PointType {3., 2., 2.},
                                     PointType {2., 3., 2.},
                                     PointType {2., 2., 3.});
  EXPECT_DOUBLE_EQ(axom::primal::intersection_volume(unit_sphere, disjoint_tet), 0.);

  // This tetrahedron contains exactly the positive octant of the unit sphere.
  const TetrahedronType octant_tet(PointType {0., 0., 0.},
                                   PointType {3., 0., 0.},
                                   PointType {0., 3., 0.},
                                   PointType {0., 0., 3.});
  EXPECT_NEAR(axom::primal::intersection_volume(unit_sphere, octant_tet), M_PI / 6., 1.e-14);

  const TetrahedronType tangent_tet(PointType {1., 0., 0.},
                                    PointType {2., 1., 0.},
                                    PointType {2., 0., 1.},
                                    PointType {2., -1., -1.});
  EXPECT_NEAR(axom::primal::intersection_volume(unit_sphere, tangent_tet), 0., 1.e-14);

  const TetrahedronType degenerate_tet(PointType {0., 0., 0.},
                                       PointType {1., 0., 0.},
                                       PointType {0., 1., 0.},
                                       PointType {1., 1., 0.});
  EXPECT_DOUBLE_EQ(axom::primal::intersection_volume(unit_sphere, degenerate_tet), 0.);
  EXPECT_DOUBLE_EQ(axom::primal::intersection_volume(SphereType(PointType::zero(), 0.), octant_tet),
                   0.);

  using Point3F = axom::primal::Point<float, 3>;
  const axom::primal::Sphere<float, 3> float_sphere(Point3F::zero(), 1.f);
  const axom::primal::Tetrahedron<float, 3> float_octant(Point3F {0.f, 0.f, 0.f},
                                                         Point3F {3.f, 0.f, 0.f},
                                                         Point3F {0.f, 3.f, 0.f},
                                                         Point3F {0.f, 0.f, 3.f});
  EXPECT_NEAR(axom::primal::intersection_volume(float_sphere, float_octant),
              static_cast<float>(M_PI / 6.),
              1.e-6f);
}

TEST(primal_intersection_volume, sphere_tet_hex_decomposition_analytic_cases)
{
  using namespace Primal3D;

  const HexahedronType unit_hex(PointType {-1., -1., -1.},
                                PointType {1., -1., -1.},
                                PointType {1., 1., -1.},
                                PointType {-1., 1., -1.},
                                PointType {-1., -1., 1.},
                                PointType {1., -1., 1.},
                                PointType {1., 1., 1.},
                                PointType {-1., 1., 1.});
  const double unit_sphere_volume = 4. * M_PI / 3.;

  struct TestCase
  {
    const char* name;
    SphereType sphere;
    double expected;
  };

  // These symmetry cases exercise face-cap, edge-wedge, and vertex-cone
  // contributions through the public sphere-tetrahedron operator.
  const std::array<TestCase, 8> cases {{
    {"face tangent", SphereType(PointType {0., 2., 0.}, 1.), 0.},
    {"face centered", SphereType(PointType {1., 0., 0.}, 1.), unit_sphere_volume / 2.},
    {"edge centered", SphereType(PointType {0., -1., 1.}, 1.), unit_sphere_volume / 4.},
    {"vertex centered", SphereType(PointType {1., -1., 1.}, 1.), unit_sphere_volume / 8.},
    {"edge tangent", SphereType(PointType {1. + std::sqrt(0.5), 1. + std::sqrt(0.5), 0.}, 1.), 0.},
    {"vertex tangent", SphereType(PointType {2., -1., 1.}, 1.), 0.},
    {"sphere contained", SphereType(PointType::zero(), 0.5), M_PI / 6.},
    {"hex contained", SphereType(PointType::zero(), 2.), unit_hex.volume()},
  }};

  for(const TestCase& test : cases)
  {
    const double volume = intersection_volume_via_hex_triangulation(test.sphere, unit_hex);
    EXPECT_NEAR(volume, test.expected, 1.e-12) << test.name;
  }

  // The oblique prism intersects the unit sphere in a pi/4 spherical wedge.
  const HexahedronType wedge_hex(PointType {0., 0., -2.},
                                 PointType {3., 3., -2.},
                                 PointType {3., 5., -2.},
                                 PointType {0., 5., -2.},
                                 PointType {0., 0., 2.},
                                 PointType {3., 3., 2.},
                                 PointType {3., 5., 2.},
                                 PointType {0., 5., 2.});
  const SphereType unit_sphere(PointType::zero(), 1.);
  EXPECT_NEAR(intersection_volume_via_hex_triangulation(unit_sphere, wedge_hex), M_PI / 6., 1.e-12);
}

TEST(primal_intersection_volume, sphere_tet_orientation_and_scale_invariance)
{
  using namespace Primal3D;

  const std::array<PointType, 4> points = {{PointType {-0.8, -0.7, -0.2},
                                            PointType {1.3, -0.4, 0.1},
                                            PointType {-0.2, 1.2, 0.3},
                                            PointType {0.1, 0.2, 1.4}}};
  const SphereType sphere(PointType {0.1, -0.05, 0.15}, 0.9);
  const double expected =
    axom::primal::intersection_volume(sphere,
                                      TetrahedronType(points[0], points[1], points[2], points[3]));

  std::array<int, 4> permutation = {{0, 1, 2, 3}};
  do
  {
    const TetrahedronType tet(points[permutation[0]],
                              points[permutation[1]],
                              points[permutation[2]],
                              points[permutation[3]]);
    EXPECT_NEAR(axom::primal::intersection_volume(sphere, tet), expected, 2.e-14);
    EXPECT_NEAR(axom::primal::intersection_volume(tet, sphere), expected, 2.e-14);
  } while(std::next_permutation(permutation.begin(), permutation.end()));

  for(const double scale : {1.e-6, 1., 1.e6})
  {
    const PointType center {0.1 * scale, -0.05 * scale, 0.15 * scale};
    std::array<PointType, 4> scaled_points;
    for(int vertex = 0; vertex < 4; ++vertex)
    {
      scaled_points[vertex] =
        PointType {points[vertex][0] * scale, points[vertex][1] * scale, points[vertex][2] * scale};
    }
    const SphereType scaled_sphere(center, 0.9 * scale);
    const TetrahedronType scaled_tet(scaled_points[0],
                                     scaled_points[1],
                                     scaled_points[2],
                                     scaled_points[3]);
    const double scaled_volume = axom::primal::intersection_volume(scaled_sphere, scaled_tet);
    const double expected_scaled = expected * scale * scale * scale;
    EXPECT_NEAR(scaled_volume, expected_scaled, 5.e-13 * axom::utilities::max(1., expected_scaled));
  }
}

TEST(primal_intersection_volume, sphere_tet_decomposition_invariance)
{
  using namespace Primal3D;

  const SphereType sphere(PointType {-0.01725, 0., 0.}, 1.0);
  const std::array<PointType, 4> points = {
    {PointType {0.9667906976744187, 0., 3.098296812907414e-16},
     PointType {1.002654107311333, 0.0384643285369352, -2.82302142880589e-16},
     PointType {1.002573643410853, 0., 4.131062417209885e-16},
     PointType {1.002573643410853, 0., -0.05063534883720874}}};
  const TetrahedronType tet(points[0], points[1], points[2], points[3]);
  const PointType centroid = PointType::midpoint(PointType::midpoint(points[0], points[1]),
                                                 PointType::midpoint(points[2], points[3]));

  double decomposed_volume = 0.;
  decomposed_volume +=
    axom::primal::intersection_volume(sphere,
                                      TetrahedronType(centroid, points[1], points[2], points[3]));
  decomposed_volume +=
    axom::primal::intersection_volume(sphere,
                                      TetrahedronType(points[0], centroid, points[2], points[3]));
  decomposed_volume +=
    axom::primal::intersection_volume(sphere,
                                      TetrahedronType(points[0], points[1], centroid, points[3]));
  decomposed_volume +=
    axom::primal::intersection_volume(sphere,
                                      TetrahedronType(points[0], points[1], points[2], centroid));

  const double volume = axom::primal::intersection_volume(sphere, tet);
  EXPECT_GT(volume, 0.);
  EXPECT_NEAR(volume, decomposed_volume, 1.e-14);
}

TEST(primal_intersection_volume, sphere_tet_small_scale_regression)
{
  using namespace Primal3D;

  // The sphere center lies close to tet faces at a scale where direct world-
  // space predicates lose precision. Unit-sphere normalization preserves the
  // half-sphere intersection and identifies tangency without spurious volume.
  const SphereType sphere(PointType {5.009999999999999e-7, 5.2e-7, 5.e-7}, 20.e-9);
  const TetrahedronType inside_tet(PointType {5.e-7, 1.e-6, 5.e-7},
                                   PointType {1.e-6, 5.e-7, 5.e-7},
                                   PointType {0., 5.e-7, 5.e-7},
                                   PointType {5.e-7, 5.e-7, 0.});
  const TetrahedronType outside_tet(PointType {0., 5.e-7, 5.e-7},
                                    PointType {5.e-7, 0., 5.e-7},
                                    PointType {1.e-6, 5.e-7, 5.e-7},
                                    PointType {5.e-7, 5.e-7, 1.e-6});
  const double half_sphere_volume = 2. * M_PI * std::pow(sphere.getRadius(), 3) / 3.;
  const double tolerance = 1.e-14 * half_sphere_volume;
  EXPECT_NEAR(axom::primal::intersection_volume(sphere, inside_tet), half_sphere_volume, tolerance);
  EXPECT_DOUBLE_EQ(axom::primal::intersection_volume(sphere, outside_tet), 0.);
}

TEST(primal_intersection_volume, sphere_tet_random_decomposition_invariance)
{
  using namespace Primal3D;

  std::mt19937 generator(928734u);
  std::uniform_real_distribution<double> coordinate(-2., 2.);
  std::uniform_real_distribution<double> radius(0.05, 1.5);
  for(int sample = 0; sample < 1000; ++sample)
  {
    std::array<PointType, 4> points;
    for(PointType& point : points)
    {
      point = PointType {coordinate(generator), coordinate(generator), coordinate(generator)};
    }
    const TetrahedronType tet(points[0], points[1], points[2], points[3]);
    if(tet.volume() < 1.e-10)
    {
      --sample;
      continue;
    }

    const SphereType sphere(
      PointType {coordinate(generator), coordinate(generator), coordinate(generator)},
      radius(generator));
    const PointType centroid = PointType::midpoint(PointType::midpoint(points[0], points[1]),
                                                   PointType::midpoint(points[2], points[3]));
    const std::array<TetrahedronType, 4> subtets = {
      {TetrahedronType(centroid, points[1], points[2], points[3]),
       TetrahedronType(points[0], centroid, points[2], points[3]),
       TetrahedronType(points[0], points[1], centroid, points[3]),
       TetrahedronType(points[0], points[1], points[2], centroid)}};

    const double volume = axom::primal::intersection_volume(sphere, tet);
    double decomposed_volume = 0.;
    for(const TetrahedronType& subtet : subtets)
    {
      decomposed_volume += axom::primal::intersection_volume(sphere, subtet);
    }
    EXPECT_NEAR(volume, decomposed_volume, 2.e-10)
      << "sample " << sample << ", radius " << sphere.getRadius();
  }
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_CUDA)
AXOM_CUDA_TEST(primal_intersection_volume, sphere_tet_cuda)
{
  check_sphere_tet_intersection_volume_on_device<axom::CUDA_EXEC<256>>();
}
  #endif

  #if defined(AXOM_USE_HIP)
TEST(primal_intersection_volume, sphere_tet_hip)
{
  check_sphere_tet_intersection_volume_on_device<axom::HIP_EXEC<256>>();
}
  #endif
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
