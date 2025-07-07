// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file hex_tet_volume.cpp
 *
 * Example that demonstrates use of Primal's intersection_volume operator to find
 * the volume of intersection between hexahedra and tetrahedra. Supports host
 * and device execution using RAJA.
 *
 * \note This example requires RAJA and Umpire.
 */

#include "axom/primal.hpp"
#include "axom/slic.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <cmath>
#include <map>
#include <string>

using HexahedronType = typename axom::primal::Hexahedron<double>;
using PointType = typename axom::primal::Point<double, 3>;
using TetrahedronType = typename axom::primal::Tetrahedron<double, 3>;

/// Choose runtime policy for RAJA
enum class RuntimePolicy
{
  seq = 0,
#ifdef AXOM_USE_OPENMP
  omp = 1,
#endif
#ifdef AXOM_USE_CUDA
  cuda = 2,
#endif
#ifdef AXOM_USE_HIP
  hip = 3
#endif
};

struct Input
{
public:
  int hexResolution {5};
  int tetResolution {5};
  RuntimePolicy policy {RuntimePolicy::seq};

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    std::string pol_str = "Sets execution space of intersection_volume operator.";
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    pol_str += "\nSet to \'seq\' or 0 to use the RAJA sequential policy.";
  #ifdef AXOM_USE_OPENMP
    pol_str += "\nSet to \'omp\' or 1 to use the RAJA OpenMP policy.";
  #endif
  #ifdef AXOM_USE_CUDA
    pol_str += "\nSet to \'cuda\' or 2 to use the RAJA CUDA policy.";
  #endif
  #ifdef AXOM_USE_HIP
    pol_str += "\nSet to \'hip\' or 3 to use the RAJA HIP policy.";
  #endif
#endif

    app.add_option("-p, --policy", this->policy, pol_str)
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(validExecPolicies));

    app.add_option("-x,--hex-resolution", hexResolution)
      ->description("Number of hexahedra to generate")
      ->capture_default_str()
      ->check(axom::CLI::PositiveNumber);

    app.add_option("-t,--tet-resolution", tetResolution)
      ->description("Number of tetrahedra to generate")
      ->capture_default_str()
      ->check(axom::CLI::PositiveNumber);

    app.get_formatter()->column_width(80);

    app.parse(argc, argv);
  }

private:
  // clang-format off
  const std::map<std::string, RuntimePolicy> validExecPolicies{
      {"seq", RuntimePolicy::seq}
    #if defined(AXOM_USE_RAJA)
      #ifdef AXOM_USE_OPENMP
    , {"omp", RuntimePolicy::omp}
      #endif
      #ifdef AXOM_USE_CUDA
    , {"cuda", RuntimePolicy::cuda}
      #endif
      #ifdef AXOM_USE_HIP
    , {"hip", RuntimePolicy::hip}
      #endif
    #endif
  };
  // clang-format on
};

// Helper function to generate a hexahedron cube with given corner point and
// edge length
HexahedronType generateCube(const PointType& point, double length)
{
  return HexahedronType(point,
                        PointType {point[0], point[1], point[2] + length},
                        PointType {point[0] + length, point[1], point[2] + length},
                        PointType {point[0] + length, point[1], point[2]},
                        PointType {point[0], point[1] + length, point[2]},
                        PointType {point[0], point[1] + length, point[2] + length},
                        PointType {point[0] + length, point[1] + length, point[2] + length},
                        PointType {point[0] + length, point[1] + length, point[2]});
}

// Function to check intersection volumes of generated hexahedra and tetrahedra
template <typename ExecSpace>
void check_intersection_volumes(const Input& params)
{
  SLIC_INFO(
    axom::fmt::format("{:-^80}",
                      axom::fmt::format("Running intersection volume check in execution space: {}",
                                        axom::execution_space<ExecSpace>::name())));

  // Get allocators
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecSpace>::allocatorID();

  // Generate hexahedra subdividing the unit cube with corner points
  // (-1,-1,-1) and (1,1,1)
  int const HEX_RESOLUTION = params.hexResolution;
  int hex_index = 0;
  int const NUM_HEXES = HEX_RESOLUTION * HEX_RESOLUTION * HEX_RESOLUTION;
  axom::Array<HexahedronType> hexes_h(NUM_HEXES, NUM_HEXES, host_allocator);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Generating {} hexahedra with hexahedra resolution set to {}",
                      NUM_HEXES,
                      HEX_RESOLUTION)));

  SLIC_INFO(axom::fmt::format("{: ^80}",
                              "Hexahedra subdivide the unit cube with corner points (-1,-1,-1) "
                              "and (1,1,1)"));

  for(int i = 0; i < HEX_RESOLUTION; i++)
  {
    for(int j = 0; j < HEX_RESOLUTION; j++)
    {
      for(int k = 0; k < HEX_RESOLUTION; k++)
      {
        double edge_length = 2.0 / HEX_RESOLUTION;
        hexes_h[hex_index] =
          generateCube(PointType {edge_length * i - 1, edge_length * j - 1, edge_length * k - 1},
                       edge_length);
        hex_index++;
      }
    }
  }

  // Generate tetrahedra from unit sphere with center (0,0,0)
  int const TET_RESOLUTION = params.tetResolution;
  int tet_index = 0;
  int const NUM_TETS = 4 * std::pow(2, TET_RESOLUTION);
  axom::Array<TetrahedronType> tets_h(NUM_TETS, NUM_TETS, host_allocator);
  double step_size = 1.0 / std::pow(2, TET_RESOLUTION);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Generating {} tetrahedra with tetrahedra resolution set to {}",
                      NUM_TETS,
                      TET_RESOLUTION)));

  SLIC_INFO(axom::fmt::format(
    "{: ^80}",
    axom::fmt::format("Tetrahedra are encapsulated by the unit sphere at the origin")));

  for(int i = 0; i < std::pow(2, TET_RESOLUTION); i++)
  {
    for(int j = 0; j <= 1; j++)
    {
      for(int k = 0; k <= 1; k++)
      {
        // Tetrahedron coordinates consist of:
        //   - two points on the unit sphere great circle where y = 0
        //   - the origin
        //   - one of the "poles" of the sphere, either (0,1,0) or (0, -1, 0)
        // The generated tetrahedra are all encapsulated by the unit cube.

        double pole_sign = j ? 1 : -1;
        double z_sign = k ? 1 : -1;

        double x1 = axom::utilities::lerp<double>(-1.0, 1.0, step_size * i);
        double z1 = std::sqrt(1 - (x1 * x1)) * z_sign;
        double x2 = axom::utilities::lerp<double>(-1.0, 1.0, step_size * (i + 1));
        double z2 = std::sqrt(1 - (x2 * x2)) * z_sign;

        tets_h[tet_index] = TetrahedronType(PointType {x1, 0, z1},
                                            PointType::zero(),
                                            PointType {x2, 0, z2},
                                            PointType {0, pole_sign, 0});
        tet_index++;
      }
    }
  }

  // Calculate expected sum of all tetrahedra volume
  axom::Array<HexahedronType> hexes_d = on_device
    ? axom::Array<HexahedronType>(hexes_h, device_allocator)
    : axom::Array<HexahedronType>();
  auto hexes_view = on_device ? hexes_d.view() : hexes_h.view();

  axom::Array<TetrahedronType> tets_d = on_device
    ? axom::Array<TetrahedronType>(tets_h, device_allocator)
    : axom::Array<TetrahedronType>();
  auto tets_view = on_device ? tets_d.view() : tets_h.view();

  axom::ReduceSum<ExecSpace, double> total_tet_vol(0.0);

  axom::for_all<ExecSpace>(
    NUM_TETS,
    AXOM_LAMBDA(axom::IndexType i) { total_tet_vol += tets_view[i].volume(); });

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Total volume of all tetrahedra is {} ", total_tet_vol.get())));

  // Calculate intersection volume for each hexahedra and tetrahedra pair.
  // Typically, a spatial index (e.g. Bounding Volume Hierarchy) can be used to
  // reduce the number of operations.
  axom::ReduceSum<ExecSpace, double> total_intersect_vol(0.0);
  constexpr double EPS = 1e-10;
  constexpr bool tryFixOrientation = true;

  // The lower of the two sizes (NUM_HEXES, NUM_TETS) is used to factor out
  // every pair of hexahedron and tetrahedron indices.
  if(NUM_HEXES > NUM_TETS)
  {
    axom::for_all<ExecSpace>(
      NUM_HEXES * NUM_TETS,
      AXOM_LAMBDA(axom::IndexType i) {
        total_intersect_vol += intersection_volume(hexes_view[i / NUM_TETS],
                                                   tets_view[i % NUM_TETS],
                                                   EPS,
                                                   tryFixOrientation);
      });
  }
  else
  {
    axom::for_all<ExecSpace>(
      NUM_HEXES * NUM_TETS,
      AXOM_LAMBDA(axom::IndexType i) {
        total_intersect_vol += intersection_volume(hexes_view[i % NUM_HEXES],
                                                   tets_view[i / NUM_HEXES],
                                                   EPS,
                                                   tryFixOrientation);
      });
  }

  SLIC_INFO(axom::fmt::format("{:-^80}",
                              axom::fmt::format("Total intersect volume between all hexahedra "
                                                "and tetrahedra is {} ",
                                                total_intersect_vol.get())));

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Difference between sums is {}",
                      std::abs(total_intersect_vol.get() - total_tet_vol.get()))));
}

int main(int argc, char** argv)
{
  // Initialize the SLIC logger
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  // Set up and parse command line arguments
  Input params;
  axom::CLI::App app {"Example of intersection volume between hexahedra and tetrahedra"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  switch(params.policy)
  {
  case RuntimePolicy::seq:
    check_intersection_volumes<axom::SEQ_EXEC>(params);
    break;
#if defined(AXOM_USE_RAJA)
  #ifdef AXOM_USE_OPENMP
  case RuntimePolicy::omp:
    check_intersection_volumes<axom::OMP_EXEC>(params);
    break;
  #endif
  #ifdef AXOM_USE_CUDA
  case RuntimePolicy::cuda:
    check_intersection_volumes<axom::CUDA_EXEC<256>>(params);
    break;
  #endif
  #ifdef AXOM_USE_HIP
  case RuntimePolicy::hip:
    check_intersection_volumes<axom::HIP_EXEC<256>>(params);
    break;
  #endif
#endif
  }
  return 0;
}
