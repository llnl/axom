// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_mesh_clipper_benchmark.cpp
 * \brief End-to-end MeshClipper benchmarks over representative geometry types.
 *
 * Background mesh and geometry construction are excluded from timing. Each
 * iteration constructs a MeshClipper and its output, matching the one-clip-per-
 * geometry lifecycle of the shaping path. Device work is synchronized before
 * an iteration ends. After timing, an independent clip checks the result
 * against a reference volume and reports per-clip work counters.
 */

#include "axom/config.hpp"

#include "axom/CLI11.hpp"
#include "axom/core.hpp"
#include "axom/klee.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest/MeshClipper.hpp"
#include "axom/quest/ShapeMesh.hpp"
#include "axom/quest/io/ProEReader.hpp"
#include "axom/quest/util/mesh_helpers.hpp"
#include "axom/quest/util/make_clipper_strategy.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"

#include "benchmark/benchmark.h"

#include <math.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace sidre = axom::sidre;
namespace slic = axom::slic;

namespace
{

using RuntimePolicy = axom::runtime_policy::Policy;
using Point3D = primal::Point<double, 3>;
using Vector3D = primal::Vector<double, 3>;
using Tet3D = primal::Tetrahedron<double, 3>;
using Hex3D = primal::Hexahedron<double, 3>;

constexpr const char* TOPO_NAME = "mesh";
constexpr const char* COORDSET_NAME = "coords";

// These dimensions mirror the default fixtures in quest_mesh_clipper.cpp.
// Those fixtures are compact enough to stay inside the [-2, 2]^3 test domain
// under rotation, and the analytic solids have volumes near 4.2. Keeping their
// sizes comparable makes differences in clipping strategy more visible than
// differences caused simply by one shape occupying much more of the mesh.
constexpr double BOX_HALF_WIDTH = 2.0;
constexpr double SPHERE_RADIUS = 1.0;
constexpr double CYLINDER_RADIUS = 0.695;
constexpr double CYLINDER_HEIGHT = 2.78;
constexpr double CONE_BASE_RADIUS = 1.23;
constexpr double CONE_TOP_RADIUS = 0.176;
constexpr double CONE_HEIGHT = 2.3;
constexpr double TET_LENGTH = 1.55;
constexpr double HEX_MEDIUM_HALF_LENGTH = 0.82;
constexpr double TET_MESH_LENGTH = 1.17;

// The established performance workload uses 51^3, or 132,651, background
// cells. Retaining that size keeps new measurements comparable while providing
// enough work for stable timings. Level 5 keeps curved-shape discretization
// error below the benchmark's validation threshold.
constexpr int DEFAULT_RESOLUTION = 51;
constexpr int DEFAULT_REFINEMENT_LEVEL = 5;
constexpr double MAX_RELATIVE_ERROR = 0.0015;
constexpr klee::TransformableGeometryProperties GEOMETRY_PROPERTIES {klee::Dimensions::Three,
                                                                     klee::LengthUnit::unspecified};

bool s_benchmark_failed = false;

std::shared_ptr<klee::CompositeOperator> makeIdentityOperator()
{
  return std::make_shared<klee::CompositeOperator>(GEOMETRY_PROPERTIES);
}

Tet3D makeTetFixture()
{
  const Point3D a {Point3D::NumericArray {0.8, 0.0, -1.0} * TET_LENGTH};
  const Point3D b {Point3D::NumericArray {-0.8, 1.0, -1.0} * TET_LENGTH};
  const Point3D c {Point3D::NumericArray {-0.8, -1.0, -1.0} * TET_LENGTH};
  const Point3D d {Point3D::NumericArray {0.0, 0.0, 1.0} * TET_LENGTH};
  return Tet3D {a, b, c, d};
}

Hex3D makeHexFixture()
{
  // Unequal side lengths prevent this fixture from reducing to a cube while
  // preserving the same approximate volume as the other analytic solids.
  constexpr double medium = HEX_MEDIUM_HALF_LENGTH;
  constexpr double large = 1.2 * HEX_MEDIUM_HALF_LENGTH;
  constexpr double small = 0.8 * HEX_MEDIUM_HALF_LENGTH;
  return Hex3D {Point3D {-large, -medium, -small},
                Point3D {large, -medium, -small},
                Point3D {large, medium, -small},
                Point3D {-large, medium, -small},
                Point3D {-large, -medium, small},
                Point3D {large, -medium, small},
                Point3D {large, medium, small},
                Point3D {-large, medium, small}};
}

axom::Array<double, 2> makeSorProfile(const std::string& shape)
{
  constexpr double sor_profile[][2] = {{-1.2, 1.1},
                                       {0.48, 1.1},
                                       {0.48, 0.77},
                                       {1.2, 0.77},
                                       {1.2, 0.44},
                                       {0.6, 0.44},
                                       {0.6, 0.33},
                                       {0.0, 0.33},
                                       {0.0, 0.55},
                                       {0.24, 0.55},
                                       {0.24, 0.77},
                                       {-1.2, 0.77}};
  constexpr int sor_point_count = sizeof(sor_profile) / sizeof(sor_profile[0]);
  const bool is_general_sor = shape == "sor";
  const int point_count = is_general_sor ? sor_point_count : 2;
  axom::Array<double, 2> profile({point_count, 2}, axom::ArrayStrideOrder::ROW);

  if(is_general_sor)
  {
    for(int i = 0; i < point_count; ++i)
    {
      profile(i, 0) = sor_profile[i][0];
      profile(i, 1) = sor_profile[i][1];
    }
  }
  else if(shape == "cyl")
  {
    profile(0, 0) = -CYLINDER_HEIGHT / 2.0;
    profile(0, 1) = CYLINDER_RADIUS;
    profile(1, 0) = CYLINDER_HEIGHT / 2.0;
    profile(1, 1) = CYLINDER_RADIUS;
  }
  else if(shape == "cone")
  {
    profile(0, 0) = -CONE_HEIGHT / 2.0;
    profile(0, 1) = CONE_BASE_RADIUS;
    profile(1, 0) = CONE_HEIGHT / 2.0;
    profile(1, 1) = CONE_TOP_RADIUS;
  }
  else
  {
    SLIC_ERROR("Unsupported surface-of-revolution shape: " << shape);
  }
  return profile;
}

double sorVolume(axom::ArrayView<const double, 2> profile)
{
  double volume = 0.0;
  for(axom::IndexType i = 0; i < profile.shape()[0] - 1; ++i)
  {
    const primal::Cone<double, 3> section(profile(i, 1),
                                          profile(i + 1, 1),
                                          profile(i + 1, 0) - profile(i, 0));
    volume += section.volume();
  }
  return volume;
}

// Primitive volumes come from the same fixture definitions used to construct
// their geometries. The cupmesh value references the fixed data file below.
double expectedVolume(const std::string& shape)
{
  if(shape == "sphere")
  {
    return 4.0 * M_PI * SPHERE_RADIUS * SPHERE_RADIUS * SPHERE_RADIUS / 3.0;
  }
  if(shape == "cone" || shape == "cyl" || shape == "sor")
  {
    const auto profile = makeSorProfile(shape);
    return sorVolume(profile.view());
  }
  if(shape == "tet")
  {
    return makeTetFixture().volume();
  }
  if(shape == "hex")
  {
    return makeHexFixture().volume();
  }
  if(shape == "plane")
  {
    // A plane through the center of the symmetric box retains half its volume.
    constexpr double box_width = 2.0 * BOX_HALF_WIDTH;
    return box_width * box_width * box_width / 2.0;
  }
  if(shape == "tetmesh")
  {
    return 4.0 * TET_MESH_LENGTH * TET_MESH_LENGTH * TET_MESH_LENGTH;
  }
  if(shape == "cupmesh")
  {
    return 5.7521;
  }

  SLIC_ERROR("Unsupported benchmark shape: " << shape);
  return 0.0;
}

void skipWithError(benchmark::State& state, const std::string& message)
{
  s_benchmark_failed = true;
  state.SkipWithError(message.c_str());
}

klee::Geometry createSphereGeometry(int refinement_level)
{
  primal::Sphere<double, 3> sphere {Point3D {0.0, 0.0, 0.0}, SPHERE_RADIUS};
  return klee::Geometry(GEOMETRY_PROPERTIES, sphere, refinement_level, makeIdentityOperator());
}

klee::Geometry createSorGeometry(const std::string& shape, int refinement_level)
{
  // Keep the axis oblique to the Cartesian background mesh so the benchmark
  // does not exercise only the cheaper axis-aligned case. Each profile row is
  // an axial position followed by its radius.
  Point3D sor_base {0.0, 0.0, 0.0};
  Vector3D sor_direction {8.0, 4.0, 2.0};
  const bool is_general_sor = shape == "sor";
  auto profile = makeSorProfile(shape);
  klee::Geometry geometry(GEOMETRY_PROPERTIES,
                          profile,
                          sor_base,
                          sor_direction,
                          refinement_level,
                          makeIdentityOperator());
  if(is_general_sor)
  {
    // Three screening stages match the established workload and exercise the
    // hierarchy without making its depth another variable between runs.
    constexpr int screen_level = 3;
    geometry.asHierarchy()["screenLevel"] = screen_level;
  }
  return geometry;
}

klee::Geometry createTetGeometry()
{
  return klee::Geometry(GEOMETRY_PROPERTIES, makeTetFixture(), makeIdentityOperator());
}

klee::Geometry createHexGeometry()
{
  return klee::Geometry(GEOMETRY_PROPERTIES, makeHexFixture(), makeIdentityOperator());
}

klee::Geometry createPlaneGeometry()
{
  // The oblique plane avoids a grid-aligned special case and passes through the
  // box center so its expected clipped volume is exactly half the domain.
  const Vector3D normal = Vector3D {1.0, 2.0, 3.0}.unitVector();
  primal::Plane<double, 3> plane {normal, Point3D {0.0, 0.0, 0.0}, true};
  return klee::Geometry(GEOMETRY_PROPERTIES, plane, {nullptr});
}

klee::Geometry createBlueprintTetGeometry(sidre::Group* mesh_group)
{
  SLIC_ASSERT(axom::mint::blueprint::isValidRootGroup(mesh_group));
  mesh_group->destroyGroup("fields");
  klee::Geometry geometry(GEOMETRY_PROPERTIES, mesh_group, TOPO_NAME, makeIdentityOperator());
  geometry.asHierarchy()["fixOrientation"] = true;
  return geometry;
}

klee::Geometry createTetMeshGeometry(sidre::DataStore& datastore)
{
  // This small volume mesh exercises the mesh-based clipping strategy without
  // requiring an external data file.
  sidre::Group* mesh_group = datastore.getRoot()->createGroup("tetmesh_geometry");
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tet_mesh(3,
                                                                  axom::mint::CellType::TET,
                                                                  mesh_group,
                                                                  TOPO_NAME,
                                                                  COORDSET_NAME);

  tet_mesh.appendNode(-TET_MESH_LENGTH, -TET_MESH_LENGTH, -TET_MESH_LENGTH);
  tet_mesh.appendNode(TET_MESH_LENGTH, -TET_MESH_LENGTH, -TET_MESH_LENGTH);
  tet_mesh.appendNode(-TET_MESH_LENGTH, TET_MESH_LENGTH, -TET_MESH_LENGTH);
  tet_mesh.appendNode(-TET_MESH_LENGTH, -TET_MESH_LENGTH, TET_MESH_LENGTH);
  tet_mesh.appendNode(TET_MESH_LENGTH, TET_MESH_LENGTH, TET_MESH_LENGTH);
  tet_mesh.appendNode(-TET_MESH_LENGTH, TET_MESH_LENGTH, TET_MESH_LENGTH);
  tet_mesh.appendNode(TET_MESH_LENGTH, TET_MESH_LENGTH, -TET_MESH_LENGTH);
  tet_mesh.appendNode(TET_MESH_LENGTH, -TET_MESH_LENGTH, TET_MESH_LENGTH);
  axom::IndexType conn_0[4] = {0, 1, 2, 3};
  // One inverted cell verifies the mesh strategy's orientation repair path.
  axom::IndexType conn_1[4] = {4, 5, 7, 6};
  axom::IndexType conn_2[4] = {1, 2, 3, 5};
  tet_mesh.appendCell(conn_0);
  tet_mesh.appendCell(conn_1);
  tet_mesh.appendCell(conn_2);

  return createBlueprintTetGeometry(tet_mesh.getSidreGroup());
}

#ifdef AXOM_DATA_DIR
void fitTetMeshToDomain(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& tet_mesh)
{
  // Center and uniformly scale the input mesh so its workload is comparable
  // to the analytic shapes inside the [-2, 2]^3 background domain.
  double* coords[] = {tet_mesh.getCoordinateArray(0),
                      tet_mesh.getCoordinateArray(1),
                      tet_mesh.getCoordinateArray(2)};
  primal::BoundingBox<double, 3> bounds;
  for(axom::IndexType i = 0; i < tet_mesh.getNumberOfNodes(); ++i)
  {
    bounds.addPoint(Point3D {coords[0][i], coords[1][i], coords[2][i]});
  }

  const Point3D center = bounds.getCentroid();
  // Dividing the box width by sqrt(3) leaves enough clearance for arbitrary
  // rotation, matching the functional fixture's conservative fit.
  const double scale = (2.0 * BOX_HALF_WIDTH) / std::sqrt(3.0) / bounds.range().array().max();
  for(axom::IndexType i = 0; i < tet_mesh.getNumberOfNodes(); ++i)
  {
    for(int dim = 0; dim < 3; ++dim)
    {
      coords[dim][i] = (coords[dim][i] - center[dim]) * scale;
    }
  }
}

klee::Geometry createCupMeshGeometry(sidre::DataStore& datastore)
{
  // cup.proe is available only when Axom has a configured data directory.
  sidre::Group* mesh_group = datastore.getRoot()->createGroup("cupmesh_geometry");
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tet_mesh(3,
                                                                  axom::mint::CellType::TET,
                                                                  mesh_group,
                                                                  TOPO_NAME,
                                                                  COORDSET_NAME);
  quest::ProEReader reader;
  const std::string path = axom::utilities::filesystem::joinPath(AXOM_DATA_DIR, "quest/cup.proe");
  reader.setFileName(path);
  const int read_status = reader.read();
  SLIC_ERROR_IF(read_status != 0, "Could not read " << path);
  reader.getMesh(&tet_mesh);
  fitTetMeshToDomain(tet_mesh);

  return createBlueprintTetGeometry(tet_mesh.getSidreGroup());
}
#endif

klee::Geometry createGeometry(const std::string& shape,
                              sidre::DataStore& datastore,
                              int refinement_level)
{
  if(shape == "sphere")
  {
    return createSphereGeometry(refinement_level);
  }
  if(shape == "cone" || shape == "cyl" || shape == "sor")
  {
    return createSorGeometry(shape, refinement_level);
  }
  if(shape == "tet")
  {
    return createTetGeometry();
  }
  if(shape == "hex")
  {
    return createHexGeometry();
  }
  if(shape == "plane")
  {
    return createPlaneGeometry();
  }
  if(shape == "tetmesh")
  {
    return createTetMeshGeometry(datastore);
  }
#ifdef AXOM_DATA_DIR
  if(shape == "cupmesh")
  {
    return createCupMeshGeometry(datastore);
  }
#endif
  SLIC_ERROR("Unsupported benchmark shape: " << shape);
  return createSphereGeometry(refinement_level);
}

struct ClipProblem
{
  // The datastore owns the Blueprint groups referenced by ShapeMesh and by
  // mesh-based geometry strategies, so it must outlive every timed clip.
  sidre::DataStore datastore;
  sidre::Group* mesh_group {nullptr};
  std::shared_ptr<quest::experimental::ShapeMesh> shape_mesh;
  std::shared_ptr<quest::experimental::MeshClipperStrategy> strategy;
};

std::unique_ptr<ClipProblem> makeProblem(const std::string& shape,
                                         RuntimePolicy policy,
                                         int resolution,
                                         int refinement_level)
{
  // All work in this routine is benchmark setup and is outside the timed loop.
  auto problem = std::make_unique<ClipProblem>();

  const int alloc_id = axom::policyToDefaultAllocatorID(policy);
  problem->mesh_group = problem->datastore.getRoot()->createGroup("mesh");
  problem->mesh_group->setDefaultAllocator(alloc_id);

  primal::BoundingBox<double, 3> bbox(Point3D {-BOX_HALF_WIDTH, -BOX_HALF_WIDTH, -BOX_HALF_WIDTH},
                                      Point3D {BOX_HALF_WIDTH, BOX_HALF_WIDTH, BOX_HALF_WIDTH});
  axom::quest::util::make_unstructured_blueprint_box_mesh_3d(problem->mesh_group,
                                                             bbox,
                                                             {resolution, resolution, resolution},
                                                             TOPO_NAME,
                                                             COORDSET_NAME,
                                                             policy);
  problem->mesh_group->createGroup("state");

  problem->shape_mesh =
    std::make_shared<quest::experimental::ShapeMesh>(policy, alloc_id, problem->mesh_group, TOPO_NAME);
  problem->shape_mesh->precomputeMeshData();

  klee::Geometry geometry = createGeometry(shape, problem->datastore, refinement_level);
  problem->strategy = quest::experimental::util::make_clipper_strategy(geometry, shape);
  return problem;
}

double checksum(const axom::Array<double>& overlap)
{
  // Validation runs after timing. Copying to host makes the same checksum path
  // work for host and device allocations without affecting measured time.
  const int host_alloc_id = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  axom::Array<double> host_overlap(overlap, host_alloc_id);
  axom::ReduceSum<axom::SEQ_EXEC, double> sum(0.0);
  auto view = host_overlap.view();
  axom::for_all<axom::SEQ_EXEC>(host_overlap.size(),
                                [=] AXOM_HOST_DEVICE(axom::IndexType i) { sum += view[i]; });
  return sum.get();
}

void synchronize(RuntimePolicy policy)
{
  // Google Benchmark stops an iteration on the host. Explicit synchronization
  // ensures asynchronous device work is charged to that iteration.
  switch(policy)
  {
  case RuntimePolicy::seq:
    axom::synchronize<axom::SEQ_EXEC>();
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    axom::synchronize<axom::OMP_EXEC>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    axom::synchronize<axom::CUDA_EXEC<256>>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    axom::synchronize<axom::HIP_EXEC<256>>();
    break;
#endif
  default:
    SLIC_ERROR("Unsupported benchmark execution policy");
  }
}

void reportResult(benchmark::State& state,
                  const std::string& shape,
                  RuntimePolicy policy,
                  const ClipProblem& problem)
{
  // Use a fresh clipper for validation and counters so reported work is always
  // for one clip, independent of the benchmark iteration count.
  quest::experimental::MeshClipper validation_clipper(*problem.shape_mesh, problem.strategy);
  axom::Array<double> validation_overlap;
  validation_clipper.clip(validation_overlap);
  synchronize(policy);

  const double volume = checksum(validation_overlap);
  const double expected = expectedVolume(shape);
  const double relative_error = std::abs(volume - expected) / expected;
  if(!std::isfinite(volume) || relative_error > MAX_RELATIVE_ERROR)
  {
    const std::string message =
      axom::fmt::format("volume {} differs from reference {} by {:.3f}% (limit is {:.3f}%)",
                        volume,
                        expected,
                        100.0 * relative_error,
                        100.0 * MAX_RELATIVE_ERROR);
    skipWithError(state, message);
    return;
  }

  state.counters["volume"] = volume;
  state.counters["volume_error_pct"] = 100.0 * relative_error;
  state.counters["cells"] = problem.shape_mesh->getCellCount();
  const conduit::Node& stats = validation_clipper.getClippingStats();
  constexpr const char* clip_counters[] =
    {"clipsCandidates", "clipsIn", "clipsOn", "clipsOut", "clipsMiss", "clipsSum"};
  for(const char* counter : clip_counters)
  {
    if(stats.has_path(counter))
    {
      state.counters[counter] = static_cast<double>(stats[counter].to_int64());
    }
  }
  if((shape == "cone" || shape == "cyl" || shape == "sor") && stats.has_path("clipsCandidates") &&
     stats.has_path("clipsSum") && stats.has_path("clipsOn"))
  {
    const double root_count = static_cast<double>(stats["clipsCandidates"].to_int64());
    const double leaf_count = static_cast<double>(stats["clipsSum"].to_int64());
    const double unresolved_count = static_cast<double>(stats["clipsOn"].to_int64());
    if(leaf_count > root_count)
    {
      // Refining one octree leaf replaces it with eight children, increasing
      // the leaf count by seven.
      state.counters["adaptiveSubdivisions"] = (leaf_count - root_count) / 7.0;
      state.counters["unresolvedLeafFraction"] = unresolved_count / leaf_count;
    }
  }
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                          problem.shape_mesh->getCellCount());
}

void clipMesh(benchmark::State& state,
              const std::string& shape,
              RuntimePolicy policy,
              int resolution,
              int refinement_level)
{
  if(axom::policyToDefaultAllocatorID(policy) == axom::INVALID_ALLOCATOR_ID)
  {
    skipWithError(state, "Execution-space allocator not available");
    return;
  }

  auto problem = makeProblem(shape, policy, resolution, refinement_level);
  for(auto _ : state)
  {
    // The background mesh and clipping strategy are shared, while each timed
    // iteration includes clipper construction and first-use output allocation.
    quest::experimental::MeshClipper clipper(*problem->shape_mesh, problem->strategy);
    axom::Array<double> overlap;
    clipper.clip(overlap);
    synchronize(policy);
    benchmark::DoNotOptimize(overlap.data());
  }

  reportResult(state, shape, policy, *problem);
}

RuntimePolicy runtimePolicy(const std::string& policy_name)
{
  if(policy_name == "seq")
  {
    return RuntimePolicy::seq;
  }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  if(policy_name == "omp")
  {
    return RuntimePolicy::omp;
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(policy_name == "cuda")
  {
    return RuntimePolicy::cuda;
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(policy_name == "hip")
  {
    return RuntimePolicy::hip;
  }
#endif

  SLIC_ERROR("Unsupported benchmark execution policy: " << policy_name);
  return RuntimePolicy::seq;
}

void registerBenchmarks(const std::vector<std::string>& shapes,
                        const std::vector<int>& resolutions,
                        const std::vector<std::string>& policies,
                        int refinement_level)
{
  for(const auto& shape : shapes)
  {
    for(int resolution : resolutions)
    {
      for(const auto& policy_name : policies)
      {
        const RuntimePolicy policy = runtimePolicy(policy_name);
        const std::string name =
          axom::fmt::format("clipMesh/{}_{}_{}", shape, resolution, policy_name);
        ::benchmark::RegisterBenchmark(name.c_str(), [=](benchmark::State& state) {
          clipMesh(state, shape, policy, resolution, refinement_level);
        });
      }
    }
  }
}

}  // namespace

int main(int argc, char** argv)
{
  slic::SimpleLogger logger;
  slic::setLoggingMsgLevel(slic::message::Warning);

  // These defaults cover every geometry strategy available in the build.
  std::vector<std::string> shapes {"sphere", "cone", "cyl", "sor", "tet", "hex", "plane", "tetmesh"};
#ifdef AXOM_DATA_DIR
  shapes.push_back("cupmesh");
#endif
  const std::set<std::string> available_shapes(shapes.begin(), shapes.end());

  std::vector<int> resolutions {DEFAULT_RESOLUTION};
  int refinement_level = DEFAULT_REFINEMENT_LEVEL;
  std::vector<std::string> policies {"seq"};
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  policies.push_back("omp");
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  policies.push_back("cuda");
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  policies.push_back("hip");
#endif
  const std::set<std::string> available_policies(policies.begin(), policies.end());

  axom::CLI::App app {"Axom MeshClipper benchmarks"};
  app.add_option("-s,--shapes", shapes)
    ->description("Shapes to benchmark")
    ->expected(-1)
    ->check(axom::CLI::IsMember(available_shapes))
    ->capture_default_str();
  app.add_option("-r,--resolutions", resolutions)
    ->description("Cubic mesh resolutions to benchmark (at least 8)")
    ->expected(-1)
    ->each([](const std::string& value) {
      // Smaller meshes leave too few cells across the compact shapes and make
      // setup costs dominate the clipping work.
      constexpr int min_resolution = 8;
      if(std::stoi(value) < min_resolution)
      {
        throw axom::CLI::ValidationError("Mesh resolutions below 8 are not supported");
      }
    })
    ->capture_default_str();
  app.add_option("--refinements", refinement_level)
    ->description("Refinement level for curved geometries (at least 5)")
    ->each([](const std::string& value) {
      // Lower levels exceed the volume-validation tolerance for curved shapes.
      constexpr int min_refinement_level = 5;
      if(std::stoi(value) < min_refinement_level)
      {
        throw axom::CLI::ValidationError("Refinement levels below 5 are not supported");
      }
    })
    ->capture_default_str();
  app.add_option("-p,--policies", policies)
    ->description("Execution policies to benchmark")
    ->expected(-1)
    ->check(axom::CLI::IsMember(available_policies))
    ->capture_default_str();
  app.allow_extras();
  CLI11_PARSE(app, argc, argv);

  // CLI11 consumes benchmark-specific options. Forward all remaining options
  // so standard Google Benchmark controls such as filters and repetitions work.
  std::vector<std::string> benchmark_args = app.remaining_for_passthrough();
  benchmark_args.insert(benchmark_args.begin(), argv[0]);
  std::vector<char*> benchmark_argv;
  benchmark_argv.reserve(benchmark_args.size());
  for(auto& arg : benchmark_args)
  {
    benchmark_argv.push_back(arg.data());
  }
  int benchmark_argc = static_cast<int>(benchmark_argv.size());
  ::benchmark::Initialize(&benchmark_argc, benchmark_argv.data());
  if(::benchmark::ReportUnrecognizedArguments(benchmark_argc, benchmark_argv.data()))
  {
    return 1;
  }

  // Avoid registering duplicate names when a value is repeated on the command
  // line, and keep output ordering deterministic.
  std::sort(shapes.begin(), shapes.end());
  shapes.erase(std::unique(shapes.begin(), shapes.end()), shapes.end());
  std::sort(resolutions.begin(), resolutions.end());
  resolutions.erase(std::unique(resolutions.begin(), resolutions.end()), resolutions.end());
  std::sort(policies.begin(), policies.end());
  policies.erase(std::unique(policies.begin(), policies.end()), policies.end());

  registerBenchmarks(shapes, resolutions, policies, refinement_level);
  ::benchmark::RunSpecifiedBenchmarks();

  return s_benchmark_failed ? 1 : 0;
}
