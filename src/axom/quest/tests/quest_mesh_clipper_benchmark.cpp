// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/CLI11.hpp"
#include "axom/core.hpp"
#include "axom/klee.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"
#include "axom/quest/MeshClipper.hpp"
#include "axom/quest/ShapeMesh.hpp"
#include "axom/quest/util/mesh_helpers.hpp"
#include "axom/quest/util/make_clipper_strategy.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"

#include "benchmark/benchmark.h"

#include <algorithm>
#include <cmath>
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

constexpr const char* topoName = "mesh";
constexpr const char* coordsetName = "coords";

klee::Geometry createSphereGeometry(int refinementLevel)
{
  klee::TransformableGeometryProperties prop {klee::Dimensions::Three, klee::LengthUnit::unspecified};
  primal::Sphere<double, 3> sphere {Point3D {0.0, 0.0, 0.0}, 1.0};
  auto compositeOp = std::make_shared<klee::CompositeOperator>(prop);
  return klee::Geometry(prop, sphere, refinementLevel, compositeOp);
}

klee::Geometry createSorGeometry(const std::string& shape, int refinementLevel)
{
  klee::TransformableGeometryProperties prop {klee::Dimensions::Three, klee::LengthUnit::unspecified};
  Point3D sorBase {0.0, 0.0, 0.0};
  Vector3D sorDirection {8.0, 4.0, 2.0};
  const bool isGeneralSor = shape == "sor";
  const int pointCount = isGeneralSor ? 12 : 2;
  axom::Array<double, 2> discreteFunction({pointCount, 2}, axom::ArrayStrideOrder::ROW);

  if(isGeneralSor)
  {
    constexpr double profile[][2] = {{-1.2, 1.1},
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
    for(int i = 0; i < pointCount; ++i)
    {
      discreteFunction(i, 0) = profile[i][0];
      discreteFunction(i, 1) = profile[i][1];
    }
  }
  else if(shape == "cyl")
  {
    constexpr double radius = 0.695;
    constexpr double height = 2.78;
    discreteFunction(0, 0) = -height / 2.0;
    discreteFunction(0, 1) = radius;
    discreteFunction(1, 0) = height / 2.0;
    discreteFunction(1, 1) = radius;
  }
  else
  {
    constexpr double baseRadius = 1.23;
    constexpr double topRadius = 0.176;
    constexpr double height = 2.3;
    discreteFunction(0, 0) = -height / 2.0;
    discreteFunction(0, 1) = baseRadius;
    discreteFunction(1, 0) = height / 2.0;
    discreteFunction(1, 1) = topRadius;
  }

  auto compositeOp = std::make_shared<klee::CompositeOperator>(prop);
  klee::Geometry geometry(prop, discreteFunction, sorBase, sorDirection, refinementLevel, compositeOp);
  if(isGeneralSor)
  {
    geometry.asHierarchy()["screenLevel"] = 3;
  }
  return geometry;
}

klee::Geometry createTetGeometry()
{
  klee::TransformableGeometryProperties prop {klee::Dimensions::Three, klee::LengthUnit::unspecified};
  constexpr double length = 1.55;
  const Point3D a {Point3D::NumericArray {0.8, 0.0, -1.0} * length};
  const Point3D b {Point3D::NumericArray {-0.8, 1.0, -1.0} * length};
  const Point3D c {Point3D::NumericArray {-0.8, -1.0, -1.0} * length};
  const Point3D d {Point3D::NumericArray {0.0, 0.0, 1.0} * length};
  primal::Tetrahedron<double, 3> tet {a, b, c, d};
  auto compositeOp = std::make_shared<klee::CompositeOperator>(prop);
  return klee::Geometry(prop, tet, compositeOp);
}

klee::Geometry createHexGeometry()
{
  klee::TransformableGeometryProperties prop {klee::Dimensions::Three, klee::LengthUnit::unspecified};
  constexpr double medium = 0.82;
  constexpr double large = 1.2 * medium;
  constexpr double small = 0.8 * medium;
  primal::Hexahedron<double, 3> hex {Point3D {-large, -medium, -small},
                                     Point3D {large, -medium, -small},
                                     Point3D {large, medium, -small},
                                     Point3D {-large, medium, -small},
                                     Point3D {-large, -medium, small},
                                     Point3D {large, -medium, small},
                                     Point3D {large, medium, small},
                                     Point3D {-large, medium, small}};
  auto compositeOp = std::make_shared<klee::CompositeOperator>(prop);
  return klee::Geometry(prop, hex, compositeOp);
}

klee::Geometry createPlaneGeometry()
{
  klee::TransformableGeometryProperties prop {klee::Dimensions::Three, klee::LengthUnit::unspecified};
  const Vector3D normal = Vector3D {1.0, 2.0, 3.0}.unitVector();
  primal::Plane<double, 3> plane {normal, Point3D {0.0, 0.0, 0.0}, true};
  return klee::Geometry(prop, plane, {nullptr});
}

klee::Geometry createTetMeshGeometry(sidre::DataStore& datastore)
{
  sidre::Group* meshGroup = datastore.getRoot()->createGroup("tetmesh_geometry");
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tetMesh(3,
                                                                 axom::mint::CellType::TET,
                                                                 meshGroup,
                                                                 topoName,
                                                                 coordsetName);

  constexpr double length = 1.17;
  tetMesh.appendNode(-length, -length, -length);
  tetMesh.appendNode(length, -length, -length);
  tetMesh.appendNode(-length, length, -length);
  tetMesh.appendNode(-length, -length, length);
  tetMesh.appendNode(length, length, length);
  tetMesh.appendNode(-length, length, length);
  tetMesh.appendNode(length, length, -length);
  tetMesh.appendNode(length, -length, length);
  axom::IndexType conn0[4] = {0, 1, 2, 3};
  axom::IndexType conn1[4] = {4, 5, 7, 6};
  axom::IndexType conn2[4] = {1, 2, 3, 5};
  tetMesh.appendCell(conn0);
  tetMesh.appendCell(conn1);
  tetMesh.appendCell(conn2);

  SLIC_ASSERT(axom::mint::blueprint::isValidRootGroup(meshGroup));
  meshGroup->destroyGroup("fields");
  klee::TransformableGeometryProperties prop {klee::Dimensions::Three, klee::LengthUnit::unspecified};
  auto compositeOp = std::make_shared<klee::CompositeOperator>(prop);
  klee::Geometry geometry(prop, tetMesh.getSidreGroup(), topoName, compositeOp);
  geometry.asHierarchy()["fixOrientation"] = true;
  return geometry;
}

#ifdef AXOM_DATA_DIR
void fitTetMeshToDomain(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& tetMesh)
{
  double* coords[] = {tetMesh.getCoordinateArray(0),
                      tetMesh.getCoordinateArray(1),
                      tetMesh.getCoordinateArray(2)};
  primal::BoundingBox<double, 3> bounds;
  for(axom::IndexType i = 0; i < tetMesh.getNumberOfNodes(); ++i)
  {
    bounds.addPoint(Point3D {coords[0][i], coords[1][i], coords[2][i]});
  }

  const Point3D center = bounds.getCentroid();
  const double scale = 4.0 / std::sqrt(3.0) / bounds.range().array().max();
  for(axom::IndexType i = 0; i < tetMesh.getNumberOfNodes(); ++i)
  {
    for(int dim = 0; dim < 3; ++dim)
    {
      coords[dim][i] = (coords[dim][i] - center[dim]) * scale;
    }
  }
}

klee::Geometry createCupMeshGeometry(sidre::DataStore& datastore)
{
  sidre::Group* meshGroup = datastore.getRoot()->createGroup("cupmesh_geometry");
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tetMesh(3,
                                                                 axom::mint::CellType::TET,
                                                                 meshGroup,
                                                                 topoName,
                                                                 coordsetName);
  quest::ProEReader reader;
  const std::string path = axom::utilities::filesystem::joinPath(AXOM_DATA_DIR, "quest/cup.proe");
  reader.setFileName(path);
  const int readStatus = reader.read();
  SLIC_ERROR_IF(readStatus != 0, "Could not read " << path);
  reader.getMesh(&tetMesh);
  fitTetMeshToDomain(tetMesh);

  SLIC_ASSERT(axom::mint::blueprint::isValidRootGroup(meshGroup));
  meshGroup->destroyGroup("fields");
  klee::TransformableGeometryProperties prop {klee::Dimensions::Three, klee::LengthUnit::unspecified};
  auto compositeOp = std::make_shared<klee::CompositeOperator>(prop);
  klee::Geometry geometry(prop, tetMesh.getSidreGroup(), topoName, compositeOp);
  geometry.asHierarchy()["fixOrientation"] = true;
  return geometry;
}
#endif

klee::Geometry createGeometry(const std::string& shape, sidre::DataStore& datastore, int refinementLevel)
{
  if(shape == "sphere")
  {
    return createSphereGeometry(refinementLevel);
  }
  if(shape == "cone" || shape == "cyl" || shape == "sor")
  {
    return createSorGeometry(shape, refinementLevel);
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
  return createSphereGeometry(refinementLevel);
}

struct ClipProblem
{
  sidre::DataStore datastore;
  sidre::Group* meshGroup {nullptr};
  std::shared_ptr<quest::experimental::ShapeMesh> shapeMesh;
  std::shared_ptr<quest::experimental::MeshClipperStrategy> strategy;
  axom::Array<double> overlap;
};

std::unique_ptr<ClipProblem> makeProblem(const std::string& shape,
                                         RuntimePolicy policy,
                                         int resolution,
                                         int refinementLevel)
{
  auto problem = std::make_unique<ClipProblem>();

  const int allocId = axom::policyToDefaultAllocatorID(policy);
  problem->meshGroup = problem->datastore.getRoot()->createGroup("mesh");
  problem->meshGroup->setDefaultAllocator(allocId);

  primal::BoundingBox<double, 3> bbox(Point3D {-2.0, -2.0, -2.0}, Point3D {2.0, 2.0, 2.0});
  axom::quest::util::make_unstructured_blueprint_box_mesh_3d(problem->meshGroup,
                                                             bbox,
                                                             {resolution, resolution, resolution},
                                                             topoName,
                                                             coordsetName,
                                                             policy);
  problem->meshGroup->createGroup("state");

  problem->shapeMesh =
    std::make_shared<quest::experimental::ShapeMesh>(policy, allocId, problem->meshGroup, topoName);
  problem->shapeMesh->precomputeMeshData();

  klee::Geometry geometry = createGeometry(shape, problem->datastore, refinementLevel);
  problem->strategy = quest::experimental::util::make_clipper_strategy(geometry, shape);
  return problem;
}

double checksum(const axom::Array<double>& overlap)
{
  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  axom::Array<double> hostOverlap(overlap, hostAllocId);
  axom::ReduceSum<axom::SEQ_EXEC, double> sum(0.0);
  auto view = hostOverlap.view();
  axom::for_all<axom::SEQ_EXEC>(hostOverlap.size(),
                                [=] AXOM_HOST_DEVICE(axom::IndexType i) { sum += view[i]; });
  return sum.get();
}

void synchronize(RuntimePolicy policy)
{
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

void clipMesh(benchmark::State& state,
              const std::string& shape,
              RuntimePolicy policy,
              int resolution,
              int refinementLevel)
{
  if(axom::policyToDefaultAllocatorID(policy) == axom::INVALID_ALLOCATOR_ID)
  {
    state.SkipWithError("Execution-space allocator not available");
    return;
  }

  auto problem = makeProblem(shape, policy, resolution, refinementLevel);
  quest::experimental::MeshClipper clipper(*problem->shapeMesh, problem->strategy);

  for(auto _ : state)
  {
    clipper.clip(problem->overlap);
    synchronize(policy);
    benchmark::DoNotOptimize(problem->overlap.data());
  }

  const double volume = checksum(problem->overlap);
  if(!std::isfinite(volume) || volume <= 0.0)
  {
    state.SkipWithError("Clipping produced an invalid volume");
    return;
  }

  state.counters["volume"] = volume;
  state.counters["cells"] = problem->shapeMesh->getCellCount();
  const conduit::Node& stats = clipper.getClippingStats();
  constexpr const char* clipCounters[] =
    {"clipsCandidates", "clipsIn", "clipsOn", "clipsOut", "clipsMiss", "clipsSum"};
  const bool statsAccumulateAcrossIterations = shape == "cone" || shape == "cyl" || shape == "sor";
  for(const char* counter : clipCounters)
  {
    if(stats.has_path(counter))
    {
      const auto flags = statsAccumulateAcrossIterations ? benchmark::Counter::kAvgIterations
                                                         : benchmark::Counter::kDefaults;
      state.counters[counter] =
        benchmark::Counter(static_cast<double>(stats[counter].to_int64()), flags);
    }
  }
  if((shape == "cone" || shape == "cyl" || shape == "sor") && stats.has_path("clipsCandidates") &&
     stats.has_path("clipsSum") && stats.has_path("clipsOn"))
  {
    const double rootCount = static_cast<double>(stats["clipsCandidates"].to_int64());
    const double leafCount = static_cast<double>(stats["clipsSum"].to_int64());
    const double unresolvedCount = static_cast<double>(stats["clipsOn"].to_int64());
    if(leafCount > rootCount)
    {
      state.counters["adaptiveSubdivisions"] =
        benchmark::Counter((leafCount - rootCount) / 7.0, benchmark::Counter::kAvgIterations);
      state.counters["unresolvedLeafFraction"] = unresolvedCount / leafCount;
    }
  }
  state.SetItemsProcessed(state.iterations() * problem->shapeMesh->getCellCount());
}

}  // namespace

int main(int argc, char** argv)
{
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Warning);

  std::vector<std::string> shapes {"sphere", "cone", "cyl", "sor", "tet", "hex", "plane", "tetmesh"};
#ifdef AXOM_DATA_DIR
  shapes.push_back("cupmesh");
#endif
  const std::set<std::string> availableShapes(shapes.begin(), shapes.end());

  std::vector<int> resolutions {51};
  int refinementLevel = 5;
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
  const std::set<std::string> availablePolicies(policies.begin(), policies.end());

  axom::CLI::App app {"Axom MeshClipper benchmarks"};
  app.add_option("-s,--shapes", shapes)
    ->description("Shapes to benchmark")
    ->expected(-1)
    ->check(axom::CLI::IsMember(availableShapes))
    ->capture_default_str();
  app.add_option("-r,--resolutions", resolutions)
    ->description("Cubic mesh resolutions to benchmark")
    ->expected(-1)
    ->check(axom::CLI::PositiveNumber)
    ->capture_default_str();
  app.add_option("--refinements", refinementLevel)
    ->description("Refinement level for sphere and SOR geometries")
    ->check(axom::CLI::NonNegativeNumber)
    ->capture_default_str();
  app.add_option("-p,--policies", policies)
    ->description("Execution policies to benchmark")
    ->expected(-1)
    ->check(axom::CLI::IsMember(availablePolicies))
    ->capture_default_str();
  app.allow_extras();
  CLI11_PARSE(app, argc, argv);

  std::vector<std::string> benchmarkArgs = app.remaining_for_passthrough();
  benchmarkArgs.insert(benchmarkArgs.begin(), argv[0]);
  std::vector<char*> benchmarkArgv;
  benchmarkArgv.reserve(benchmarkArgs.size());
  for(auto& arg : benchmarkArgs)
  {
    benchmarkArgv.push_back(arg.data());
  }
  int benchmarkArgc = static_cast<int>(benchmarkArgv.size());
  ::benchmark::Initialize(&benchmarkArgc, benchmarkArgv.data());
  if(::benchmark::ReportUnrecognizedArguments(benchmarkArgc, benchmarkArgv.data()))
  {
    slic::finalize();
    return 1;
  }

  std::sort(shapes.begin(), shapes.end());
  shapes.erase(std::unique(shapes.begin(), shapes.end()), shapes.end());
  std::sort(resolutions.begin(), resolutions.end());
  resolutions.erase(std::unique(resolutions.begin(), resolutions.end()), resolutions.end());
  std::sort(policies.begin(), policies.end());
  policies.erase(std::unique(policies.begin(), policies.end()), policies.end());

  for(const auto& shape : shapes)
  {
    for(int resolution : resolutions)
    {
      for(const auto& policyName : policies)
      {
        RuntimePolicy policy = RuntimePolicy::seq;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
        if(policyName == "omp")
        {
          policy = RuntimePolicy::omp;
        }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
        if(policyName == "cuda")
        {
          policy = RuntimePolicy::cuda;
        }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
        if(policyName == "hip")
        {
          policy = RuntimePolicy::hip;
        }
#endif

        const std::string name =
          axom::fmt::format("clipMesh/{}_ares_{}_{}", shape, resolution, policyName);
        ::benchmark::RegisterBenchmark(name.c_str(), [=](benchmark::State& state) {
          clipMesh(state, shape, policy, resolution, refinementLevel);
        });
      }
    }
  }

  ::benchmark::RunSpecifiedBenchmarks();

  slic::finalize();
  return 0;
}
