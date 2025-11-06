// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: device_spatial_indexes.cpp
///
/// This example uses a spatial index, the linear BVH tree from Axom's spin
/// component, in addition to RAJA and Umpire based kernels for a highly
//  efficient performance-portable self-intersection algorithm.
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "../patch/hip_patch.hpp"

#include "axom/core.hpp"
#include "axom/slic.hpp"

#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#else
  #error This example requires axom to be configured with RAJA support
#endif

#ifdef AXOM_USE_UMPIRE
  #include "umpire/Umpire.hpp"
#else
  #error This example requires axom to be configured with Umpire support
#endif

#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <memory>
#include <array>

using RuntimePolicy = axom::runtime_policy::Policy;

//-----------------------------------------------------------------------------
/// Basic RAII utility class for initializing and finalizing slic logger
//-----------------------------------------------------------------------------
struct BasicLogger
{
  BasicLogger()
  {
    namespace slic = axom::slic;

    // Initialize the SLIC logger
    slic::initialize();
    slic::setLoggingMsgLevel(slic::message::Debug);

    // Customize logging levels and formatting
    const std::string slicFormatStr = "[lesson_04: <LEVEL>] <MESSAGE> \n";

    slic::addStreamToMsgLevel(new slic::GenericOutputStream(&std::cerr), slic::message::Error);
    slic::addStreamToMsgLevel(new slic::GenericOutputStream(&std::cerr, slicFormatStr),
                              slic::message::Warning);

    auto* compactStream = new slic::GenericOutputStream(&std::cout, slicFormatStr);
    slic::addStreamToMsgLevel(compactStream, slic::message::Info);
    slic::addStreamToMsgLevel(compactStream, slic::message::Debug);
  }

  ~BasicLogger() { axom::slic::finalize(); }
};

//-----------------------------------------------------------------------------
/// Struct to help with parsing and storing command line args
//-----------------------------------------------------------------------------

struct Input
{
  static const std::map<std::string, RuntimePolicy> s_validPolicies;

  std::string mesh_file {""};
  bool verboseOutput {false};
  double weldThreshold {1e-6};
  double intersectionThreshold {1e-08};
  RuntimePolicy policy {RuntimePolicy::seq};

  void parse(int argc, char** argv, axom::CLI::App& app);
  bool isVerbose() const { return verboseOutput; }
};

void Input::parse(int argc, char** argv, axom::CLI::App& app)
{
  app.add_option("-i, --infile", mesh_file)
    ->description("The input STL mesh file")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_flag("-v,--verbose", verboseOutput)
    ->description("Increase logging verbosity?")
    ->capture_default_str();

  app.add_option("--weld-threshold", weldThreshold)
    ->description(
      "Threshold to use when welding vertices.\n"
      "Will skip if not strictly positive.")
    ->capture_default_str();

  app.add_option("--intersection-threshold", intersectionThreshold)
    ->description("Threshold to use when testing for intersecting triangles")
    ->capture_default_str();

  std::stringstream pol_sstr;
  pol_sstr << "Set runtime policy for intersection-based sampling method.";
  pol_sstr << "\nSet to 'seq' or 0 to use the RAJA sequential policy.";
#ifdef AXOM_USE_OPENMP
  pol_sstr << "\nSet to 'omp' or 1 to use the RAJA OpenMP policy.";
#endif
#ifdef AXOM_USE_CUDA
  pol_sstr << "\nSet to 'cuda' or 2 to use the RAJA CUDA policy.";
#endif
#ifdef AXOM_USE_HIP
  pol_sstr << "\nSet to 'hip' or 3 to use the RAJA HIP policy.";
#endif

  app.add_option("-p, --policy", policy)
    ->description(pol_sstr.str())
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

  app.get_formatter()->column_width(40);

  app.parse(argc, argv);  // Could throw an exception

  // Output parsed information
  SLIC_INFO(axom::fmt::format(
    R"(
     Parsed parameters:
      * STL mesh: '{}'
      * Threshold for welding: {}
      * Skip welding: {}
      * Threshold for intersections: {}
      * Verbose logging: {}
      * Runtime execution policy: '{}'
      )",
    mesh_file,
    weldThreshold,
    (weldThreshold <= 0.),
    intersectionThreshold,
    verboseOutput,
    axom::runtime_policy::policyToName(policy)));
}

//-----------------------------------------------------------------------------
/// Basic triangle mesh to be used in our application
//-----------------------------------------------------------------------------
struct TriangleMesh
{
  using Point = axom::primal::Point<double, 3>;
  using Triangle = axom::primal::Triangle<double, 3>;
  using BoundingBox = axom::primal::BoundingBox<double, 3>;

  axom::IndexType numTriangles() const { return m_triangles.size(); }
  axom::Array<Triangle>& triangles() { return m_triangles; }
  const axom::Array<Triangle>& triangles() const { return m_triangles; }

  BoundingBox& meshBoundingBox() { return m_meshBoundingBox; }
  const BoundingBox& meshBoundingBox() const { return m_meshBoundingBox; }

  axom::Array<BoundingBox>& triangleBoundingBoxes() { return m_triangleBoundingBoxes; }
  const axom::Array<BoundingBox>& triangleBoundingBoxes() const { return m_triangleBoundingBoxes; }

  axom::Array<Triangle> m_triangles;
  axom::Array<BoundingBox> m_triangleBoundingBoxes;
  BoundingBox m_meshBoundingBox;
};

TriangleMesh makeTriangleMesh(const std::string& stl_mesh_path, double weldThreshold)
{
  TriangleMesh triMesh;

  // load STL mesh into a mint unstructured mesh
  auto* surface_mesh =
    new axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>(3, axom::mint::TRIANGLE);
  {
    axom::utilities::Timer timer(true);

    auto reader = std::make_unique<axom::quest::STLReader>();
    reader->setFileName(stl_mesh_path);
    reader->read();
    reader->getMesh(surface_mesh);

    timer.stop();
    SLIC_INFO(axom::fmt::format("Loading the mesh took {:4.3} seconds.", timer.elapsedTimeInSec()));
  }

  // optionally weld triangle mesh
  if(weldThreshold > 0.)
  {
    axom::utilities::Timer timer(true);
    axom::quest::weldTriMeshVertices(&surface_mesh, weldThreshold);
    timer.stop();

    SLIC_INFO(axom::fmt::format("Vertex welding took {:4.3} seconds.", timer.elapsedTimeInSec()));
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "After welding, mesh has {:L} vertices and {:L} triangles.",
                                surface_mesh->getNumberOfNodes(),
                                surface_mesh->getNumberOfCells()));
  }

  // extract triangles into an axom::Array
  const int numCells = surface_mesh->getNumberOfCells();
  triMesh.m_triangles.reserve(numCells);
  {
    TriangleMesh::Triangle tri;
    std::array<axom::IndexType, 3> triCell;
    for(int i = 0; i < numCells; ++i)
    {
      surface_mesh->getCellNodeIDs(i, triCell.data());
      surface_mesh->getNode(triCell[0], tri[0].data());
      surface_mesh->getNode(triCell[1], tri[1].data());
      surface_mesh->getNode(triCell[2], tri[2].data());

      triMesh.m_triangles.emplace_back(tri);
    }
  }

  delete surface_mesh;
  surface_mesh = nullptr;

  // compute and store triangle bounding boxes and mesh bounding box
  triMesh.m_triangleBoundingBoxes.reserve(numCells);
  for(const auto& tri : triMesh.triangles())
  {
    triMesh.m_triangleBoundingBoxes.emplace_back(axom::primal::compute_bounding_box(tri));
    triMesh.m_meshBoundingBox.addBox(triMesh.m_triangleBoundingBoxes.back());
  }

  SLIC_INFO(axom::fmt::format("Mesh bounding box is {}.", triMesh.meshBoundingBox()));

  return triMesh;
}

using IndexPair = std::pair<axom::IndexType, axom::IndexType>;

template <typename ExecSpace>
axom::Array<IndexPair> findIntersectionsBVH(const TriangleMesh& triMesh,
                                            double tol,
                                            bool verboseOutput = false)
{
  SLIC_INFO("Running naive intersection algorithm in execution Space: "
            << axom::execution_space<ExecSpace>::name());

  using TriangleArray = axom::Array<typename TriangleMesh::Triangle>;
  using BBoxArray = axom::Array<typename TriangleMesh::BoundingBox>;
  using IndexArray = axom::Array<axom::IndexType>;
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

  axom::Array<IndexPair> intersectionPairs;

  // Get ids of necessary allocators
  const int host_allocator = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  const int kernel_allocator = on_device
    ? axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    : axom::execution_space<ExecSpace>::allocatorID();

  // Copy the triangles to the device, if necessary
  // Either way, tris_v will be a view w/ data in the correct space
  auto& tris_h = triMesh.triangles();
  TriangleArray tris_d = on_device ? TriangleArray(tris_h, kernel_allocator) : TriangleArray();
  auto tris_v = on_device ? tris_d.view() : tris_h.view();

  // Copy the bboxes to the device, if necessary
  // Either way, bbox_v will be a view w/ data in the correct space
  auto& bbox_h = triMesh.triangleBoundingBoxes();
  BBoxArray bbox_d = on_device ? BBoxArray(bbox_h, kernel_allocator) : BBoxArray();
  auto bbox_v = on_device ? bbox_d.view() : bbox_h.view();

  axom::utilities::Timer timer;

  // Initialize a BVH tree over the triangle bounding boxes
  timer.start();
  axom::spin::BVH<3, ExecSpace, double> bvh;
  bvh.setAllocatorID(kernel_allocator);
  bvh.initialize(bbox_v, bbox_v.size());
  timer.stop();
  SLIC_INFO_IF(
    verboseOutput,
    axom::fmt::format("0: Initializing BVH took {:4.3} seconds.", timer.elapsedTimeInSec()));

  // Search for intersecting bounding boxes of triangles;
  // result is returned as CSR arrays for candidate data
  timer.start();
  IndexArray offsets_d(bbox_v.size(), bbox_v.size(), kernel_allocator);
  IndexArray counts_d(bbox_v.size(), bbox_v.size(), kernel_allocator);
  IndexArray candidates_d(0, 0, kernel_allocator);

  auto offsets_v = offsets_d.view();
  auto counts_v = counts_d.view();
  bvh.findBoundingBoxes(offsets_v, counts_v, candidates_d, bbox_v.size(), bbox_v);

  timer.stop();
  SLIC_INFO_IF(verboseOutput,
               axom::fmt::format("1: Querying candidate bounding boxes took {:4.3} seconds.",
                                 timer.elapsedTimeInSec()));

  // Expand candidate list into corresponding arrays of indices
  // Only keep results where candidate has a greater index than triangle
  // and when both are non-degenerate
  IndexArray indices_d(axom::ArrayOptions::Uninitialized {},
                       candidates_d.size(),
                       candidates_d.size(),
                       kernel_allocator);

  IndexArray validCandidates_d(axom::ArrayOptions::Uninitialized {},
                               candidates_d.size(),
                               candidates_d.size(),
                               kernel_allocator);

  axom::IndexType numCandidates {};
  timer.start();
  {
    const int totalTriangles = triMesh.numTriangles();

    IndexArray numValidCandidates_d(1, 1, kernel_allocator);
    auto* numValidCandidates_p = numValidCandidates_d.data();

    auto indices_v = indices_d.view();
    auto validCandidates_v = validCandidates_d.view();
    auto candidates_v = candidates_d.view();

    // Compute a device bool array of validity flags
    axom::Array<bool> is_valid_d(axom::ArrayOptions::Uninitialized {},
                                 bbox_v.size(),
                                 bbox_v.size(),
                                 kernel_allocator);
    auto is_valid_v = is_valid_d.view();

    axom::for_all<ExecSpace>(
      totalTriangles,
      AXOM_LAMBDA(axom::IndexType i) { is_valid_v[i] = !tris_v[i].degenerate(); });

    // Keep pairs of valid triangles whose bounding boxes overlap
    axom::for_all<ExecSpace>(
      totalTriangles,
      AXOM_LAMBDA(axom::IndexType i) {
        for(int j = 0; j < counts_v[i]; j++)
        {
          const axom::IndexType potential = candidates_v[offsets_v[i] + j];
          if(i < potential && is_valid_v[i] && is_valid_v[potential])
          {
            const auto idx =
              axom::atomicAdd<axom::auto_atomic>(numValidCandidates_p, axom::IndexType {1});
            indices_v[idx] = i;
            validCandidates_v[idx] = potential;
          }
        }
      });

    axom::copy(&numCandidates, numValidCandidates_p, sizeof(axom::IndexType));
  }
  timer.stop();
  SLIC_INFO_IF(verboseOutput,
               axom::fmt::format("2: Filtering invalid candidates took {:4.3} seconds.",
                                 timer.elapsedTimeInSec()));

  // Iterate through valid candidates to find actual intersections
  IndexArray intersect_d[2] = {
    IndexArray(axom::ArrayOptions::Uninitialized {}, numCandidates, numCandidates, kernel_allocator),
    IndexArray(axom::ArrayOptions::Uninitialized {}, numCandidates, numCandidates, kernel_allocator)};
  axom::IndexType numIntersections {};
  timer.start();
  {
    auto intersect1_v = intersect_d[0].view();
    auto intersect2_v = intersect_d[1].view();

    IndexArray numIntersections_d(1, 1, kernel_allocator);
    auto* numIntersections_p = numIntersections_d.data();

    auto indices_v = indices_d.view();
    auto validCandidates_v = validCandidates_d.view();

    // Perform triangle-triangle tests
    axom::for_all<ExecSpace>(
      numCandidates,
      AXOM_LAMBDA(axom::IndexType i) {
        constexpr bool includeBoundaries = false;
        const auto index = indices_v[i];
        const auto candidate = validCandidates_v[i];
        if(axom::primal::intersect(tris_v[index], tris_v[candidate], includeBoundaries, tol))
        {
          const auto idx =
            axom::atomicAdd<axom::auto_atomic>(numIntersections_p, axom::IndexType {1});
          intersect1_v[idx] = index;
          intersect2_v[idx] = candidate;
        }
      });

    axom::copy(&numIntersections, numIntersections_p, sizeof(axom::IndexType));
  }
  intersect_d[0].resize(numIntersections);
  intersect_d[1].resize(numIntersections);

  timer.stop();
  SLIC_INFO_IF(verboseOutput,
               axom::fmt::format("3: Finding actual intersections took {:4.3} seconds.",
                                 timer.elapsedTimeInSec()));

  SLIC_INFO_IF(verboseOutput,
               axom::fmt::format(
                 axom::utilities::locale(),
                 R"(Stats for self-intersection query
    -- Number of mesh triangles {:L}
    -- Total possible candidates {:L}
    -- Candidates from BVH query {:L}
    -- Potential candidates after filtering {:L}
    -- Actual intersections {:L}
    )",
                 triMesh.numTriangles(),
                 [](std::int64_t sz) { return sz * (sz - 1) / 2.; }(triMesh.numTriangles()),
                 candidates_d.size(),
                 numCandidates,
                 numIntersections));

  // copy results back to host and into return vector
  IndexArray intersect_h[2] = {on_device ? IndexArray(intersect_d[0], host_allocator) : IndexArray(),
                               on_device ? IndexArray(intersect_d[1], host_allocator) : IndexArray()};

  auto intersect1_h_v = on_device ? intersect_h[0].view() : intersect_d[0].view();
  auto intersect2_h_v = on_device ? intersect_h[1].view() : intersect_d[1].view();

  for(axom::IndexType idx = 0; idx < numIntersections; ++idx)
  {
    intersectionPairs.emplace_back(std::make_pair(intersect1_h_v[idx], intersect2_h_v[idx]));
  }

  return intersectionPairs;
}

int main(int argc, char** argv)
{
  // Initialize logger; use RAII so it will finalize at the end of the application
  BasicLogger logger;

  // Parse the command line arguments
  Input params;
  {
    axom::CLI::App app {"Naive triangle mesh intersection tester"};
    try
    {
      params.parse(argc, argv, app);
    }
    catch(const axom::CLI::ParseError& e)
    {
      return app.exit(e);
    }
  }

  // Update the logging level based on verbosity flag
  axom::slic::setLoggingMsgLevel(params.isVerbose() ? axom::slic::message::Debug
                                                    : axom::slic::message::Info);

  // Load STL mesh into local TriangleMesh struct
  SLIC_INFO(axom::fmt::format("Reading file: '{}'...\n", params.mesh_file));
  TriangleMesh mesh = makeTriangleMesh(params.mesh_file, params.weldThreshold);

  // Check for self-intersections; results are returned as an array of index pairs
  axom::Array<IndexPair> intersectionPairs;
  axom::utilities::Timer timer(true);
  switch(params.policy)
  {
#ifdef AXOM_USE_OPENMP
  case RuntimePolicy::omp:
    intersectionPairs =
      findIntersectionsBVH<axom::OMP_EXEC>(mesh, params.intersectionThreshold, params.isVerbose());
    break;
#endif
#ifdef AXOM_USE_CUDA
  case RuntimePolicy::cuda:
    intersectionPairs = findIntersectionsBVH<axom::CUDA_EXEC<256>>(mesh,
                                                                   params.intersectionThreshold,
                                                                   params.isVerbose());
    break;
#endif
#ifdef AXOM_USE_HIP
  case RuntimePolicy::hip:
    intersectionPairs = findIntersectionsBVH<axom::HIP_EXEC<256>>(mesh,
                                                                  params.intersectionThreshold,
                                                                  params.isVerbose());
    break;
#endif
  default:  // RuntimePolicy::seq
    intersectionPairs =
      findIntersectionsBVH<axom::SEQ_EXEC>(mesh, params.intersectionThreshold, params.isVerbose());
    break;
  }
  timer.stop();

  SLIC_INFO(axom::fmt::format("Computing intersections {} took {:4.3} seconds.",
                              "with a BVH tree",
                              timer.elapsedTimeInSec()));
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Mesh had {:L} intersection pairs",
                              intersectionPairs.size()));

  // print first few pairs
  const int numIntersections = intersectionPairs.size();
  if(numIntersections > 0 && params.isVerbose())
  {
    constexpr int MAX_PRINT = 20;
    if(numIntersections > MAX_PRINT)
    {
      intersectionPairs.resize(MAX_PRINT);
      SLIC_INFO(axom::fmt::format("First {} intersection pairs: {} ...\n",
                                  MAX_PRINT,
                                  axom::fmt::join(intersectionPairs, ", ")));
    }
    else
    {
      SLIC_INFO(
        axom::fmt::format("Intersection pairs: {}\n", axom::fmt::join(intersectionPairs, ", ")));
    }
  }

  return 0;
}
