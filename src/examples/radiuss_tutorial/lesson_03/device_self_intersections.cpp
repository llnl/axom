// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: device_self_intersections.cpp
///
/// This example ports the naive self-intersection algorithm to additional
/// execution and memory spaces w/ the help of RAJA and Umpire
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
    const std::string slicFormatStr = "[lesson_03: <LEVEL>] <MESSAGE> \n";

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
  bool useBoundingBoxes {false};
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

  app.add_flag("--use-bounding-boxes", useBoundingBoxes)
    ->description("Use bounding boxes to accelerate intersection query?")
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
      * Use bounding boxes to accelerate query: {}
      * Runtime execution policy: '{}'
      )",
    mesh_file,
    weldThreshold,
    (weldThreshold <= 0.),
    intersectionThreshold,
    verboseOutput,
    useBoundingBoxes,
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

  axom::IndexType numDegenerateTriangles() const
  {
    return static_cast<axom::IndexType>(
      std::count_if(m_degeneracies.begin(), m_degeneracies.end(), [](bool b) { return b == true; }));
  }

  bool isTriangleDegenerate(axom::IndexType idx) const { return m_degeneracies[idx]; }

  BoundingBox& meshBoundingBox() { return m_meshBoundingBox; }
  const BoundingBox& meshBoundingBox() const { return m_meshBoundingBox; }

  axom::Array<BoundingBox>& triangleBoundingBoxes() { return m_triangleBoundingBoxes; }
  const axom::Array<BoundingBox>& triangleBoundingBoxes() const { return m_triangleBoundingBoxes; }

  axom::Array<Triangle> m_triangles;
  axom::Array<bool> m_degeneracies;
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

  // keep track of degenerate triangles
  triMesh.m_degeneracies.reserve(numCells);
  for(int i = 0; i < numCells; ++i)
  {
    const bool is_degenerate = triMesh.m_triangles[i].degenerate();
    triMesh.m_degeneracies.push_back(is_degenerate);
  }

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Mesh has {:L} degenerate triangles.",
                              triMesh.numDegenerateTriangles()));

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
axom::Array<IndexPair> naiveFindIntersections(const TriangleMesh& triMesh,
                                              double tol,
                                              const bool useBoundingBoxes,
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

  // copy indices of non-degenerate triangles into a new array
  // we expect this to be fast enough to perform on the host
  // valid_v will be a view to the array in the proper space
  axom::utilities::Timer timer;
  timer.start();
  const auto totalTriangles = triMesh.numTriangles();
  IndexArray valid_h(0, totalTriangles, host_allocator);
  for(axom::IndexType idx = 0; idx < totalTriangles; ++idx)
  {
    if(!triMesh.isTriangleDegenerate(idx))
    {
      valid_h.push_back(idx);
    }
  }
  auto valid_d = on_device ? IndexArray(valid_h, kernel_allocator) : IndexArray();
  auto valid_v = on_device ? valid_d.view() : valid_h.view();

  timer.stop();
  SLIC_INFO_IF(verboseOutput,
               axom::fmt::format("Computing non-degenerate triangles took {:4.3} seconds.",
                                 timer.elapsedTimeInSec()));

  // lambda to check if two triangles w/ given indices intersect
  auto trianglesIntersect = AXOM_LAMBDA(axom::IndexType idx1, axom::IndexType idx2)
  {
    constexpr bool includeBoundaries = false;  // only use triangle interiors

    return (!useBoundingBoxes || axom::primal::intersect(bbox_v[idx1], bbox_v[idx2])) &&
      axom::primal::intersect(tris_v[idx1], tris_v[idx2], includeBoundaries, tol);
  };

  // Phase I: Find the number of intersecting pairs of triangles
  const auto validCount = valid_v.size();
  RAJA::RangeSegment row_range(0, validCount);
  RAJA::RangeSegment col_range(0, validCount);

  using KERNEL_POL = typename axom::internal::nested_for_exec<ExecSpace>::loop2d_policy;

  axom::ReduceSum<ExecSpace, int> numIntersect(0);

  // Compute the number of intersections
  timer.start();
  RAJA::kernel<KERNEL_POL>(
    RAJA::make_tuple(col_range, row_range),
    AXOM_LAMBDA(int col, int row) {
      if(row < col && trianglesIntersect(valid_v[row], valid_v[col]))
      {
        numIntersect += 1;
      }
    });
  timer.stop();
  SLIC_INFO_IF(
    verboseOutput,
    axom::fmt::format("Finding intersection count took {:4.3} seconds.", timer.elapsedTimeInSec()));

  // Phase II: If there are intersections, we need to extract their indices to an array
  if(numIntersect.get() > 0)
  {
    const auto numIndices = 2 * numIntersect.get();

    // allocate an array for the intersections (without initializing it)
    IndexArray intersections_d(axom::ArrayOptions::Uninitialized {},
                               numIndices,
                               numIndices,
                               kernel_allocator);

    // allocate a counter; we'll use atomic operations to increment it
    auto counter_d = IndexArray(1, 1, kernel_allocator);
    counter_d.fill(0);
    auto* counter_p = counter_d.data();

    // populate array with intersections
    timer.start();
    auto intersections_v = intersections_d.view();
    RAJA::kernel<KERNEL_POL>(
      RAJA::make_tuple(col_range, row_range),
      AXOM_LAMBDA(int col, int row) {
        if(row < col && trianglesIntersect(valid_v[row], valid_v[col]))
        {
          const auto idx = axom::atomicAdd<axom::auto_atomic>(counter_p, axom::IndexType {2});
          intersections_v[idx + 0] = valid_v[row];
          intersections_v[idx + 1] = valid_v[col];
        }
      });
    timer.stop();
    SLIC_INFO_IF(verboseOutput,
                 axom::fmt::format("Inserting intersection pairs took {:4.3} seconds.",
                                   timer.elapsedTimeInSec()));

    // copy intersections back to host, if necessary
    auto intersections_h = on_device ? IndexArray(intersections_d, host_allocator) : IndexArray();
    // in either event, intersections_h_v will be a valid host array
    auto interections_h_v = on_device ? intersections_h.view() : intersections_d.view();

    // copy the results into the return vector
    timer.start();
    for(axom::IndexType idx = 0; idx < numIndices; idx += 2)
    {
      intersectionPairs.emplace_back(std::make_pair(interections_h_v[idx], interections_h_v[idx + 1]));
    }
    timer.stop();
    SLIC_INFO_IF(
      verboseOutput,
      axom::fmt::format("Copying back to array took {:4.3} seconds.", timer.elapsedTimeInSec()));
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
    intersectionPairs = naiveFindIntersections<axom::OMP_EXEC>(mesh,
                                                               params.intersectionThreshold,
                                                               params.useBoundingBoxes,
                                                               params.isVerbose());
    break;
#endif
#ifdef AXOM_USE_CUDA
  case RuntimePolicy::cuda:
    intersectionPairs = naiveFindIntersections<axom::CUDA_EXEC<256>>(mesh,
                                                                     params.intersectionThreshold,
                                                                     params.useBoundingBoxes,
                                                                     params.isVerbose());
    break;
#endif
#ifdef AXOM_USE_HIP
  case RuntimePolicy::hip:
    intersectionPairs = naiveFindIntersections<axom::HIP_EXEC<256>>(mesh,
                                                                    params.intersectionThreshold,
                                                                    params.useBoundingBoxes,
                                                                    params.isVerbose());
    break;
#endif
  default:  // RuntimePolicy::seq
    intersectionPairs = naiveFindIntersections<axom::SEQ_EXEC>(mesh,
                                                               params.intersectionThreshold,
                                                               params.useBoundingBoxes,
                                                               params.isVerbose());
    break;
  }
  timer.stop();

  SLIC_INFO(
    axom::fmt::format("Computing intersections {} took {:4.3} seconds.",
                      params.useBoundingBoxes ? "with bounding boxes" : "without bounding boxes",
                      timer.elapsedTimeInSec()));
  SLIC_INFO(axom::fmt::format("Mesh had {} intersection pairs", intersectionPairs.size()));

  SLIC_INFO_IF(
    intersectionPairs.size() > 0 && params.isVerbose(),
    axom::fmt::format("Intersecting pairs: {}\n", axom::fmt::join(intersectionPairs, ", ")));

  return 0;
}
