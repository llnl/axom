// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: quest_candidates_examples.cpp
///
/// This example takes as input two Blueprint unstructured hex meshes, and
/// finds the candidates of intersection between the meshes using a
/// spatial index, either a Bounding Volume Hierarchy or an Implicit Grid.
/// The example supports HIP and CUDA execution through RAJA.
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#ifdef AXOM_USE_UMPIRE
  #include "umpire/Umpire.hpp"
#else
  #error This example requires axom to be configured with Umpire support
#endif

#ifdef AXOM_USE_CONDUIT
  #include "conduit_relay.hpp"
  #include "conduit_blueprint.hpp"
#else
  #error This example requires axom to be configured with Conduit support
#endif

#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <memory>
#include <array>

using seq_exec = axom::SEQ_EXEC;

using UMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
using IndexPair = std::pair<axom::IndexType, axom::IndexType>;
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
    const std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";

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
  static const std::set<std::string> s_validMethods;

  std::string mesh_file_first {""};
  std::string mesh_file_second {""};
  std::string method {"bvh"};
  int resolution {0};
  bool verboseOutput {false};
  RuntimePolicy policy {RuntimePolicy::seq};
  std::string annotationMode {"none"};

  void parse(int argc, char** argv, axom::CLI::App& app);
  bool isVerbose() const { return verboseOutput; }
};

void Input::parse(int argc, char** argv, axom::CLI::App& app)
{
  app.add_option("-i, --infile", mesh_file_first)
    ->description("The first input Blueprint mesh file to insert into spatial index")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_option("-q, --queryfile", mesh_file_second)
    ->description("The second input Blueprint mesh file to query spatial index")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app
    .add_option("-r,--resolution",
                resolution,
                "With '-m implicit', set resolution of implicit grid. \n"
                "Set to less than 1 to use the implicit grid spatial index\n"
                "with a resolution of the cube root of the number of\n"
                "hexes.")
    ->capture_default_str();

  app.add_flag("-v,--verbose", verboseOutput)
    ->description("Increase logging verbosity")
    ->capture_default_str();

  app.add_option("-p, --policy", policy)
    ->description(
      "Execution policy."
      "\nSet to 'seq' or 0 to use the RAJA sequential policy."
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
      "\nSet to 'omp' or 1 to use the RAJA openmp policy."
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
      "\nSet to 'cuda' or 2 to use the RAJA cuda policy."
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
      "\nSet to 'hip' or 3 to use the RAJA hip policy."
#endif
      )
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

  app
    .add_option("-m, --method",
                method,
                "Method to use. \n"
                "Set to 'bvh' to use the bounding volume hierarchy spatial index.\n"
                "Set to 'implicit' to use the implicit grid spatial index.")
    ->capture_default_str()
    ->check(axom::CLI::IsMember {Input::s_validMethods});

#ifdef AXOM_USE_CALIPER
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif

  app.get_formatter()->column_width(40);

  app.parse(argc, argv);  // Could throw an exception

  // Output parsed information
  SLIC_INFO(axom::fmt::format(
    R"(
     Parsed parameters:
      * First Blueprint mesh to insert: '{}'
      * Second Blueprint mesh to query: '{}'
      * Verbose logging: {}
      * Spatial method: '{}'
      * Resolution: '{}'
      * Runtime execution policy: '{}'
      )",
    mesh_file_first,
    mesh_file_second,
    verboseOutput,
    method == "bvh" ? "Bounding Volume Hierarchy (BVH)" : "Implicit Grid",
    method == "bvh" ? "Not Applicable" : std::to_string(resolution),
    policy ==
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
        RuntimePolicy::omp
      ? "omp"
      : policy ==
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
        RuntimePolicy::cuda
      ? "cuda"
      : policy ==
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
        RuntimePolicy::hip
      ? "hip"
      : policy ==
#endif
        RuntimePolicy::seq
      ? "seq"
      : "policy not valid"));
}

const std::set<std::string> Input::s_validMethods({
  "bvh",
  "implicit",
});

//-----------------------------------------------------------------------------
/// Basic hexahedron mesh to be used in our application
//-----------------------------------------------------------------------------
struct HexMesh
{
  using Point = axom::primal::Point<double, 3>;
  using Hexahedron = axom::primal::Hexahedron<double, 3>;
  using BoundingBox = axom::primal::BoundingBox<double, 3>;

  axom::IndexType numHexes() const { return m_hexes.size(); }
  axom::Array<Hexahedron>& hexes() { return m_hexes; }
  const axom::Array<Hexahedron>& hexes() const { return m_hexes; }

  BoundingBox& meshBoundingBox() { return m_meshBoundingBox; }
  const BoundingBox& meshBoundingBox() const { return m_meshBoundingBox; }

  axom::Array<BoundingBox>& hexBoundingBoxes() { return m_hexBoundingBoxes; }
  const axom::Array<BoundingBox>& hexBoundingBoxes() const { return m_hexBoundingBoxes; }

  axom::Array<Hexahedron> m_hexes;
  axom::Array<BoundingBox> m_hexBoundingBoxes;
  BoundingBox m_meshBoundingBox;
};

template <typename ExecSpace>
HexMesh loadBlueprintHexMesh(const std::string& mesh_path, bool verboseOutput = false)
{
  AXOM_ANNOTATE_SCOPE("load Blueprint hexahedron mesh");

  using HexArray = axom::Array<typename HexMesh::Hexahedron>;
  using BBoxArray = axom::Array<typename HexMesh::BoundingBox>;
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

  // Get ids of necessary allocators
  const int host_allocator = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  const int kernel_allocator = on_device
    ? axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    : axom::execution_space<ExecSpace>::allocatorID();

  HexMesh hexMesh;

  // Load Blueprint mesh into Conduit node
  conduit::Node n_load;
  conduit::relay::io::blueprint::read_mesh(mesh_path, n_load);

  // Check if Blueprint mesh conforms
  conduit::Node n_info;
  if(conduit::blueprint::verify("mesh", n_load, n_info) == false)
  {
    n_info.print();
    SLIC_ERROR("Mesh verification has failed!");
  }

  // Verify this is a hexahedral mesh
  std::string shape = n_load[0]["topologies/topo/elements/shape"].as_string();
  if(shape != "hex")
  {
    SLIC_ERROR("A hex mesh was expected!");
  }

  const int HEX_OFFSET = 8;

  int num_nodes = (n_load[0]["coordsets/coords/values/x"]).dtype().number_of_elements();

  int connectivity_size =
    (n_load[0]["topologies/topo/elements/connectivity"]).dtype().number_of_elements();

  // Sanity check for number of cells
  int cell_calc_from_nodes = std::round(std::pow(std::pow(num_nodes, 1.0 / 3.0) - 1, 3));
  int cell_calc_from_connectivity = connectivity_size / HEX_OFFSET;
  if(cell_calc_from_nodes != cell_calc_from_connectivity)
  {
    SLIC_ERROR("Number of cells is not expected!\n"
               << "First calculation is " << cell_calc_from_nodes << " and second calculation is "
               << cell_calc_from_connectivity);
  }

  // extract hexes into an axom::Array
  auto x_vals_h = axom::ArrayView<double>(n_load[0]["coordsets/coords/values/x"].value(), num_nodes);
  auto y_vals_h = axom::ArrayView<double>(n_load[0]["coordsets/coords/values/y"].value(), num_nodes);
  auto z_vals_h = axom::ArrayView<double>(n_load[0]["coordsets/coords/values/z"].value(), num_nodes);

  // Move xyz values onto device
  axom::Array<double> x_vals_d =
    on_device ? axom::Array<double>(x_vals_h, kernel_allocator) : axom::Array<double>();
  auto x_vals_view = on_device ? x_vals_d.view() : x_vals_h;

  axom::Array<double> y_vals_d =
    on_device ? axom::Array<double>(y_vals_h, kernel_allocator) : axom::Array<double>();
  auto y_vals_view = on_device ? y_vals_d.view() : y_vals_h;

  axom::Array<double> z_vals_d =
    on_device ? axom::Array<double>(z_vals_h, kernel_allocator) : axom::Array<double>();
  auto z_vals_view = on_device ? z_vals_d.view() : z_vals_h;

  // Move connectivity information onto device
  auto connectivity_h =
    axom::ArrayView<int>(n_load[0]["topologies/topo/elements/connectivity"].value(),
                         connectivity_size);

  axom::Array<int> connectivity_d =
    on_device ? axom::Array<int>(connectivity_h, kernel_allocator) : axom::Array<int>();
  auto connectivity_view = on_device ? connectivity_d.view() : connectivity_h;

  // Initialize hex elements and bounding boxes
  const int numCells = connectivity_size / HEX_OFFSET;

  hexMesh.m_hexes = HexArray(numCells, numCells, kernel_allocator);
  auto m_hexes_v = (hexMesh.m_hexes).view();

  hexMesh.m_hexBoundingBoxes = BBoxArray(numCells, numCells, kernel_allocator);
  auto m_hexBoundingBoxes_v = (hexMesh.m_hexBoundingBoxes).view();

  axom::for_all<ExecSpace>(
    numCells,
    AXOM_LAMBDA(axom::IndexType icell) {
      HexMesh::Hexahedron hex;
      HexMesh::Point hexPoints[HEX_OFFSET];
      for(int j = 0; j < HEX_OFFSET; j++)
      {
        int offset = icell * HEX_OFFSET;
        hexPoints[j] = HexMesh::Point({x_vals_view[connectivity_view[offset + j]],
                                       y_vals_view[connectivity_view[offset + j]],
                                       z_vals_view[connectivity_view[offset + j]]});
      }
      hex = HexMesh::Hexahedron(hexPoints[0],
                                hexPoints[1],
                                hexPoints[2],
                                hexPoints[3],
                                hexPoints[4],
                                hexPoints[5],
                                hexPoints[6],
                                hexPoints[7]);
      m_hexes_v[icell] = hex;
      m_hexBoundingBoxes_v[icell] = axom::primal::compute_bounding_box(hex);
    });

  // Initialize mesh's bounding box on the host
  BBoxArray hexBoundingBoxes_h =
    on_device ? BBoxArray(hexMesh.m_hexBoundingBoxes, host_allocator) : hexMesh.m_hexBoundingBoxes;
  for(const auto& hexbb : hexBoundingBoxes_h)
  {
    hexMesh.m_meshBoundingBox.addBox(hexbb);
  }

  SLIC_INFO(axom::fmt::format("Mesh bounding box is {}.\n", hexMesh.meshBoundingBox()));

  // Optional verbose output that writes Blueprint mesh to vtk
  if(verboseOutput)
  {
    UMesh* mesh = new UMesh(3, axom::mint::HEX);

    // Append mesh nodes
    for(int i = 0; i < num_nodes; i++)
    {
      mesh->appendNode(x_vals_h[i], y_vals_h[i], z_vals_h[i]);
    }

    // Append mesh cells
    for(int i = 0; i < numCells; i++)
    {
      const axom::IndexType cell[] = {
        connectivity_h[i * HEX_OFFSET],
        connectivity_h[(i * HEX_OFFSET) + 1],
        connectivity_h[(i * HEX_OFFSET) + 2],
        connectivity_h[(i * HEX_OFFSET) + 3],
        connectivity_h[(i * HEX_OFFSET) + 4],
        connectivity_h[(i * HEX_OFFSET) + 5],
        connectivity_h[(i * HEX_OFFSET) + 6],
        connectivity_h[(i * HEX_OFFSET) + 7],
      };

      mesh->appendCell(cell);
    }

    // Write out to vtk for test viewing
    SLIC_INFO("Writing out Blueprint mesh to test.vtk for debugging...");
    axom::mint::write_vtk(mesh, "test.vtk");

    delete mesh;
    mesh = nullptr;
  }  // end of verbose output

  return hexMesh;
}  // end of loadBlueprintHexMesh

template <typename ExecSpace>
std::vector<IndexPair> findCandidatesBVH(const HexMesh& insertMesh, const HexMesh& queryMesh)
{
  AXOM_ANNOTATE_BEGIN("initializing BVH");

  SLIC_INFO("Running BVH candidates algorithm in execution Space: "
            << axom::execution_space<ExecSpace>::name());

  using IndexArray = axom::Array<axom::IndexType>;
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

  std::vector<IndexPair> candidatePairs;

  // Get ids of necessary allocators
  const int host_allocator = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  const int kernel_allocator = on_device
    ? axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    : axom::execution_space<ExecSpace>::allocatorID();

  // Get the insert and query bounding boxes
  auto insert_bbox_v = (insertMesh.hexBoundingBoxes()).view();
  auto query_bbox_v = (queryMesh.hexBoundingBoxes()).view();

  // Initialize a BVH tree over the insert mesh bounding boxes
  axom::spin::BVH<3, ExecSpace, double> bvh;
  bvh.setAllocatorID(kernel_allocator);
  bvh.initialize(insert_bbox_v, insert_bbox_v.size());
  SLIC_INFO(axom::fmt::format("0: Initialized BVH."));
  AXOM_ANNOTATE_END("initializing BVH");

  // Search for candidate bounding boxes of hexes to query;
  AXOM_ANNOTATE_BEGIN("query candidates");
  IndexArray offsets_d(query_bbox_v.size(), query_bbox_v.size(), kernel_allocator);
  IndexArray counts_d(query_bbox_v.size(), query_bbox_v.size(), kernel_allocator);
  IndexArray candidates_d(0, 0, kernel_allocator);

  auto offsets_v = offsets_d.view();
  auto counts_v = counts_d.view();
  bvh.findBoundingBoxes(offsets_v, counts_v, candidates_d, query_bbox_v.size(), query_bbox_v);

  SLIC_INFO(axom::fmt::format("1: Queried candidate bounding boxes."));
  AXOM_ANNOTATE_END("query candidates");

  AXOM_ANNOTATE_BEGIN("write candidate pairs");

  // Initialize candidate pairs on device.
  auto candidates_v = candidates_d.view();
  IndexArray firstPair_d(candidates_v.size(), candidates_v.size(), kernel_allocator);
  IndexArray secondPair_d(candidates_v.size(), candidates_v.size(), kernel_allocator);
  auto first_pair_v = firstPair_d.view();
  auto second_pair_v = secondPair_d.view();

  axom::for_all<ExecSpace>(
    query_bbox_v.size(),
    AXOM_LAMBDA(axom::IndexType icell) {
      axom::IndexType offset = offsets_v[icell];

      for(int j = 0; j < counts_v[icell]; j++)
      {
        int pair_index = offset + j;
        first_pair_v[pair_index] = icell;
        second_pair_v[pair_index] = candidates_v[pair_index];
      }
    });

  SLIC_INFO(axom::fmt::format("2: Initialized candidate pairs (on device)."));
  AXOM_ANNOTATE_END("write candidate pairs");

  // copy pairs back to host and into return array
  AXOM_ANNOTATE_BEGIN("copy pairs to host");

  IndexArray candidates_h[2] = {on_device ? IndexArray(firstPair_d, host_allocator) : IndexArray(),
                                on_device ? IndexArray(secondPair_d, host_allocator) : IndexArray()};

  auto candidate1_h_v = on_device ? candidates_h[0].view() : firstPair_d.view();
  auto candidate2_h_v = on_device ? candidates_h[1].view() : secondPair_d.view();

  candidatePairs.reserve(candidates_v.size());
  for(axom::IndexType idx = 0; idx < candidates_v.size(); ++idx)
  {
    candidatePairs.emplace_back(std::make_pair(candidate1_h_v[idx], candidate2_h_v[idx]));
  }

  SLIC_INFO(axom::fmt::format("3: Moved candidate pairs to host."));

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              R"(Stats for query
    -- Number of insert-BVH mesh hexes {:L}
    -- Number of query mesh hexes {:L}
    -- Total possible candidates {:L}
    -- Candidates from BVH query {:L}
    )",
                              insertMesh.numHexes(),
                              queryMesh.numHexes(),
                              1.0 * insertMesh.numHexes() * queryMesh.numHexes(),
                              candidatePairs.size()));
  AXOM_ANNOTATE_END("copy pairs to host");

  return candidatePairs;
}  // end of findCandidatesBVH for Blueprint Meshes

template <typename ExecSpace>
std::vector<IndexPair> findCandidatesImplicit(const HexMesh& insertMesh,
                                              const HexMesh& queryMesh,
                                              int resolution)
{
  AXOM_ANNOTATE_BEGIN("initializing implicit grid");

  std::vector<IndexPair> candidatePairs;

  SLIC_INFO("Running Implicit Grid candidates algorithm in execution Space: "
            << axom::execution_space<ExecSpace>::name());

  using IndexArray = axom::Array<int>;
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

  // Get ids of necessary allocators
  const int host_allocator = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  const int kernel_allocator = on_device
    ? axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    : axom::execution_space<ExecSpace>::allocatorID();

  // Get the insert and query bounding boxes
  auto insert_bbox_v = (insertMesh.hexBoundingBoxes()).view();
  auto query_bbox_v = (queryMesh.hexBoundingBoxes()).view();

  // Bounding box of entire insert mesh
  HexMesh::BoundingBox insert_mesh_bbox_h = insertMesh.meshBoundingBox();

  // If given resolution is less than one, use the cube root of the
  // number of hexes
  if(resolution < 1)
  {
    resolution = (int)(1 + std::pow(insertMesh.numHexes(), 1 / 3.));
  }

  const axom::primal::Point<int, 3> resolutions(resolution);

  axom::spin::ImplicitGrid<3, ExecSpace, int> gridIndex(insert_mesh_bbox_h,
                                                        &resolutions,
                                                        insertMesh.numHexes(),
                                                        kernel_allocator);
  gridIndex.insert(insertMesh.numHexes(), insert_bbox_v.data());
  SLIC_INFO(axom::fmt::format("0: Initialized Implicit Grid."));
  AXOM_ANNOTATE_END("initializing implicit grid");

  AXOM_ANNOTATE_BEGIN("query candidates");
  IndexArray offsets_d(query_bbox_v.size(), query_bbox_v.size(), kernel_allocator);
  IndexArray counts_d(query_bbox_v.size(), query_bbox_v.size(), kernel_allocator);

  auto offsets_v = offsets_d.view();
  auto counts_v = counts_d.view();

  // Object to query the candidates
  const auto grid_device = gridIndex.getQueryObject();

  axom::ReduceSum<ExecSpace, int> totalCandidatePairs(0);

  // First pass: get number of bounding box candidates for each query bounding
  // box. Logic here mirrors the quest_bvh_two_pass.cpp example. See also
  // the "Device Traversal API" for BVH.
  axom::for_all<ExecSpace>(
    queryMesh.numHexes(),
    AXOM_LAMBDA(int icell) {
      int count = 0;

      // Define a function that is called on every candidate reached during
      // traversal. The function below simply counts the number of candidates
      // that intersect with the given query bounding box.
      auto isBBIntersect = [&](int candidateIdx) {
        if(axom::primal::intersect(insert_bbox_v[candidateIdx], query_bbox_v[icell]))
        {
          count++;
        }
      };

      // Call visitCandidates to iterate through the candidates
      grid_device.visitCandidates(query_bbox_v[icell], isBBIntersect);

      // Store the number of intersections for each query bounding box
      counts_v[icell] = count;
      totalCandidatePairs += count;
    });

  SLIC_INFO(axom::fmt::format("1: Queried candidate bounding boxes."));
  AXOM_ANNOTATE_END("query candidates");

  AXOM_ANNOTATE_BEGIN("write candidate pairs");

  // Generate offsets from the counts
  // (Note: exclusive scan for offsets in candidate array
  axom::exclusive_scan<ExecSpace>(axom::ArrayView<int>(counts_v.data(), queryMesh.numHexes()),
                                  axom::ArrayView<int>(offsets_v.data(), queryMesh.numHexes()));

  // Initialize candidatePairs to return.
  // Allocate arrays for candidate pairs
  IndexArray firstPair_d(totalCandidatePairs.get(), totalCandidatePairs.get(), kernel_allocator);
  IndexArray secondPair_d(totalCandidatePairs.get(), totalCandidatePairs.get(), kernel_allocator);
  auto first_pair_v = firstPair_d.view();
  auto second_pair_v = secondPair_d.view();

  // Second pass: fill candidates array pairs on device
  axom::for_all<ExecSpace>(
    queryMesh.numHexes(),
    AXOM_LAMBDA(axom::IndexType icell) {
      axom::IndexType offset = offsets_v[icell];

      // Store the intersection candidate
      auto fillCandidates = [&](int candidateIdx) {
        if(axom::primal::intersect(insert_bbox_v[candidateIdx], query_bbox_v[icell]))
        {
          first_pair_v[offset] = icell;
          second_pair_v[offset] = candidateIdx;
          offset++;
        }
      };

      grid_device.visitCandidates(query_bbox_v[icell], fillCandidates);
    });

  SLIC_INFO(axom::fmt::format("2: Initialized candidate pairs (on device)."));
  AXOM_ANNOTATE_END("write candidate pairs");

  // copy results back to host and into return vector
  AXOM_ANNOTATE_BEGIN("copy pairs to host");

  IndexArray candidates_h[2] = {on_device ? IndexArray(firstPair_d, host_allocator) : IndexArray(),
                                on_device ? IndexArray(secondPair_d, host_allocator) : IndexArray()};

  auto candidate1_h_v = on_device ? candidates_h[0].view() : firstPair_d.view();
  auto candidate2_h_v = on_device ? candidates_h[1].view() : secondPair_d.view();

  candidatePairs.reserve(totalCandidatePairs.get());
  for(axom::IndexType idx = 0; idx < totalCandidatePairs.get(); ++idx)
  {
    candidatePairs.emplace_back(std::make_pair(candidate1_h_v[idx], candidate2_h_v[idx]));
  }

  SLIC_INFO(axom::fmt::format("3: Moved candidate pairs to host."));

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              R"(Stats for query
    -- Number of insert mesh hexes {:L}
    -- Number of query mesh hexes {:L}
    -- Total possible candidates {:L}
    -- Candidates from Implicit Grid query {:L}
    )",
                              insertMesh.numHexes(),
                              queryMesh.numHexes(),
                              1.0 * insertMesh.numHexes() * queryMesh.numHexes(),
                              candidatePairs.size()));
  AXOM_ANNOTATE_END("copy pairs to host");

  return candidatePairs;
}

int main(int argc, char** argv)
{
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  using omp_exec = axom::OMP_EXEC;
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  constexpr int CUDA_BLK_SZ = 256;
  using cuda_exec = axom::CUDA_EXEC<CUDA_BLK_SZ>;
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  constexpr int HIP_BLK_SZ = 256;
  using hip_exec = axom::HIP_EXEC<HIP_BLK_SZ>;
#endif

  // Initialize logger; use RAII so it will finalize at the end of the application
  BasicLogger logger;

  // Parse the command line arguments
  Input params;
  {
    axom::CLI::App app {"Blueprint Hex BVH mesh candidate tester"};
    try
    {
      params.parse(argc, argv, app);
    }
    catch(const axom::CLI::ParseError& e)
    {
      return app.exit(e);
    }
  }

  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(params.annotationMode);
  AXOM_ANNOTATE_SCOPE("quest candidates example");

  AXOM_ANNOTATE_BEGIN("load Blueprint meshes");

  // Update the logging level based on verbosity flag
  axom::slic::setLoggingMsgLevel(params.isVerbose() ? axom::slic::message::Debug
                                                    : axom::slic::message::Info);

  // Load Blueprint mesh to insert into spatial index
  SLIC_INFO(axom::fmt::format("Reading Blueprint file to insert: '{}'...\n", params.mesh_file_first));

  HexMesh insert_mesh;

  switch(params.policy)
  {
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  case RuntimePolicy::omp:
    insert_mesh = loadBlueprintHexMesh<omp_exec>(params.mesh_file_first, params.isVerbose());
    break;
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  case RuntimePolicy::cuda:
    insert_mesh = loadBlueprintHexMesh<cuda_exec>(params.mesh_file_first, params.isVerbose());
    break;
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  case RuntimePolicy::hip:
    insert_mesh = loadBlueprintHexMesh<hip_exec>(params.mesh_file_first, params.isVerbose());
    break;
#endif
  default:  // RuntimePolicy::seq
    insert_mesh = loadBlueprintHexMesh<seq_exec>(params.mesh_file_first, params.isVerbose());
    break;
  }

  // Load Blueprint mesh for querying spatial index
  SLIC_INFO(axom::fmt::format("Reading Blueprint file to query: '{}'...\n", params.mesh_file_second));

  HexMesh query_mesh;

  switch(params.policy)
  {
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  case RuntimePolicy::omp:
    query_mesh = loadBlueprintHexMesh<omp_exec>(params.mesh_file_second, params.isVerbose());
    break;
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  case RuntimePolicy::cuda:
    query_mesh = loadBlueprintHexMesh<cuda_exec>(params.mesh_file_second, params.isVerbose());
    break;
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  case RuntimePolicy::hip:
    query_mesh = loadBlueprintHexMesh<hip_exec>(params.mesh_file_second, params.isVerbose());
    break;
#endif
  default:  // RuntimePolicy::seq
    query_mesh = loadBlueprintHexMesh<seq_exec>(params.mesh_file_second, params.isVerbose());
    break;
  }

  SLIC_INFO(axom::fmt::format("Finished reading in Blueprint files."));
  AXOM_ANNOTATE_END("load Blueprint meshes");

  AXOM_ANNOTATE_BEGIN("find candidates");
  // Check for candidates; results are returned as an array of index pairs
  std::vector<IndexPair> candidatePairs;

  if(params.method == "bvh")
  {
    switch(params.policy)
    {
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
    case RuntimePolicy::omp:
      candidatePairs = findCandidatesBVH<omp_exec>(insert_mesh, query_mesh);
      break;
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
    case RuntimePolicy::cuda:
      candidatePairs = findCandidatesBVH<cuda_exec>(insert_mesh, query_mesh);
      break;
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
    case RuntimePolicy::hip:
      candidatePairs = findCandidatesBVH<hip_exec>(insert_mesh, query_mesh);
      break;
#endif
    default:  // RuntimePolicy::seq
      candidatePairs = findCandidatesBVH<seq_exec>(insert_mesh, query_mesh);
      break;
    }
  }
  // Implicit Grid
  else
  {
    switch(params.policy)
    {
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
    case RuntimePolicy::omp:
      candidatePairs = findCandidatesImplicit<omp_exec>(insert_mesh, query_mesh, params.resolution);
      break;
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
    case RuntimePolicy::cuda:
      candidatePairs = findCandidatesImplicit<cuda_exec>(insert_mesh, query_mesh, params.resolution);
      break;
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
    case RuntimePolicy::hip:
      candidatePairs = findCandidatesImplicit<hip_exec>(insert_mesh, query_mesh, params.resolution);
      break;
#endif
    default:  // RuntimePolicy::seq
      candidatePairs = findCandidatesImplicit<seq_exec>(insert_mesh, query_mesh, params.resolution);
      break;
    }
  }

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Mesh had {:L} candidates pairs",
                              candidatePairs.size()));

  // print first few pairs
  const int numCandidates = candidatePairs.size();
  if(numCandidates > 0 && params.isVerbose())
  {
    constexpr int MAX_PRINT = 20;
    if(numCandidates > MAX_PRINT)
    {
      candidatePairs.resize(MAX_PRINT);
      SLIC_INFO(axom::fmt::format("First {} candidate pairs: {} ...\n",
                                  MAX_PRINT,
                                  axom::fmt::join(candidatePairs, ", ")));
    }
    else
    {
      SLIC_INFO(axom::fmt::format("Candidate pairs: {}\n", axom::fmt::join(candidatePairs, ", ")));
    }

    // Write out candidate pairs
    SLIC_INFO("Writing out candidate pairs...");
    std::ofstream outf("candidates.txt");

    outf << candidatePairs.size() << " candidate pairs:" << std::endl;
    for(unsigned long i = 0; i < candidatePairs.size(); ++i)
    {
      outf << candidatePairs[i].first << " " << candidatePairs[i].second << std::endl;
    }
  }
  AXOM_ANNOTATE_END("find candidates");

  return 0;
}
