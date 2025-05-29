// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file spin_bvh_two_pass.cpp
 *  \brief This example code demonstrates the device kernel query capabilities
 *   of the BVH class in Axom's spin component.
 *
 *  This file shows how to run customized traversal queries on a bounding volume
 *  hierarchy with a device kernel-copyable object, by demonstrating its use to
 *  generate a set of candidate mesh self-intersections.
 */

// Axom includes
#include "axom/core.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/slic.hpp"
#include "axom/quest.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

namespace mint = axom::mint;
namespace primal = axom::primal;
namespace spin = axom::spin;
namespace slic = axom::slic;

using IndexType = axom::IndexType;
using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

enum class ExecPolicy
{
  CPU,
  OpenMP,
  CUDA,
  HIP
};

const std::map<std::string, ExecPolicy> validExecPolicies {{"seq", ExecPolicy::CPU},
#ifdef AXOM_USE_OPENMP
                                                           {"omp", ExecPolicy::OpenMP},
#endif
#ifdef AXOM_USE_CUDA
                                                           {"cuda", ExecPolicy::CUDA}
#endif
#ifdef AXOM_USE_HIP
                                                           {"hip", ExecPolicy::HIP}
#endif
};

//------------------------------------------------------------------------------
void initialize_logger()
{
  // initialize logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  // setup the logstreams
  std::string fmt = "";
  slic::LogStream* logStream = nullptr;

  fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);

  // register stream objects with the logger
  slic::addStreamToAllMsgLevels(logStream);
}

//------------------------------------------------------------------------------
void finalize_logger()
{
  slic::flushStreams();
  slic::finalize();
}

/*!
 * \brief Runs the "broad-phase" part of a collision detection query on a mesh.
 *
 *  In most implementations of collision detection, the search is divided into
 *  two separate phases: a "broad phase" which utilizes a fast-but-conservative
 *  collision test to quickly prune the search space of potential collisions,
 *  followed by a "narrow phase" which uses a more accurate collision test on
 *  the remaining candidate collisions.
 *
 *  The below implements the broad-phase query with a bounding volume hierarchy
 *  of bounding boxes constructed from the input mesh's zones, and returns an
 *  output of potentially-intersecting mesh elements as pairs of indices.
 *
 *  The subsequent narrow-phase implementation can then be run on the returned
 *  pairs of indices. What that narrow-phase check constitutes depends on the
 *  mesh elements; for example, triangle-triangle tests if run on a triangle
 *  surface mesh, or polyhedral intersection tests for a mesh of arbitrary
 *  polyhedral zones.
 *
 * \param [in] mesh the input mesh
 * \param [out] firstPair first index of pairs of potentially colliding zones
 * \param [out] secondPair second index of pairs of potentially colliding zones
 */
template <typename ExecSpace>
void find_collisions_broadphase(const mint::Mesh* mesh,
                                axom::Array<IndexType>& firstPair,
                                axom::Array<IndexType>& secondPair)
{
  using PointType = axom::primal::Point<double, 3>;
  using BoxType = axom::primal::BoundingBox<double, 3>;

  const int allocatorId = axom::execution_space<ExecSpace>::allocatorID();

  const int ncells = mesh->getNumberOfCells();

  axom::Array<BoxType> aabbs(ncells, ncells, allocatorId);
  const auto v_aabbs = aabbs.view();

  // Initialize the bounding box for each cell
  mint::for_all_cells<ExecSpace, mint::xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx, axom::numerics::Matrix<double> & coords, const IndexType* nodeIds) {
      AXOM_UNUSED_VAR(nodeIds);
      int numNodes = coords.getNumColumns();
      BoxType aabb;

      for(IndexType inode = 0; inode < numNodes; ++inode)
      {
        const double* node = coords.getColumn(inode);
        PointType vtx {node[mint::X_COORDINATE], node[mint::Y_COORDINATE], node[mint::Z_COORDINATE]};
        aabb.addPoint(vtx);
      }  // END for all cells nodes

      v_aabbs[cellIdx] = aabb;
    });

  // _bvh_traverse_init_start
  // Initialize BVH
  spin::BVH<3, ExecSpace, double> bvh;
  bvh.setAllocatorID(allocatorId);
  bvh.initialize(v_aabbs, v_aabbs.size());

  // Create a traverser object from the BVH
  const auto bvh_device = bvh.getTraverser();
  // _bvh_traverse_init_end

  // Define a simple intersection predicate for traversing the BVH.
  // We traverse any nodes whose bounding boxes intersect the query
  // bounding box.
  // _bvh_traverse_predicate_start
  auto bbIsect = [] AXOM_HOST_DEVICE(const BoxType& queryBbox, const BoxType& bvhBbox) -> bool {
    return queryBbox.intersectsWith(bvhBbox);
  };
  // _bvh_traverse_predicate_end

  axom::ReduceSum<ExecSpace, IndexType> total_count_reduce(0);

  // _bvh_traverse_first_pass_start
  // Allocate arrays for offsets and counts
  axom::Array<IndexType> offsets(ncells, ncells, allocatorId);
  axom::Array<IndexType> counts(ncells, ncells, allocatorId);
  const auto v_counts = counts.view();
  const auto v_offsets = offsets.view();

  // First pass: get number of bounding box collisions for each surface element
  axom::for_all<ExecSpace>(
    ncells,
    AXOM_LAMBDA(IndexType icell) {
      IndexType count = 0;

      // Define a function that is called on every leaf node reached during
      // traversal. The function below simply counts the number of candidate
      // collisions with the given query object.
      auto countCollisions = [&](std::int32_t currentNode, const std::int32_t* leafNodes) {
        AXOM_UNUSED_VAR(leafNodes);
        if(currentNode > icell)
        {
          count++;
        }
      };

      // Call traverse_tree() to run the counting query.
      bvh_device.traverse_tree(v_aabbs[icell], countCollisions, bbIsect);

      // Afterwards, we can store the number of collisions for each surface
      // element, as well as an overall count of intersections.
      v_counts[icell] = count;
      total_count_reduce += count;
    });

  // Generate offsets
  axom::exclusive_scan<ExecSpace>(counts, offsets);

  // _bvh_traverse_first_pass_end

  IndexType ncollisions = total_count_reduce.get();

  SLIC_INFO(axom::fmt::format("Found {} candidate collisions.", ncollisions));

  // _bvh_traverse_second_pass_start
  // Allocate arrays for intersection pairs
  firstPair = axom::Array<IndexType>(ncollisions, ncollisions, allocatorId);
  secondPair = axom::Array<IndexType>(ncollisions, ncollisions, allocatorId);
  const auto v_first_pair = firstPair.view();
  const auto v_second_pair = secondPair.view();

  // Second pass: fill broad-phase collisions array
  axom::for_all<ExecSpace>(
    ncells,
    AXOM_LAMBDA(IndexType icell) {
      IndexType offset = v_offsets[icell];

      // Define a leaf node function that stores the intersection candidate.
      auto fillCollisions = [&](std::int32_t currentNode, const std::int32_t* leafs) {
        if(currentNode > icell)
        {
          v_first_pair[offset] = icell;
          v_second_pair[offset] = leafs[currentNode];
          offset++;
        }
      };

      // Call traverse_tree() a second time to run the counting query.
      bvh_device.traverse_tree(v_aabbs[icell], fillCollisions, bbIsect);
    });
  // _bvh_traverse_second_pass_end
}

/*!
 * \brief Runs the "narrow-phase" part of a collision detection query on a mesh,
 *  given a set of candidate pairs from the broad-phase search.
 *
 *  This assumes a surface mesh of triangular elements.
 *
 * \param [in] mesh the input mesh
 * \param [in] inFirstPair first index of pairs of potentially colliding zones
 * \param [in] inSecondPair second index of pairs of potentially colliding zones
 * \param [in] outFirstPair first index of pairs of actual colliding zones
 * \param [in] outSecondPair second index of pairs of actual colliding zones
 */
template <typename ExecSpace>
void find_collisions_narrowphase(const mint::Mesh* mesh,
                                 axom::Array<IndexType>& inFirstPair,
                                 axom::Array<IndexType>& inSecondPair,
                                 axom::Array<IndexType>& outFirstPair,
                                 axom::Array<IndexType>& outSecondPair)
{
  using PointType = axom::primal::Point<double, 3>;
  using TriangleType = axom::primal::Triangle<double, 3>;

  int allocatorId = axom::execution_space<ExecSpace>::allocatorID();

  const int ncells = mesh->getNumberOfCells();

  axom::Array<TriangleType> triangles(ncells, ncells, allocatorId);
  const auto v_triangles = triangles.view();

  // Create an array with our surface mesh's triangles
  mint::for_all_cells<ExecSpace, mint::xargs::coords>(
    mesh,
    AXOM_LAMBDA(IndexType cellIdx, axom::numerics::Matrix<double> & coords, const IndexType* nodeIds) {
      AXOM_UNUSED_VAR(nodeIds);
      TriangleType tri;

      int numNodes = coords.getNumColumns();
      SLIC_ASSERT(numNodes == 3);
      for(IndexType inode = 0; inode < numNodes; ++inode)
      {
        const double* node = coords.getColumn(inode);
        PointType vtx {node[mint::X_COORDINATE], node[mint::Y_COORDINATE], node[mint::Z_COORDINATE]};
        tri[inode] = vtx;
      }  // END for all cells nodes

      v_triangles[cellIdx] = tri;
    });

  // Allocate our output arrays - we'll set the initial capacity to the number
  // of candidates.
  IndexType ncandidates = inFirstPair.size();
  outFirstPair = axom::Array<IndexType>(ncandidates, ncandidates, allocatorId);
  outSecondPair = axom::Array<IndexType>(ncandidates, ncandidates, allocatorId);

  axom::Array<IndexType> counter(1, 1, allocatorId);

  axom::ReduceSum<ExecSpace, IndexType> total_count_reduce(0);

  const auto v_counter = counter.view();
  const auto v_inFirstPair = inFirstPair.view();
  const auto v_inSecondPair = inSecondPair.view();
  const auto v_outFirstPair = outFirstPair.view();
  const auto v_outSecondPair = outSecondPair.view();

  // For each candidate pair, we'll run a triangle-triangle intersection check.
  axom::for_all<ExecSpace>(
    ncandidates,
    AXOM_LAMBDA(IndexType idx) {
      IndexType firstIdx = v_inFirstPair[idx];
      IndexType secondIdx = v_inSecondPair[idx];
      if(primal::intersect(v_triangles[firstIdx], v_triangles[secondIdx], false))
      {
        auto outIdx = axom::atomicAdd<ExecSpace>(&v_counter[0], IndexType {1});
        // Store actually-intersecting triangle pairs sequentially in output
        // array.
        v_outFirstPair[outIdx] = firstIdx;
        v_outSecondPair[outIdx] = secondIdx;
        total_count_reduce += 1;
      }
    });

  IndexType numValidIsects = total_count_reduce.get();

  SLIC_INFO(axom::fmt::format("Found {} actual collisions.", numValidIsects));

  // Resize to the total number of intersecting triangles.
  outFirstPair.resize(numValidIsects);
  outSecondPair.resize(numValidIsects);
}

struct Arguments
{
  std::string file_name;
  ExecPolicy exec_space {ExecPolicy::CPU};

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-f,--file", this->file_name, "specifies the input mesh file")
      ->check(axom::CLI::ExistingFile)
      ->required();

    std::string pol_info = "Sets execution space of the BVH two-pass example.\n";
    pol_info += "Set to \'seq\' to use sequential execution policy.";
#ifdef AXOM_USE_OPENMP
    pol_info += "\nSet to \'omp\' to use an OpenMP execution policy.";
#endif
#ifdef AXOM_USE_CUDA
    pol_info += "\nSet to \'cuda\' to use a CUDA GPU execution policy.";
#endif
#ifdef AXOM_USE_HIP
    pol_info += "\nSet to \'hip\' to use a HIP GPU execution policy.";
#endif
    app.add_option("-e, --exec_space", this->exec_space, pol_info)
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(validExecPolicies));

    app.get_formatter()->column_width(40);

    // could throw an exception
    app.parse(argc, argv);

    slic::flushStreams();
  }
};

int main(int argc, char** argv)
{
  initialize_logger();
  Arguments args;
  axom::CLI::App app {"Demo of the BVH 2-pass algorithm"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    retval = app.exit(e);
    finalize_logger();
    return retval;
  }
  std::unique_ptr<UMesh> surface_mesh;

  // Read file
  SLIC_INFO("Reading file: '" << args.file_name << "'...\n");
  {
    axom::quest::STLReader reader;
    reader.setFileName(args.file_name);
    reader.read();

    // Get surface mesh
    surface_mesh.reset(new UMesh(3, mint::TRIANGLE));
    reader.getMesh(surface_mesh.get());
  }

  SLIC_INFO("Mesh has " << surface_mesh->getNumberOfNodes() << " vertices and "
                        << surface_mesh->getNumberOfCells() << " triangles.");

  axom::Array<IndexType> candFirstPair, candSecondPair;
  axom::Array<IndexType> firstPair, secondPair;

  switch(args.exec_space)
  {
  case ExecPolicy::CPU:
    find_collisions_broadphase<axom::SEQ_EXEC>(surface_mesh.get(), candFirstPair, candSecondPair);
    find_collisions_narrowphase<axom::SEQ_EXEC>(surface_mesh.get(),
                                                candFirstPair,
                                                candSecondPair,
                                                firstPair,
                                                secondPair);
    break;
#if defined(AXOM_USE_RAJA)
  #if defined(AXOM_USE_OPENMP)
  case ExecPolicy::OpenMP:
    find_collisions_broadphase<axom::OMP_EXEC>(surface_mesh.get(), candFirstPair, candSecondPair);
    find_collisions_narrowphase<axom::OMP_EXEC>(surface_mesh.get(),
                                                candFirstPair,
                                                candSecondPair,
                                                firstPair,
                                                secondPair);
    break;
  #endif
  #if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  case ExecPolicy::CUDA:
    find_collisions_broadphase<axom::CUDA_EXEC<256>>(surface_mesh.get(), candFirstPair, candSecondPair);
    find_collisions_narrowphase<axom::CUDA_EXEC<256>>(surface_mesh.get(),
                                                      candFirstPair,
                                                      candSecondPair,
                                                      firstPair,
                                                      secondPair);
    break;
  #endif
  #if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_HIP)
  case ExecPolicy::HIP:
    find_collisions_broadphase<axom::HIP_EXEC<256>>(surface_mesh.get(), candFirstPair, candSecondPair);
    find_collisions_narrowphase<axom::HIP_EXEC<256>>(surface_mesh.get(),
                                                     candFirstPair,
                                                     candSecondPair,
                                                     firstPair,
                                                     secondPair);
    break;
  #endif
#endif
  default:
    SLIC_ERROR("Unsupported execution space.");
    return 1;
  }

  finalize_logger();
}
