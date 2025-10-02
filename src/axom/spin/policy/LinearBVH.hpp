// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_POLICY_LINEARBVH_HPP_
#define AXOM_SPIN_POLICY_LINEARBVH_HPP_

// axom core includes
#include "axom/core/Types.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/AnnotationMacros.hpp"
#include "axom/core/numerics/floating_point_limits.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Vector.hpp"

// linear bvh includes
#include "axom/spin/internal/linear_bvh/RadixTree.hpp"
#include "axom/spin/internal/linear_bvh/build_radix_tree.hpp"
#include "axom/spin/internal/linear_bvh/bvh_traverse.hpp"
#include "axom/spin/internal/linear_bvh/bvh_vtkio.hpp"

// C/C++ includes
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace axom
{
namespace spin
{
namespace policy
{
namespace lbvh = internal::linear_bvh;

/*
 * \brief Interface for traversing a BVH tree.
 *
 * \brief Traverse a tree to perform user-specified actions at the
 * leaves, while limiting the search to branches satisfying a
 * user-provided predicate.
 *
 * To initiate traversals, use \a traverse_tree.  It requires
 * -# an action functor to call at the leaves,
 * -# a predicate functor to decide whether to descend a branch and
 * -# some data to pass to the functors
 *
 * Both functors should only access memory that's available to the
 * execution space.  For example, GPU execution should access only
 * device and unified memory.
 *
 * \see internal::linear_bvh::bvh_traverse
 */
template <typename FloatType, int NDIMS>
class LinearBVHTraverser
{
public:
  using BoxType = primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;

  LinearBVHTraverser(axom::ArrayView<const BoxType> bboxes,
                     axom::ArrayView<const std::int32_t> inner_node_children,
                     axom::ArrayView<const std::int32_t> leaf_nodes)
    : m_inner_nodes(bboxes)
    , m_inner_node_children(inner_node_children)
    , m_leaf_nodes(leaf_nodes)
  { }

  template <typename LeafAction, typename Predicate>
  AXOM_HOST_DEVICE void traverse_tree(const PointType& p, LeafAction&& lf, Predicate&& predicate) const
  {
    auto traversePref = [](const BoxType& l, const BoxType& r, const PointType& p) {
      double sqDistL = primal::squared_distance(p, l.getCentroid());
      // If the right bbox is not valid, return max. Otherwise, the invalid right
      // bbox might actually win when we should ignore it.
      double sqDistR = r.isValid() ? primal::squared_distance(p, r.getCentroid())
                                   : axom::numerics::floating_point_limits<double>::max();
      return sqDistL > sqDistR;
    };

    lbvh::bvh_traverse(m_inner_nodes, m_inner_node_children, m_leaf_nodes, p, predicate, lf, traversePref);
  }

  /*
   * Functors \a lf and \a predicate should access only memory compatible
   * with the execution space.  For example, GPU execution should access
   * only device and unified memory.
   */
  template <typename Primitive, typename LeafAction, typename Predicate>
  AXOM_HOST_DEVICE void traverse_tree(const Primitive& p, LeafAction&& lf, Predicate&& predicate) const
  {
    auto noTraversePref = [](const BoxType& l, const BoxType& r, const Primitive& p) {
      AXOM_UNUSED_VAR(l);
      AXOM_UNUSED_VAR(r);
      AXOM_UNUSED_VAR(p);
      return false;
    };

    lbvh::bvh_traverse(m_inner_nodes, m_inner_node_children, m_leaf_nodes, p, predicate, lf, noTraversePref);
  }

private:
  axom::ArrayView<const BoxType> m_inner_nodes;  // BVH bins including leafs
  axom::ArrayView<const std::int32_t> m_inner_node_children;
  axom::ArrayView<const std::int32_t> m_leaf_nodes;  // leaf data
};

/*!
 * \brief LinearBVH provides a policy for a BVH implementation which supports
 *  parallel linear construction on both CPU and GPU.
 *
 * \note The internal node data is organized into two arrays, one containing
 *  bounding boxes of the internal nodes' two child nodes, and the other
 *  containing indices of the two child nodes (scaled by two if inner node,
 *  ones-complement if leaf node).
 *
 */
template <typename FloatType, int NDIMS, typename ExecSpace>
class LinearBVH
{
public:
  using TraverserType = LinearBVHTraverser<FloatType, NDIMS>;
  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;

  LinearBVH() = default;

  /*!
   * \brief Builds a linear BVH with the given bounding boxes as leaf nodes.
   *
   * \param [in] boxes the bounding boxes for each leaf node
   * \param [in] numBoxes the number of bounding boxes
   * \param [in] scaleFactor scale factor applied to each bounding box before insertion into the BVH
   */
  template <typename BoxIndexable>
  void buildImpl(const BoxIndexable boxes, IndexType numBoxes, FloatType scaleFactor, int allocatorID);

  /*!
   * \brief Performs a traversal to find the candidates for each query primitive.
   *
   * \param [in] predicate traversal predicate functor for bin check.
   * \param [out] offsets array of offsets into the candidate array for each query primitive
   * \param [out] counts array of candidate counts for each query primitive
   * \param [out] candidates array of the potential candidates for intersection with the BVH
   * \param [in] numObjs the number of user-supplied query primitives
   * \param [in] objs array of primitives to query against the BVH
   *
   * \return total_count the total count of candidates for all query primitives.
   */
  template <typename PrimitiveType, typename Predicate, typename PrimitiveIndexable>
  axom::Array<IndexType> findCandidatesImpl(Predicate&& predicate,
                                            const axom::ArrayView<IndexType> offsets,
                                            const axom::ArrayView<IndexType> counts,
                                            IndexType numObjs,
                                            PrimitiveIndexable objs,
                                            int allocatorID) const;

  void writeVtkFileImpl(const std::string& fileName) const;

  BoundingBoxType getBoundsImpl() const { return m_bounds; }

  TraverserType getTraverserImpl() const
  {
    return TraverserType(m_inner_nodes.view(), m_inner_node_children.view(), m_leaf_nodes.view());
  }

private:
  void allocate(std::int32_t size, int allocID)
  {
    AXOM_ANNOTATE_SCOPE("LinearBVH::allocate");
    IndexType numInnerNodes = (size - 1) * 2;
    // Need to allocate this uninitialized, since primal::BoundingBox is
    // considered non-trivially-copyable on GCC 4.9.3
    m_inner_nodes = axom::Array<BoundingBoxType>(axom::ArrayOptions::Uninitialized {},
                                                 numInnerNodes,
                                                 numInnerNodes,
                                                 allocID);
    m_inner_node_children = axom::Array<std::int32_t>(numInnerNodes, numInnerNodes, allocID);
    m_leaf_nodes = axom::Array<std::int32_t>(size, size, allocID);
  }

  bool m_initialized {false};
  axom::Array<BoundingBoxType> m_inner_nodes;  // BVH bins including leafs
  axom::Array<std::int32_t> m_inner_node_children;
  axom::Array<std::int32_t> m_leaf_nodes;  // leaf data
  primal::BoundingBox<FloatType, NDIMS> m_bounds;
};

template <typename FloatType, int NDIMS, typename ExecSpace>
template <typename BoxIndexable>
void LinearBVH<FloatType, NDIMS, ExecSpace>::buildImpl(const BoxIndexable boxes,
                                                       IndexType numBoxes,
                                                       FloatType scaleFactor,
                                                       int allocatorID)
{
  AXOM_ANNOTATE_SCOPE("LinearBVH::buildImpl");

  // STEP 1: Build a RadixTree consisting of the bounding boxes, sorted
  // by their corresponding morton code.
  lbvh::RadixTree<FloatType, NDIMS> radix_tree;
  primal::BoundingBox<FloatType, NDIMS> global_bounds;
  lbvh::build_radix_tree<ExecSpace>(boxes, numBoxes, global_bounds, radix_tree, scaleFactor, allocatorID);

  // STEP 2: emit the BVH data-structure from the radix tree
  m_bounds = global_bounds;
  allocate(numBoxes, allocatorID);

  // STEP 3: emit the BVH
  const std::int32_t size = radix_tree.m_size;
  AXOM_UNUSED_VAR(size);
  const std::int32_t inner_size = radix_tree.m_inner_size;
  SLIC_ASSERT(inner_size == size - 1);

  const auto lchildren_ptr = radix_tree.m_left_children.view();
  const auto rchildren_ptr = radix_tree.m_right_children.view();

  const auto leaf_aabb_ptr = radix_tree.m_leaf_aabbs.view();
  const auto inner_aabb_ptr = radix_tree.m_inner_aabbs.view();

  const auto bvh_inner_nodes = m_inner_nodes.view();
  const auto bvh_inner_node_children = m_inner_node_children.view();

  AXOM_ANNOTATE_BEGIN("emit_bvh_parents");
  for_all<ExecSpace>(
    inner_size,
    AXOM_LAMBDA(std::int32_t node) {
      BoundingBoxType l_aabb, r_aabb;

      std::int32_t lchild = lchildren_ptr[node];
      if(lchild >= inner_size)
      {
        l_aabb = leaf_aabb_ptr[lchild - inner_size];
        lchild = -(lchild - inner_size + 1);
      }
      else
      {
        l_aabb = inner_aabb_ptr[lchild];
        // do the offset now
        lchild *= 2;
      }

      std::int32_t rchild = rchildren_ptr[node];
      if(rchild >= inner_size)
      {
        r_aabb = leaf_aabb_ptr[rchild - inner_size];
        rchild = -(rchild - inner_size + 1);
      }
      else
      {
        r_aabb = inner_aabb_ptr[rchild];
        // do the offset now
        rchild *= 2;
      }

      const std::int32_t out_offset = node * 2;
      bvh_inner_nodes[out_offset + 0] = l_aabb;
      bvh_inner_nodes[out_offset + 1] = r_aabb;

      bvh_inner_node_children[out_offset + 0] = lchild;
      bvh_inner_node_children[out_offset + 1] = rchild;
    });
  AXOM_ANNOTATE_END("emit_bvh_parents");

  m_leaf_nodes = std::move(radix_tree.m_leafs);

  m_initialized = true;
}

template <typename FloatType, int NDIMS, typename ExecSpace>
template <typename PrimitiveType, typename Predicate, typename PrimitiveIndexable>
axom::Array<IndexType> LinearBVH<FloatType, NDIMS, ExecSpace>::findCandidatesImpl(
  Predicate&& predicate,
  const axom::ArrayView<IndexType> offsets,
  const axom::ArrayView<IndexType> counts,
  IndexType numObjs,
  PrimitiveIndexable objs,
  int allocatorID) const

{
  AXOM_ANNOTATE_SCOPE("LinearBVH::findCandidatesImpl");

  SLIC_ERROR_IF(offsets.size() != numObjs, "offsets length not equal to numObjs");
  SLIC_ERROR_IF(counts.size() != numObjs, "counts length not equal to numObjs");
  SLIC_ASSERT(m_initialized);

  const auto inner_nodes = m_inner_nodes.view();
  const auto inner_node_children = m_inner_node_children.view();
  const auto leaf_nodes = m_leaf_nodes.view();

  auto noTraversePref = [] AXOM_HOST_DEVICE(const BoundingBoxType&,
                                            const BoundingBoxType&,
                                            const PrimitiveType&) { return false; };

#if defined(AXOM_USE_RAJA)
  // STEP 1: count number of candidates for each query point
  axom::ReduceSum<ExecSpace, IndexType> total_count_reduce(0);

  AXOM_ANNOTATE_BEGIN("PASS[1]:count_traversal");
  for_all<ExecSpace>(
    numObjs,
    AXOM_LAMBDA(IndexType i) {
      std::int32_t count = 0;
      PrimitiveType primitive {objs[i]};

      auto leafAction = [&count](std::int32_t AXOM_UNUSED_PARAM(current_node),
                                 const std::int32_t* AXOM_UNUSED_PARAM(leaf_nodes)) { count++; };

      lbvh::bvh_traverse(inner_nodes,
                         inner_node_children,
                         leaf_nodes,
                         primitive,
                         predicate,
                         leafAction,
                         noTraversePref);

      counts[i] = count;
      total_count_reduce += count;
    });
  AXOM_ANNOTATE_END("PASS[1]:count_traversal");

  // STEP 2: exclusive scan to get offsets in candidate array for each query
  AXOM_ANNOTATE_BEGIN("exclusive_scan");
  axom::exclusive_scan<ExecSpace>(axom::ArrayView<IndexType>(counts.data(), numObjs),
                                  axom::ArrayView<IndexType>(offsets.data(), numObjs));

  AXOM_ANNOTATE_END("exclusive_scan");
  IndexType total_candidates = total_count_reduce.get();

  // STEP 3: allocate memory for all candidates
  AXOM_ANNOTATE_BEGIN("allocate_candidates");
  auto candidates = axom::Array<IndexType>(total_candidates, total_candidates, allocatorID);
  AXOM_ANNOTATE_END("allocate_candidates");
  const auto candidates_v = candidates.view();

  // STEP 4: fill in candidates for each point
  AXOM_ANNOTATE_BEGIN("PASS[2]:fill_traversal");
  for_all<ExecSpace>(
    numObjs,
    AXOM_LAMBDA(IndexType i) {
      std::int32_t offset = offsets[i];

      PrimitiveType obj {objs[i]};
      auto leafAction = [&offset, candidates_v](std::int32_t current_node, const std::int32_t* leafs) {
        candidates_v[offset] = leafs[current_node];
        offset++;
      };

      lbvh::bvh_traverse(inner_nodes,
                         inner_node_children,
                         leaf_nodes,
                         obj,
                         predicate,
                         leafAction,
                         noTraversePref);
    });
  AXOM_ANNOTATE_END("PASS[2]:fill_traversal");

  return candidates;
#else  // CPU-only and no RAJA: do single traversal
  AXOM_UNUSED_VAR(allocatorID);

  axom::Array<IndexType> search_candidates;
  int current_offset = 0;

  // STEP 1: do single-pass traversal with std::vector for candidates
  AXOM_ANNOTATE_BEGIN("PASS[1]:fill_traversal");
  for_all<ExecSpace>(numObjs, [&](IndexType i) {
    int matching_leaves = 0;
    PrimitiveType obj {objs[i]};
    offsets[i] = current_offset;

    auto leafAction = [&](std::int32_t current_node, const std::int32_t* leafs) {
      search_candidates.emplace_back(leafs[current_node]);
      matching_leaves++;
      current_offset++;
    };

    lbvh::bvh_traverse(inner_nodes,
                       inner_node_children,
                       leaf_nodes,
                       obj,
                       predicate,
                       leafAction,
                       noTraversePref);
    counts[i] = matching_leaves;
  });
  AXOM_ANNOTATE_END("PASS[1]:fill_traversal");

  SLIC_ASSERT(current_offset == static_cast<IndexType>(search_candidates.size()));

  return search_candidates;
#endif
}

template <typename FloatType, int NDIMS, typename ExecSpace>
void LinearBVH<FloatType, NDIMS, ExecSpace>::writeVtkFileImpl(const std::string& fileName) const
{
  std::ostringstream nodes;
  std::ostringstream cells;
  std::ostringstream levels;

  // STEP 0: Write VTK header
  std::ofstream ofs;
  ofs.open(fileName.c_str());
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " BVHTree \n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  // STEP 1: write root
  std::int32_t numPoints = 0;
  std::int32_t numBins = 0;
  lbvh::write_root(m_bounds, numPoints, numBins, nodes, cells, levels);

  // STEP 2: traverse the BVH and dump each bin
  constexpr std::int32_t ROOT = 0;
  lbvh::write_recursive<FloatType, NDIMS>(m_inner_nodes,
                                          m_inner_node_children,
                                          ROOT,
                                          1,
                                          numPoints,
                                          numBins,
                                          nodes,
                                          cells,
                                          levels);

  // STEP 3: write nodes
  ofs << "POINTS " << numPoints << " double\n";
  ofs << nodes.str() << std::endl;

  // STEP 4: write cells
  const std::int32_t nnodes = (NDIMS == 2) ? 4 : 8;
  ofs << "CELLS " << numBins << " " << numBins * (nnodes + 1) << std::endl;
  ofs << cells.str() << std::endl;

  // STEP 5: write cell types
  ofs << "CELL_TYPES " << numBins << std::endl;
  const std::int32_t cellType = (NDIMS == 2) ? 9 : 12;
  for(std::int32_t i = 0; i < numBins; ++i)
  {
    ofs << cellType << std::endl;
  }

  // STEP 6: dump level information
  ofs << "CELL_DATA " << numBins << std::endl;
  ofs << "SCALARS level int\n";
  ofs << "LOOKUP_TABLE default\n";
  ofs << levels.str() << std::endl;
  ofs << std::endl;

  // STEP 7: close file
  ofs.close();
}

}  // namespace policy
}  // namespace spin
}  // namespace axom
#endif  // AXOM_SPIN_POLICY_LINEARBVH_HPP_
