// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

// axom core includes
#include "axom/core/Types.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/AnnotationMacros.hpp"
#include "axom/core/numerics/floating_point_limits.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Vector.hpp"

// linear bvh includes
#include "axom/spin/SpaceFillingCurve.hpp"
#include "axom/spin/BucketedBVHGrouping.hpp"
#include "axom/spin/internal/linear_bvh/BVHNode.hpp"
#include "axom/spin/internal/linear_bvh/RadixTree.hpp"
#include "axom/spin/internal/linear_bvh/build_radix_tree.hpp"
#include "axom/spin/internal/linear_bvh/bvh_traverse.hpp"
#include "axom/spin/internal/linear_bvh/bvh_vtkio.hpp"

// C/C++ includes
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace axom::spin::policy
{
namespace lbvh = internal::linear_bvh;

template <typename FloatType, int NDIMS>
using BVH2Node = lbvh::BVH2Node<FloatType, NDIMS>;

/*!
 * \brief Tags the leaf representations supported by LinearBVH.
 */
struct LinearBVHType
{
  /// A leaf stores the index of exactly one input primitive.
  struct SinglePrimitiveLeaf
  { };

  /// A leaf stores a Morton-local bucket of input primitive indices.
  struct BucketedPrimitiveLeaves
  { };
};

namespace detail
{

/*!
 * \brief Contains optional indirection data for bucketed BVH leaves.
 *
 * \tparam EnableBucketIndirection When true, store the mapping from BVH
 *         leaves to primitive buckets. When false, remain empty.
 */
template <bool EnableBucketIndirection>
struct BucketLeafData
{ };

template <>
struct BucketLeafData<true>
{
  // Maps each linear-BVH leaf slot to its bucket ID.
  axom::Array<std::int32_t> m_leaf_bucket_ids;

  // Maps bucket IDs to half-open ranges in m_bucket_primitive_ids.
  axom::Array<IndexType> m_bucket_offsets;

  // Stores the original primitive IDs for every bucket, contiguously.
  axom::Array<std::int32_t> m_bucket_primitive_ids;
};

}  // namespace detail

/*
 * \brief Interface for a BVH tree through a traversal operation (which 
 *  searches a tree based on a user-provided predicate) or a reduce
 *  operation (which calculates a value at each node by adding its children).
 *
 * \brief Traverse a tree to perform user-specified actions at the
 * leaves, while limiting the search to branches satisfying a
 * user-provided predicate.
 * 
 * To initiate traversals, use \a traverse_tree. It requires
 * -# an action functor to call at the leaves,
 * -# a predicate functor to decide whether to descend a branch and
 * -# some data to pass to the functors
 * 
 * \brief Reduce a tree by invoking a user-specified leaf action at 
 * each leaf node, then iterating up the tree so the value at each node
 * is a sum of the value for its children.
 *
 * To initiate reductions, use \a reduce_tree. It requires
 * -# an action functor to call at the leaves,
 *
 * All functors should only access memory that's available to the
 * execution space.  For example, GPU execution should access only
 * device and unified memory.
 *
 * \see internal::linear_bvh::bvh_traverse
 */
template <typename FloatType, int NDIMS, typename LeafType = LinearBVHType::SinglePrimitiveLeaf>
class LinearBVHTraverser
{
public:
  static_assert(std::is_same_v<LeafType, LinearBVHType::SinglePrimitiveLeaf>,
                "Only LinearBVHType::SinglePrimitiveLeaf is currently supported.");

  using BoxType = primal::BoundingBox<FloatType, NDIMS>;
  using PointType = primal::Point<FloatType, NDIMS>;
  using BVHNode = BVH2Node<FloatType, NDIMS>;

  LinearBVHTraverser(axom::ArrayView<const BVHNode> nodes,
                     axom::ArrayView<const std::int32_t> leaf_nodes)
    : m_inner_nodes(nodes)
    , m_leaf_nodes(leaf_nodes)
  { }

  /*
   * Functors \a leaf_action and \a predicate should access only memory compatible
   * with the execution space.  For example, GPU execution should access
   * only device and unified memory.
   */
  template <typename Primitive, typename LeafAction, typename Predicate>
  AXOM_HOST_DEVICE void traverse_tree(const Primitive& p,
                                      LeafAction&& leaf_action,
                                      Predicate&& predicate) const
  {
    auto traversePref = [](const BoxType& l, const BoxType& r, const Primitive& p) {
      return LinearBVHTraverser::traverseClosestFirst(l, r, p);
    };

    lbvh::BVHStack stack;

    lbvh::bvh_traverse(m_inner_nodes, m_leaf_nodes, p, stack, predicate, leaf_action, traversePref);
  }

  template <typename ExecSpace, typename Primitive, typename LeafAction, typename Predicate>
  AXOM_HOST_DEVICE void traverseTreeShared(const Primitive& p,
                                           LeafAction&& leaf_action,
                                           Predicate&& predicate) const
  {
    auto noTraversePref = [](const BoxType& l, const BoxType& r, const Primitive& p) {
      return LinearBVHTraverser::traverseClosestFirst(l, r, p);
    };

    constexpr int BlockSize = axom::execution_space<ExecSpace>::BlockSize;

#ifdef AXOM_DEVICE_CODE
    lbvh::SharedBVHStack<BlockSize> stack;
#else
    AXOM_UNUSED_VAR(BlockSize);
    lbvh::BVHStack stack;
#endif

    lbvh::bvh_traverse(m_inner_nodes, m_leaf_nodes, p, stack, predicate, leaf_action, noTraversePref);
  }

  /*!
   * \brief Iterate over the tree, invoking the leaf action at each leaf node to
   *        produce a value and then iterate back up the tree, combining nodes
   *        using a "+" reduction. Return the Array that contains values for
   *        all tree nodes.
   *
   * \param leaf_action The function to invoke on a leaf node to make its data.
   * \param allocatorID The allocator to use to allocate array data.
   *
   * \return An Array that contains the reduced data for all nodes in the BVH.
   */
  template <typename ExecSpace, typename ValueType, typename LeafAction>
  axom::Array<ValueType> reduce_tree(LeafAction&& leaf_action,
                                     int allocatorID = axom::getDefaultAllocatorID()) const
  {
    // Make a field over all of the nodes (the return field).
    const auto num_node_slots = 2 * m_inner_nodes.size();
    axom::Array<ValueType> reducedField(num_node_slots, num_node_slots, allocatorID);

    if constexpr(std::is_same_v<ExecSpace, axom::SEQ_EXEC>)
    {
      reduce_recursion(std::forward<LeafAction>(leaf_action), reducedField.view(), 0);
      reduce_recursion(std::forward<LeafAction>(leaf_action), reducedField.view(), 1);
    }
    else
    {
      // Do it in 2 stages.

      // Make a field for just the leaf data. Compute it in parallel.
      axom::Array<ValueType> leafField(m_leaf_nodes.size(), m_leaf_nodes.size(), allocatorID);
      auto leafFieldView = leafField.view();
      const std::int32_t* leaf_nodes_data = m_leaf_nodes.data();
      axom::for_all<ExecSpace>(m_leaf_nodes.size(), [&](axom::IndexType currentNode) {
        const auto idx = leaf_nodes_data[currentNode];
        leafFieldView[idx] = leaf_action(static_cast<std::int32_t>(currentNode), leaf_nodes_data);
      });

      // Return the precomputed values in the reduction.
      auto returnLeafValue = [&](std::int32_t currentNode, const std::int32_t* leafNodes) {
        const auto idx = leafNodes[currentNode];
        return leafFieldView[idx];
      };

      // TODO: Replace this with GPU-compatible code.
      reduce_recursion(returnLeafValue, reducedField.view(), 0);
      reduce_recursion(returnLeafValue, reducedField.view(), 1);
    }

    return reducedField;
  }

private:
  template <typename PrimitiveType>
  AXOM_HOST_DEVICE static bool traverseClosestFirst(const BoxType& l,
                                                    const BoxType& r,
                                                    const PrimitiveType& p)
  {
    if constexpr(std::is_same_v<PrimitiveType, PointType>)
    {
      double sqDistL = primal::squared_distance(p, l.getCentroid());
      // If the right bbox is not valid, return max. Otherwise, the invalid right
      // bbox might actually win when we should ignore it.
      double sqDistR = r.isValid() ? primal::squared_distance(p, r.getCentroid())
                                   : axom::numerics::floating_point_limits<double>::max();
      return sqDistL > sqDistR;
    }
    else
    {
      return false;
    }
  }

  /*!
   * \brief This is a helper method used in reduce_tree.
   *
   * \param leaf_action The function to invoke on a leaf node to make its data.
   * \param node_data The view that contains the traversal order for leaf nodes.
   * \param current_node The current node.
   */
  template <typename ValueType, typename LeafAction>
  void reduce_recursion(LeafAction&& leaf_action,
                        axom::ArrayView<ValueType> node_data,
                        std::int32_t current_node) const
  {
    const auto& node = m_inner_nodes[current_node / 2];
    auto child_index = current_node % 2 == 0 ? node.left_child : node.right_child;

    if(child_index >= 0)
    {
      child_index *= 2;
    }

    // Check if node is a leaf
    if(child_index < 0)
    {
      node_data[current_node] = leaf_action(-child_index - 1, m_leaf_nodes.data());

      return;
    }

    // Populate children
    reduce_recursion(std::forward<LeafAction>(leaf_action), node_data, child_index + 0);
    reduce_recursion(std::forward<LeafAction>(leaf_action), node_data, child_index + 1);

    // Sum to get value for current node
    node_data[current_node] = node_data[child_index + 0] + node_data[child_index + 1];
  }

  axom::ArrayView<const BVHNode> m_inner_nodes;      // BVH bins including leafs
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
template <typename FloatType, int NDIMS, typename ExecSpace, typename LeafType = LinearBVHType::SinglePrimitiveLeaf>
class LinearBVH
{
public:
  using TraverserType = LinearBVHTraverser<FloatType, NDIMS, LeafType>;
  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;
  using BVHNode = BVH2Node<FloatType, NDIMS>;

  LinearBVH() = default;

  /// Returns whether this BVH has been built.
  bool isInitialized() const { return m_initialized; }

  /*!
   * \brief Builds a linear BVH with the given bounding boxes as leaf nodes.
   *
   * \param [in] boxes the bounding boxes for each leaf node
   * \param [in] numBoxes the number of bounding boxes
   * \param [in] scaleFactor scale factor applied to each bounding box before insertion into the BVH
   * \param [in] allocatorID allocator used to store the BVH
   * \param [in] ordering curve used to order box centroids. Defaults to Morton.
   */
  template <typename BoxIndexable>
  void buildImpl(const BoxIndexable boxes,
                 IndexType numBoxes,
                 FloatType scaleFactor,
                 int allocatorID,
                 SpaceFillingCurve ordering = SpaceFillingCurve::Morton);

  /*!
   * \brief Builds a bucketed linear BVH from spatially sorted primitives.
   *
   * \param [in] boxes the bounding boxes for each leaf node
   * \param [in] numBoxes the number of bounding boxes
   * \param [in] maxBucketSize maximum number of primitives in a bucket. The
   *               grouping policy enforces this upper bound.
   * \param [in] scaleFactor scale factor applied to each bounding box before insertion into the BVH
   * \param [in] allocatorID allocator for BVH and bucket storage
   * \param [in] ordering curve used to order box centroids.
   * \param [in] grouping policy used to form buckets from spatially sorted primitives.
   * \pre This LinearBVH requires BucketedPrimitiveLeaves.
   */
  template <typename BoxIndexable>
  void buildImpl(const BoxIndexable boxes,
                 IndexType numBoxes,
                 IndexType maxBucketSize,
                 FloatType scaleFactor,
                 int allocatorID,
                 SpaceFillingCurve ordering = SpaceFillingCurve::Morton,
                 BucketedBVHGrouping grouping = BucketedBVHGrouping::SharedPrefix);

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
    if constexpr(std::is_same_v<LeafType, LinearBVHType::SinglePrimitiveLeaf>)
    {
      return TraverserType(m_inner_nodes.view(), m_leaf_nodes.view());
    }
    else
    {
      return TraverserType(m_inner_nodes.view(),
                           m_leaf_storage.m_leaf_bucket_ids.view(),
                           m_leaf_storage.m_bucket_offsets.view(),
                           m_leaf_storage.m_bucket_primitive_ids.view());
    }
  }

private:
  /*!
   * \brief Emit LinearBVH node and leaf storage from a completed radix tree.
   *
   * \param[in,out] radix_tree Intermediate tree whose leaf-ID array is moved
   *            into this BVH.
   * \param[in] bounds Bounds of the original primitive boxes.
   * \param[in] allocatorID Allocator for the emitted internal-node storage.
   */
  void emitRadixTree(lbvh::RadixTree<FloatType, NDIMS>& radix_tree,
                     const BoundingBoxType& bounds,
                     int allocatorID);

  void allocate(std::int32_t size, int allocID)
  {
    AXOM_ANNOTATE_SCOPE("LinearBVH::allocate");
    IndexType numInnerNodes = size - 1;
    // Need to allocate this uninitialized, since primal::BoundingBox is
    // considered non-trivially-copyable on GCC 4.9.3
    m_inner_nodes =
      axom::Array<BVHNode>(axom::ArrayOptions::Uninitialized {}, numInnerNodes, numInnerNodes, allocID);
    m_leaf_nodes = axom::Array<std::int32_t>(size, size, allocID);
  }

  bool m_initialized {false};
  axom::Array<BVHNode> m_inner_nodes;      // BVH bins including leafs
  axom::Array<std::int32_t> m_leaf_nodes;  // leaf data
  primal::BoundingBox<FloatType, NDIMS> m_bounds;

  // Optional storage for bucketed primitives,
  // Resolves to an empty type unless this LinearBVH uses bucketed primitives.
  [[no_unique_address]] detail::BucketLeafData<
    std::is_same_v<LeafType, LinearBVHType::BucketedPrimitiveLeaves>> m_leaf_storage;
};

template <typename FloatType, int NDIMS, typename ExecSpace, typename LeafType>
template <typename BoxIndexable>
void LinearBVH<FloatType, NDIMS, ExecSpace, LeafType>::buildImpl(const BoxIndexable boxes,
                                                                 IndexType numBoxes,
                                                                 FloatType scaleFactor,
                                                                 int allocatorID,
                                                                 SpaceFillingCurve ordering)
{
  static_assert(std::is_same_v<LeafType, LinearBVHType::SinglePrimitiveLeaf>,
                "This build overload requires LinearBVHType::SinglePrimitiveLeaf.");

  AXOM_ANNOTATE_SCOPE("LinearBVH::buildImpl");

  // STEP 1: Sort the boxes by spatial key and build a RadixTree.
  SLIC_ASSERT(numBoxes <= std::numeric_limits<std::int32_t>::max());
  const auto numBoxesInt = static_cast<std::int32_t>(numBoxes);

  lbvh::RadixTree<FloatType, NDIMS> radix_tree;
  primal::BoundingBox<FloatType, NDIMS> global_bounds;
  lbvh::build_radix_tree<ExecSpace>(boxes,
                                    numBoxesInt,
                                    global_bounds,
                                    radix_tree,
                                    scaleFactor,
                                    allocatorID,
                                    ordering);

  // STEP 2: emit the BVH data-structure from the radix tree via LinearBVH::allocate
  // STEP 3: emit the BVH
  emitRadixTree(radix_tree, global_bounds, allocatorID);
}

template <typename FloatType, int NDIMS, typename ExecSpace, typename LeafType>
template <typename BoxIndexable>
void LinearBVH<FloatType, NDIMS, ExecSpace, LeafType>::buildImpl(const BoxIndexable boxes,
                                                                 IndexType numBoxes,
                                                                 IndexType maxBucketSize,
                                                                 FloatType scaleFactor,
                                                                 int allocatorID,
                                                                 SpaceFillingCurve ordering,
                                                                 BucketedBVHGrouping grouping)
{
  static_assert(std::is_same_v<LeafType, LinearBVHType::BucketedPrimitiveLeaves>,
                "This build overload requires LinearBVHType::BucketedPrimitiveLeaves.");

  AXOM_ANNOTATE_SCOPE("LinearBVH::buildBucketedImpl");
  SLIC_ASSERT(numBoxes <= std::numeric_limits<std::int32_t>::max());
  SLIC_ASSERT(maxBucketSize > 0);
  SLIC_ASSERT(maxBucketSize <= std::numeric_limits<std::int32_t>::max());
  const auto numBoxesInt = static_cast<std::int32_t>(numBoxes);
  const auto maxBucketSizeInt = static_cast<std::int32_t>(maxBucketSize);

  lbvh::RadixTree<FloatType, NDIMS> radix_tree;
  primal::BoundingBox<FloatType, NDIMS> global_bounds;
  axom::Array<IndexType> bucket_offsets;
  axom::Array<std::int32_t> bucket_primitive_ids;
  lbvh::build_bucketed_radix_tree<ExecSpace>(boxes,
                                             numBoxesInt,
                                             global_bounds,
                                             radix_tree,
                                             bucket_offsets,
                                             bucket_primitive_ids,
                                             maxBucketSizeInt,
                                             scaleFactor,
                                             allocatorID,
                                             ordering,
                                             grouping);

  emitRadixTree(radix_tree, global_bounds, allocatorID);
  m_leaf_storage.m_leaf_bucket_ids = std::move(m_leaf_nodes);
  m_leaf_storage.m_bucket_offsets = std::move(bucket_offsets);
  m_leaf_storage.m_bucket_primitive_ids = std::move(bucket_primitive_ids);
}

template <typename FloatType, int NDIMS, typename ExecSpace, typename LeafType>
void LinearBVH<FloatType, NDIMS, ExecSpace, LeafType>::emitRadixTree(
  lbvh::RadixTree<FloatType, NDIMS>& radix_tree,
  const BoundingBoxType& bounds,
  int allocatorID)
{
  // Emit the LinearBVH representation from a completed radix tree.
  m_bounds = bounds;
  allocate(radix_tree.m_size, allocatorID);

  const std::int32_t size = radix_tree.m_size;
  AXOM_UNUSED_VAR(size);
  const std::int32_t inner_size = radix_tree.m_inner_size;
  SLIC_ASSERT(inner_size == size - 1);

  const auto lchildren_ptr = radix_tree.m_left_children.view();
  const auto rchildren_ptr = radix_tree.m_right_children.view();

  const auto leaf_aabb_ptr = radix_tree.m_leaf_aabbs.view();
  const auto inner_aabb_ptr = radix_tree.m_inner_aabbs.view();

  const auto bvh_inner_nodes = m_inner_nodes.view();

  AXOM_ANNOTATE_BEGIN("emit_bvh_parents");
  for_all<ExecSpace>(inner_size, [=] AXOM_HOST_DEVICE(std::int32_t node) {
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
    }

    bvh_inner_nodes[node].left = l_aabb;
    bvh_inner_nodes[node].right = r_aabb;

    bvh_inner_nodes[node].left_child = lchild;
    bvh_inner_nodes[node].right_child = rchild;
  });
  AXOM_ANNOTATE_END("emit_bvh_parents");

  m_leaf_nodes = std::move(radix_tree.m_leafs);

  m_initialized = true;
}

template <typename FloatType, int NDIMS, typename ExecSpace, typename LeafType>
template <typename PrimitiveType, typename Predicate, typename PrimitiveIndexable>
axom::Array<IndexType> LinearBVH<FloatType, NDIMS, ExecSpace, LeafType>::findCandidatesImpl(
  Predicate&& predicate,
  const axom::ArrayView<IndexType> offsets,
  const axom::ArrayView<IndexType> counts,
  IndexType numObjs,
  PrimitiveIndexable objs,
  int allocatorID) const

{
  static_assert(std::is_same_v<LeafType, LinearBVHType::SinglePrimitiveLeaf>,
                "findCandidatesImpl() is not implemented for bucketed leaves.");

  AXOM_ANNOTATE_SCOPE("LinearBVH::findCandidatesImpl");

  SLIC_ERROR_IF(offsets.size() != numObjs, "offsets length not equal to numObjs");
  SLIC_ERROR_IF(counts.size() != numObjs, "counts length not equal to numObjs");
  SLIC_ASSERT(m_initialized);

  TraverserType tree_view = this->getTraverserImpl();

#if defined(AXOM_USE_RAJA)
  // STEP 1: count number of candidates for each query point
  axom::ReduceSum<ExecSpace, IndexType> total_count_reduce(0);

  AXOM_ANNOTATE_BEGIN("PASS[1]:count_traversal");
  for_all<ExecSpace>(numObjs, [=] AXOM_HOST_DEVICE(IndexType i) {
    std::int32_t count = 0;
    PrimitiveType primitive {objs[i]};

    auto leafAction = [&count](std::int32_t AXOM_UNUSED_PARAM(current_node),
                               const std::int32_t* AXOM_UNUSED_PARAM(leaf_nodes)) { count++; };

    tree_view.traverse_tree(primitive, leafAction, predicate);

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
  for_all<ExecSpace>(numObjs, [=] AXOM_HOST_DEVICE(IndexType i) {
    std::int32_t offset = offsets[i];

    PrimitiveType obj {objs[i]};
    auto leafAction = [&offset, candidates_v](std::int32_t current_node, const std::int32_t* leafs) {
      candidates_v[offset] = leafs[current_node];
      offset++;
    };

    tree_view.traverse_tree(obj, leafAction, predicate);
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
    IndexType matching_leaves = 0;
    PrimitiveType obj {objs[i]};
    offsets[i] = current_offset;

    auto leafAction = [&](std::int32_t current_node, const std::int32_t* leafs) {
      search_candidates.emplace_back(leafs[current_node]);
      matching_leaves++;
      current_offset++;
    };

    tree_view.traverse_tree(obj, leafAction, predicate);

    counts[i] = matching_leaves;
  });
  AXOM_ANNOTATE_END("PASS[1]:fill_traversal");

  SLIC_ASSERT(current_offset == static_cast<IndexType>(search_candidates.size()));

  return search_candidates;
#endif
}

template <typename FloatType, int NDIMS, typename ExecSpace, typename LeafType>
void LinearBVH<FloatType, NDIMS, ExecSpace, LeafType>::writeVtkFileImpl(const std::string& fileName) const
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
  lbvh::write_recursive<FloatType, NDIMS>(m_inner_nodes, ROOT, 1, numPoints, numBins, nodes, cells, levels);

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

///@{
//! @name Specializations for LinearBVHTraverser<FloatType, NDIMS, LinearBVHType::BucketedPrimitiveLeaves>

/*!
 * \brief A device-copyable view of the original primitive IDs in one BVH leaf bucket.
 */
class LinearBVHLeafRange
{
public:
  AXOM_HOST_DEVICE LinearBVHLeafRange(const std::int32_t* primitive_ids, IndexType size)
    : m_primitive_ids(primitive_ids)
    , m_size(size)
  { }

  AXOM_HOST_DEVICE IndexType size() const { return m_size; }

  AXOM_HOST_DEVICE std::int32_t operator[](IndexType index) const { return m_primitive_ids[index]; }

private:
  const std::int32_t* m_primitive_ids;
  IndexType m_size;
};

/*!
 * \brief Traverses or reduces a linear BVH whose leaves contain primitive buckets.
 *
 * This specialization adapts the single-primitive traverser so each leaf action
 * receives a LinearBVHLeafRange for the primitives in that leaf's bucket.
 * Bucketed leaves do not support shared-memory traversal.
 */
template <typename FloatType, int NDIMS>
class LinearBVHTraverser<FloatType, NDIMS, LinearBVHType::BucketedPrimitiveLeaves>
{
public:
  using BVHNode = BVH2Node<FloatType, NDIMS>;
  using SinglePrimitiveTraverser = LinearBVHTraverser<FloatType, NDIMS>;

  LinearBVHTraverser(axom::ArrayView<const BVHNode> nodes,
                     axom::ArrayView<const std::int32_t> leaf_bucket_ids,
                     axom::ArrayView<const IndexType> bucket_offsets,
                     axom::ArrayView<const std::int32_t> bucket_primitive_ids)
    : m_traverser(nodes, leaf_bucket_ids)
    , m_bucket_offsets(bucket_offsets)
    , m_bucket_primitive_ids(bucket_primitive_ids)
  { }

  /*!
   * \brief Traverse the tree and invoke an action for each intersected leaf bucket.
   *
   * \param[in] p The primitive to use for traversal.
   * \param[in] leaf_action Action invoked with a LinearBVHLeafRange for each
   *            visited leaf bucket.
   * \param[in] predicate Predicate that determines whether to descend a branch.
   *
   * \note Functors must access only memory compatible with the execution space.
   */
  template <typename Primitive, typename LeafAction, typename Predicate>
  AXOM_HOST_DEVICE void traverse_tree(const Primitive& p,
                                      LeafAction&& leaf_action,
                                      Predicate&& predicate) const
  {
    auto action = std::forward<LeafAction>(leaf_action);
    const auto bucket_offsets = m_bucket_offsets;
    const auto bucket_primitive_ids = m_bucket_primitive_ids;
    auto bucket_action = [=] AXOM_HOST_DEVICE(std::int32_t leaf_slot,
                                              const std::int32_t* leaf_bucket_ids) mutable {
      const auto bucket_id = leaf_bucket_ids[leaf_slot];
      const IndexType begin = bucket_offsets[bucket_id];
      action(LinearBVHLeafRange(bucket_primitive_ids.data() + begin,
                                bucket_offsets[bucket_id + 1] - begin));
    };

    m_traverser.traverse_tree(p, bucket_action, std::forward<Predicate>(predicate));
  }

  template <typename ExecSpace, typename Primitive, typename LeafAction, typename Predicate>
  AXOM_HOST_DEVICE void traverseTreeShared(const Primitive&, LeafAction&&, Predicate&&) const
  {
    static_assert(!std::is_same_v<ExecSpace, ExecSpace>,
                  "traverseTreeShared() is not supported for bucketed leaves.");
  }

  /*!
   * \brief Invoke a leaf action for each primitive bucket, then sum the values
   *        upward through the tree.
   *
   * \param[in] leaf_action Action invoked with a LinearBVHLeafRange for each
   *            leaf bucket. It returns the value used in the reduction.
   * \param[in] allocatorID Allocator used for the returned array.
   *
   * \return An array of reduced values for all BVH nodes.
   */
  template <typename ExecSpace, typename ValueType, typename LeafAction>
  axom::Array<ValueType> reduce_tree(LeafAction&& leaf_action,
                                     int allocatorID = axom::getDefaultAllocatorID()) const
  {
    auto action = std::forward<LeafAction>(leaf_action);
    const auto bucket_offsets = m_bucket_offsets;
    const auto bucket_primitive_ids = m_bucket_primitive_ids;
    auto bucket_action = [=](std::int32_t leaf_slot,
                             const std::int32_t* leaf_bucket_ids) mutable -> ValueType {
      const auto bucket_id = leaf_bucket_ids[leaf_slot];
      const IndexType begin = bucket_offsets[bucket_id];
      return action(LinearBVHLeafRange(bucket_primitive_ids.data() + begin,
                                       bucket_offsets[bucket_id + 1] - begin));
    };

    return m_traverser.template reduce_tree<ExecSpace, ValueType>(bucket_action, allocatorID);
  }

private:
  SinglePrimitiveTraverser m_traverser;
  axom::ArrayView<const IndexType> m_bucket_offsets;
  axom::ArrayView<const std::int32_t> m_bucket_primitive_ids;
};
///@}

}  // namespace axom::spin::policy
