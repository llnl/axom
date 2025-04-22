// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_TRAVERSE_HPP_
#define AXOM_SPIN_BVH_TRAVERSE_HPP_

#include "axom/config.hpp"       // compile-time definitions
#include "axom/core/Macros.hpp"  // for AXOM_HOST_DEVICE
#include "axom/core/Types.hpp"   // for axom types
#include "axom/slic.hpp"         // for SLIC macros

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
/*!
 * \brief Checks if the node corresponding to the given node ID is a leaf node.
 * \param [in] nodeIdx index of the BVH node in query.
 * \return status true if the node is a leaf node, otherwise, false.
 */
AXOM_HOST_DEVICE
inline bool leaf_node(const std::int32_t& nodeIdx) { return (nodeIdx < 0); }

/*!
 * \brief Generic BVH traversal routine.
 *
 * \param [in] inner_nodes pointer to the BVH bins.
 * \param [in] inner_node_children pointer to pairs of child indices.
 * \param [in] leaf_nodes pointer to the leaf node IDs.
 * \param [in] p the primitive in query, e.g., a point, ray, etc.
 * \param [in] B functor that defines the check for the bins
 * \param [in] A functor that defines the leaf action
 * \param [in] Comp functor used for determining which child node to traverse
 *  down first if both bins are to be traversed
 *
 * \note The supplied functor `B` is expected to take the following two
 *  arguments:
 *    (1) The supplied primitive, p
 *    (2) a primal::BoundingBox< FloatType, NDIMS > of the BVH bin
 *
 * \note The supplied functor `Comp` is expected to take the following three
 *  arguments:
 *    (1) The left child bounding box
 *    (2) The right child bounding box
 *    (3) The primitive being queried
 *  It should return true if the primitive is closer to the right child bounding
 *  box (indicating a swap is necessary) and false if the primitive is closer to
 *  the left child bounding box.
 *
 * \see BVHData for the details on the internal data layout of the BVH.
 *
 * \note Moreover, the functor `B` returns a boolean status that indicates
 *  if the specified traversal predicate is satisfied.
 *
 * \note Functors A, B and Comp may access only memory available in
 * the execution space.  For example, GPU execution may access only
 * device and unified memory.
 *
 */
template <int NDIMS, typename FloatType, typename PrimitiveType, typename InBinCheck, typename LeafAction, typename TraversePref>
AXOM_HOST_DEVICE inline void bvh_traverse(
  axom::ArrayView<const primal::BoundingBox<FloatType, NDIMS>> inner_nodes,
  axom::ArrayView<const std::int32_t> inner_node_children,
  axom::ArrayView<const std::int32_t> leaf_nodes,
  const PrimitiveType& p,
  InBinCheck&& B,
  LeafAction&& A,
  TraversePref&& Comp)
{
  using BBoxType = primal::BoundingBox<FloatType, NDIMS>;

  // setup stack
  constexpr std::int32_t STACK_SIZE = 64;
  constexpr std::int32_t BARRIER = -2000000000;
  std::int32_t todo[STACK_SIZE];
  std::int32_t stackptr = 0;
  todo[stackptr] = BARRIER;

  std::int32_t found_leaf = 0;
  std::int32_t current_node = 0;

  while(current_node != BARRIER)
  {
    // Traverse until we hit a leaf node or the barrier.
    while(!leaf_node(current_node))
    {
      BBoxType left_bin = inner_nodes[current_node + 0];
      BBoxType right_bin = inner_nodes[current_node + 1];
      const bool in_left = left_bin.isValid() ? B(p, left_bin) : false;
      const bool in_right = right_bin.isValid() ? B(p, right_bin) : false;
      std::int32_t l_child = inner_node_children[current_node + 0];
      std::int32_t r_child = inner_node_children[current_node + 1];
      bool swap = Comp(left_bin, right_bin, p);

      if(!in_left && !in_right)
      {
        // pop the stack and continue
        current_node = todo[stackptr];
        stackptr--;
      }
      else
      {
        current_node = (in_left) ? l_child : r_child;
        if(in_left && in_right)
        {
          if(swap)
          {
            // Ensure we go down the closer of the two children.
            // (For a user-defined meaning of "closer")
            axom::utilities::swap(current_node, r_child);
          }

          stackptr++;
          todo[stackptr] = r_child;
        }
      }  // END else

      if(leaf_node(current_node) && !leaf_node(found_leaf))
      {
        // Save this leaf and continue traversing
        found_leaf = current_node;
        if(current_node != BARRIER)
        {
          current_node = todo[stackptr];
          stackptr--;
        }
      }
    }  // END while

    // After the traversal, each thread may have found:
    // - two leaf nodes (found_leaf=l1, current_node=l2)
    // - one leaf node (found_leaf=l1, current_node=BARRIER)
    // - no leaf nodes (found_leaf=0, current_node=BARRIER)
    while(leaf_node(found_leaf) && found_leaf != BARRIER)
    {
      int leaf_idx = -found_leaf - 1;
      A(leaf_idx, leaf_nodes.data());
      found_leaf = current_node;
      if(leaf_node(current_node) && current_node != BARRIER)
      {
        // pop the stack and continue
        current_node = todo[stackptr];
        stackptr--;
      }
    }
    found_leaf = 0;
  }  // END while
}

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_SPIN_BVH_TRAVERSE_HPP_ */
