// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"       // compile-time definitions
#include "axom/core/Macros.hpp"  // for AXOM_HOST_DEVICE
#include "axom/core/Types.hpp"   // for axom types
#include "axom/slic.hpp"         // for SLIC macros
#include "axom/spin/internal/linear_bvh/BVHNode.hpp"

#include <type_traits>  // For template magic
#include <utility>

namespace axom
{
namespace spin
{
namespace internal
{
/*!
 * \brief Checks if provided InBinCheck can accept a third positional argument representing the
 *        index of `bvhBox` in the traverser's internal node list.
 * 
 * This is necessary because typical usage of a BVH traverser requires no access to the internal
 *  node layout as implemented in LinearBVH, requiring a two-argument InBin predicate.
 * However, `LinearBVHTraverser::reduce_tree` returns an array of values that are associated with 
 *  each node, and a three-argument InBin predicate allows indexing into this array in the check.
 * 
 * \sa quest::fast_approximate_winding_number
 */
template <typename InBinCheck, typename PointType, typename BoxType, typename IndexType>
AXOM_HOST_DEVICE inline bool invoke_InBinCheck(InBinCheck&& B,
                                               PointType&& p,
                                               BoxType&& bvhBox,
                                               IndexType&& node_idx)
{
  // Case 1: InBinCheck is invokable as B( p, bvhBox, node_idx )
  if constexpr(std::is_invocable_v<InBinCheck, PointType, BoxType, IndexType>)
  {
    return std::forward<InBinCheck>(B)(std::forward<PointType>(p),
                                       std::forward<BoxType>(bvhBox),
                                       std::forward<IndexType>(node_idx));
  }
  // Case 2: InBinCheck is invokable as B( p, bvhBox ), which is more common
  else if constexpr(std::is_invocable_v<InBinCheck, PointType, BoxType>)
  {
    return std::forward<InBinCheck>(B)(std::forward<PointType>(p), std::forward<BoxType>(bvhBox));
  }
  // Case 3: Neither works
  else
  {
    static_assert(std::is_invocable_v<InBinCheck, PointType, BoxType, IndexType> ||
                    std::is_invocable_v<InBinCheck, PointType, BoxType>,
                  "InBinCheck must be callable as B(p, bvhBox, node_idx) or B(p, bvhBox)");
    return false;
  }
}

namespace linear_bvh
{
/*!
 * \brief Checks if the node corresponding to the given node ID is a leaf node.
 * \param [in] nodeIdx index of the BVH node in query.
 * \return status true if the node is a leaf node, otherwise, false.
 */
AXOM_HOST_DEVICE
inline bool leaf_node(const std::int32_t& nodeIdx) { return (nodeIdx < 0); }

struct BVHStack
{
public:
  constexpr static std::int32_t STACK_SIZE = 64;
  constexpr static std::int32_t BARRIER = -2000000000;

  using LocalStack = std::int32_t[STACK_SIZE];
  AXOM_HOST_DEVICE BVHStack() { }

  AXOM_HOST_DEVICE void setLocalStack(LocalStack& local_stack)
  {
    stack_ptr = 0;
    stack = &(local_stack[0]);
    stack[stack_ptr] = BARRIER;
  }

  AXOM_HOST_DEVICE std::int32_t pop()
  {
    std::int32_t top = stack[stack_ptr];
    stack_ptr--;
    return top;
  }

  AXOM_HOST_DEVICE void push(std::int32_t value)
  {
    stack_ptr++;
    stack[stack_ptr] = value;
  }

private:
  std::int32_t stack_ptr {-1};
  std::int32_t* stack {nullptr};
};

template <int BlockSize>
struct SharedBVHStack
{
public:
  constexpr static std::int32_t CHUNK_SIZE = 4;
  constexpr static std::int32_t STACK_SIZE = 16;
  constexpr static std::int32_t SHMEM_SIZE_PER_THREAD = CHUNK_SIZE * 2;
  constexpr static std::int32_t BARRIER = -2000000000;

  struct Chunk
  {
    std::int32_t values[CHUNK_SIZE];
  };

  using LocalStack = Chunk[STACK_SIZE];

  AXOM_DEVICE static int* Get_Shared_Mem_Buffer()
  {
    __shared__ int shmem_buf[SHMEM_SIZE_PER_THREAD * BlockSize];
    return &(shmem_buf[0]);
  }

  AXOM_DEVICE SharedBVHStack()
    : s_block_dim(blockDim.x)
    , s_thread_id(threadIdx.x)
    , s_stack(SharedBVHStack::Get_Shared_Mem_Buffer())
  { }

  AXOM_HOST_DEVICE void setLocalStack(LocalStack& local_stack) { g_stack = &(local_stack[0]); }

  AXOM_DEVICE std::int32_t pop()
  {
    // Shared is empty, try to refill a chunk.
    if(s_ptr == 0 && g_ptr > 0)
    {
      --g_ptr;
      Chunk g_top_chunk = g_stack[g_ptr];
      for(int i = 0; i < CHUNK_SIZE; i++)
      {
        shared_stack(s_ptr + i) = g_top_chunk.values[i];
      }
      s_ptr += CHUNK_SIZE;
    }
    if(s_ptr > 0)
    {
      // Can pop directly from shared.
      --s_ptr;
      std::int32_t top = shared_stack(s_ptr);
      return top;
    }
    else
    {
      // Empty stack, return barrier.
      return BARRIER;
    }
  }

  AXOM_DEVICE void push(std::int32_t value)
  {
    if(s_ptr == 2 * CHUNK_SIZE)
    {
      // At capacity. Take bottom values and push onto global memory stack.
      Chunk g_bottom_chunk;
      for(int i = 0; i < CHUNK_SIZE; i++)
      {
        g_bottom_chunk.values[i] = shared_stack(i);
      }
      g_stack[g_ptr] = g_bottom_chunk;
      g_ptr++;
      // Move remaining stack values down.
      for(int i = 0; i < CHUNK_SIZE; i++)
      {
        shared_stack(i) = shared_stack(i + CHUNK_SIZE);
      }
      s_ptr -= CHUNK_SIZE;
    }
    assert(s_ptr < 2 * CHUNK_SIZE);
    // Push value onto shared stack.
    shared_stack(s_ptr) = value;
    s_ptr++;
  }

private:
  AXOM_DEVICE std::int32_t& shared_stack(int index)
  {
    return s_stack[index * s_block_dim + s_thread_id];
  }

  std::int16_t s_block_dim;
  std::int16_t s_thread_id;
  std::int32_t s_ptr {0};
  std::int32_t* s_stack;
  std::int32_t g_ptr {0};
  Chunk* g_stack;
};

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
template <int NDIMS,
          typename FloatType,
          typename PrimitiveType,
          typename TraverseStack,
          typename InBinCheck,
          typename LeafAction,
          typename TraversePref>
AXOM_HOST_DEVICE inline void bvh_traverse(axom::ArrayView<const BVH2Node<FloatType, NDIMS>> inner_nodes,
                                          axom::ArrayView<const std::int32_t> leaf_nodes,
                                          const PrimitiveType& p,
                                          TraverseStack& stack,
                                          InBinCheck&& B,
                                          LeafAction&& A,
                                          TraversePref&& Comp)
{
  using BBoxType = primal::BoundingBox<FloatType, NDIMS>;
  using BVHNode = BVH2Node<FloatType, NDIMS>;

  // setup stack
  typename TraverseStack::LocalStack local_mem;
  stack.setLocalStack(local_mem);

  std::int32_t found_leaf = 0;
  std::int32_t current_node = 0;

  while(current_node != TraverseStack::BARRIER)
  {
    // Traverse until we hit a leaf node or the barrier.
    while(!leaf_node(current_node))
    {
      BVHNode curr_node;
      curr_node = inner_nodes[current_node];
      BBoxType left_bin = curr_node.left;
      BBoxType right_bin = curr_node.right;
      const bool in_left =
        left_bin.isValid() ? invoke_InBinCheck(B, p, left_bin, current_node + 0) : false;
      const bool in_right =
        right_bin.isValid() ? invoke_InBinCheck(B, p, right_bin, current_node + 1) : false;
      std::int32_t l_child = curr_node.left_child;
      std::int32_t r_child = curr_node.right_child;
      bool swap = Comp(left_bin, right_bin, p);

      if(!in_left && !in_right)
      {
        // pop the stack and continue
        current_node = stack.pop();
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

          stack.push(r_child);
        }
      }  // END else

      if(leaf_node(current_node) && !leaf_node(found_leaf))
      {
        // Save this leaf and continue traversing
        found_leaf = current_node;
        if(current_node != TraverseStack::BARRIER)
        {
          current_node = stack.pop();
        }
      }
    }  // END while

    // After the traversal, each thread may have found:
    // - two leaf nodes (found_leaf=l1, current_node=l2)
    // - one leaf node (found_leaf=l1, current_node=BARRIER)
    // - no leaf nodes (found_leaf=0, current_node=BARRIER)
    while(leaf_node(found_leaf) && found_leaf != TraverseStack::BARRIER)
    {
      int leaf_idx = -found_leaf - 1;
      A(leaf_idx, leaf_nodes.data());
      found_leaf = current_node;
      if(leaf_node(current_node) && current_node != TraverseStack::BARRIER)
      {
        // pop the stack and continue
        current_node = stack.pop();
      }
    }
    found_leaf = 0;
  }  // END while
}

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */
