// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BUILD_RADIX_TREE_H_
#define AXOM_SPIN_BUILD_RADIX_TREE_H_

#include "axom/config.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/atomics.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/reductions.hpp"
#include "axom/core/execution/sorts.hpp"
#include "axom/core/AnnotationMacros.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/BitUtilities.hpp"
#include "axom/core/NumericLimits.hpp"

#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"

#include "axom/spin/internal/linear_bvh/RadixTree.hpp"
#include "axom/spin/MortonIndex.hpp"

#include <atomic>

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA)
  // NOTE: uses the cub installation that is  bundled with RAJA
  #include "cub/device/device_radix_sort.cuh"
#endif

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
//------------------------------------------------------------------------------
//Returns 30 bit morton code for coordinates point is expected to be between [0,1]
template <typename FloatType, int Dims>
static inline AXOM_HOST_DEVICE std::uint32_t morton32_encode(const primal::Vector<FloatType, Dims>& point)
{
  using PointType = primal::Point<std::int32_t, Dims>;
  using axom::utilities::clampVal;

  //for a float, take the first 10 bits. Note, 2^10 = 1024
  constexpr int NUM_BITS_PER_DIM = 32 / Dims;
  constexpr FloatType FLOAT_TO_INT = 1 << NUM_BITS_PER_DIM;
  constexpr FloatType FLOAT_CEILING = FLOAT_TO_INT - 1;
  constexpr FloatType FLOAT_ZERO {0.};

  std::int32_t int_coords[Dims];
  for(int i = 0; i < Dims; i++)
  {
    int_coords[i] =
      static_cast<std::int32_t>(clampVal(point[i] * FLOAT_TO_INT, FLOAT_ZERO, FLOAT_CEILING));
  }

  return convertPointToMorton<std::uint32_t>(PointType(int_coords));
}

//------------------------------------------------------------------------------
//Returns 63 bit morton code for a 3D point with coordinates in the unit cube
static inline AXOM_HOST_DEVICE std::uint64_t morton64_encode(float x, float y, float z = 0.0)
{
  using PointType = primal::Point<std::int64_t, 3>;
  using axom::utilities::clampVal;

  //take the first 21 bits. Note, 2^21= 2097152.0f
  constexpr float F_21 = 2097152.0f;
  return convertPointToMorton<std::uint64_t>(
    PointType {static_cast<std::int64_t>(clampVal(x * F_21, 0.f, F_21)),
               static_cast<std::int64_t>(clampVal(y * F_21, 0.f, F_21)),
               static_cast<std::int64_t>(clampVal(z * F_21, 0.f, F_21))});
}

template <typename ExecSpace, typename BoxIndexable, typename FloatType, int NDIMS>
void transform_boxes(const BoxIndexable boxes,
                     ArrayView<primal::BoundingBox<FloatType, NDIMS>> aabbs,
                     std::int32_t size,
                     FloatType scale_factor)
{
  AXOM_ANNOTATE_SCOPE("transform_boxes");

  for_all<ExecSpace>(
    size,
    AXOM_LAMBDA(std::int32_t i) {
      aabbs[i] = boxes[i];
      aabbs[i].scale(scale_factor);
    });
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
primal::BoundingBox<FloatType, NDIMS> reduce(ArrayView<const primal::BoundingBox<FloatType, NDIMS>> aabbs,
                                             std::int32_t size)
{
  AXOM_ANNOTATE_SCOPE("reduce_abbs");

#ifdef AXOM_USE_RAJA
  primal::Point<FloatType, NDIMS> min_pt, max_pt;

  FloatType infinity = axom::numeric_limits<FloatType>::max();
  FloatType neg_infinity = axom::numeric_limits<FloatType>::lowest();

  for(int dim = 0; dim < NDIMS; dim++)
  {
    axom::ReduceMin<ExecSpace, FloatType> min_coord(infinity);
    axom::ReduceMax<ExecSpace, FloatType> max_coord(neg_infinity);

    for_all<ExecSpace>(
      size,
      AXOM_LAMBDA(std::int32_t i) {
        const primal::BoundingBox<FloatType, NDIMS>& aabb = aabbs[i];
        min_coord.min(aabb.getMin()[dim]);
        max_coord.max(aabb.getMax()[dim]);
      });

    min_pt[dim] = min_coord.get();
    max_pt[dim] = max_coord.get();
  }

  return primal::BoundingBox<FloatType, NDIMS>(min_pt, max_pt);
#else
  static_assert(std::is_same<ExecSpace, SEQ_EXEC>::value, "Only SEQ_EXEC supported without RAJA");

  primal::BoundingBox<FloatType, NDIMS> global_bounds;

  for_all<ExecSpace>(size, [&](std::int32_t i) { global_bounds.addBox(aabbs[i]); });

  return global_bounds;
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
void get_mcodes(ArrayView<const primal::BoundingBox<FloatType, NDIMS>> aabbs,
                std::int32_t size,
                const primal::BoundingBox<FloatType, NDIMS>& bounds,
                const ArrayView<std::uint32_t> mcodes)
{
  AXOM_ANNOTATE_SCOPE("get_mcodes");

  primal::Vector<FloatType, NDIMS> extent, inv_extent, min_coord;

  extent = primal::Vector<FloatType, NDIMS>(bounds.getMax());
  extent -= primal::Vector<FloatType, NDIMS>(bounds.getMin());
  min_coord = primal::Vector<FloatType, NDIMS>(bounds.getMin());

  for(int i = 0; i < NDIMS; ++i)
  {
    inv_extent[i] = utilities::isNearlyEqual<FloatType>(extent[i], .0f) ? 0.f : 1.f / extent[i];
  }

  for_all<ExecSpace>(
    size,
    AXOM_LAMBDA(std::int32_t i) {
      const primal::BoundingBox<FloatType, NDIMS>& aabb = aabbs[i];

      // get the center and normalize it
      primal::Vector<FloatType, NDIMS> centroid(aabb.getCentroid());
      centroid =
        primal::Vector<FloatType, NDIMS>((centroid - min_coord).array() * inv_extent.array());
      mcodes[i] = morton32_encode(centroid);
    });
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename IntType>
void array_counting(ArrayView<IntType> iterator,
                    const IntType& size,
                    const IntType& start,
                    const IntType& step)
{
  AXOM_ANNOTATE_SCOPE("array_counting");

  for_all<ExecSpace>(size, AXOM_LAMBDA(std::int32_t i) { iterator[i] = start + i * step; });
}

//------------------------------------------------------------------------------
//
// reorder and array based on a new set of indices.
// array   [a,b,c]
// indices [1,0,2]
// result  [b,a,c]
//
template <typename ExecSpace, typename T>
void reorder(const ArrayView<const std::int32_t> indices,
             Array<T>& array,
             std::int32_t size,
             int allocatorID)
{
  AXOM_ANNOTATE_SCOPE("reorder");

  Array<T> temp = Array<T>(ArrayOptions::Uninitialized {}, size, size, allocatorID);
  const auto array_v = array.view();
  const auto temp_v = temp.view();

  for_all<ExecSpace>(
    size,
    AXOM_LAMBDA(std::int32_t i) {
      std::int32_t in_idx = indices[i];
      temp_v[i] = array_v[in_idx];
    });

  array = std::move(temp);
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_RAJA) && \
  ((RAJA_VERSION_MAJOR > 0) || ((RAJA_VERSION_MAJOR == 0) && (RAJA_VERSION_MINOR >= 12)))

template <typename ExecSpace>
void sort_mcodes(ArrayView<std::uint32_t> mcodes, std::int32_t size, ArrayView<std::int32_t> iter)
{
  AXOM_ANNOTATE_SCOPE("sort_mcodes");

  array_counting<ExecSpace>(iter, size, 0, 1);

  {
    AXOM_ANNOTATE_SCOPE("raja_stable_sort");
    ArrayView<std::uint32_t> mcode_view(mcodes.data(), size);
    ArrayView<std::int32_t> iter_view(iter.data(), size);
    axom::stable_sort_pairs<ExecSpace>(mcode_view, iter_view);
  }
}

#else

// fall back to std::stable_sort
template <typename ExecSpace>
void sort_mcodes(Array<std::uint32_t>& mcodes, std::int32_t size, const ArrayView<std::int32_t> iter)
{
  AXOM_ANNOTATE_SCOPE("sort_mcodes");

  array_counting<ExecSpace>(iter, size, 0, 1);

  {
    AXOM_ANNOTATE_SCOPE("cpu_sort");

    std::stable_sort(iter.begin(), iter.begin() + size, [&](std::int32_t i1, std::int32_t i2) {
      return mcodes[i1] < mcodes[i2];
    });
  }

  const int allocID = axom::execution_space<ExecSpace>::allocatorID();
  reorder<ExecSpace>(iter, mcodes, size, allocID);
}

#endif /* RAJA version 0.12.0 and above */

//------------------------------------------------------------------------------
template <typename IntType, typename MCType>
AXOM_HOST_DEVICE IntType
delta(const IntType& a, const IntType& b, const IntType& inner_size, axom::ArrayView<MCType> mcodes)
{
  bool tie = false;
  bool out_of_range = (b < 0 || b > inner_size);
  //still make the call but with a valid adderss
  const std::int32_t bb = (out_of_range) ? 0 : b;
  const std::uint32_t acode = mcodes[a];
  const std::uint32_t bcode = mcodes[bb];
  //use xor to find where they differ
  std::uint32_t exor = acode ^ bcode;
  tie = (exor == 0);
  //break the tie, a and b must always differ
  exor = tie ? std::uint32_t(a) ^ std::uint32_t(bb) : exor;
  std::int32_t count = axom::utilities::countl_zero(exor);
  if(tie)
  {
    count += 32;
  }
  count = (out_of_range) ? -1 : count;
  return count;
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
void build_tree(RadixTree<FloatType, NDIMS>& data)
{
  AXOM_ANNOTATE_SCOPE("build_tree");

  // http://research.nvidia.com/sites/default/files/publications/karras2012hpg_paper.pdf

  // Pointers and vars are redeclared because I have a faint memory
  // of a huge amount of pain and suffering due to cuda
  // lambda captures of pointers inside a struct. Bad memories
  // of random segfaults ........ be warned
  const std::int32_t inner_size = data.m_inner_size;
  const auto lchildren_ptr = data.m_left_children.view();
  const auto rchildren_ptr = data.m_right_children.view();
  const auto parent_ptr = data.m_parents.view();
  const auto mcodes_ptr = data.m_mcodes.view();

  for_all<ExecSpace>(
    inner_size,
    AXOM_LAMBDA(std::int32_t i) {
      //determine range direction
      std::int32_t d =
        0 > (delta(i, i + 1, inner_size, mcodes_ptr) - delta(i, i - 1, inner_size, mcodes_ptr)) ? -1
                                                                                                : 1;

      //find upper bound for the length of the range
      std::int32_t min_delta = delta(i, i - d, inner_size, mcodes_ptr);
      std::int32_t lmax = 2;
      while(delta(i, i + lmax * d, inner_size, mcodes_ptr) > min_delta)
      {
        lmax *= 2;
      }

      //binary search to find the lower bound
      std::int32_t l = 0;
      for(std::int32_t t = lmax / 2; t >= 1; t /= 2)
      {
        if(delta(i, i + (l + t) * d, inner_size, mcodes_ptr) > min_delta)
        {
          l += t;
        }
      }

      std::int32_t j = i + l * d;
      std::int32_t delta_node = delta(i, j, inner_size, mcodes_ptr);
      std::int32_t s = 0;
      FloatType div_factor = 2.f;
      //find the split postition using a binary search
      for(std::int32_t t = (std::int32_t)ceil(float32(l) / div_factor);;
          div_factor *= 2, t = (std::int32_t)ceil(float32(l) / div_factor))
      {
        if(delta(i, i + (s + t) * d, inner_size, mcodes_ptr) > delta_node)
        {
          s += t;
        }
        if(t == 1)
        {
          break;
        }
      }

      std::int32_t split = i + s * d + utilities::min(d, 0);
      // assign parent/child pointers
      if(utilities::min(i, j) == split)
      {
        //leaf
        parent_ptr[split + inner_size] = i;
        lchildren_ptr[i] = split + inner_size;
      }
      else
      {
        //inner node
        parent_ptr[split] = i;
        lchildren_ptr[i] = split;
      }

      if(utilities::max(i, j) == split + 1)
      {
        //leaf
        parent_ptr[split + inner_size + 1] = i;
        rchildren_ptr[i] = split + inner_size + 1;
      }
      else
      {
        parent_ptr[split + 1] = i;
        rchildren_ptr[i] = split + 1;
      }

      if(i == 0)
      {
        // flag the root
        parent_ptr[0] = -1;
      }
    });
}

//------------------------------------------------------------------------------
// Fetches a bounding box value synchronized with another thread's store.
// On the CPU, this is achieved with an acquire fence. This is only really
//   needed for non-x86 architectures with a weaker memory model (Power, ARM)
// On the GPU, we poll the bounding box values for a non-sentinel value.
template <typename ExecSpace, typename BBoxType>
AXOM_HOST_DEVICE static inline BBoxType sync_load(const BBoxType& box)
{
#ifdef AXOM_DEVICE_CODE

  using FloatType = typename BBoxType::CoordType;
  using PointType = typename BBoxType::PointType;

  constexpr int NDIMS = PointType::DIMENSION;

  PointType min_pt {BBoxType::InvalidMin};
  PointType max_pt {BBoxType::InvalidMax};

  #ifdef SPIN_BVH_DEBUG_MEMORY_HAZARD
  int nreads = 0;  // number of extra reads needed for a non-sentinel value
  #endif
  for(int dim = 0; dim < NDIMS; dim++)
  {
    // Cast to volatile so reads always hit L2$ or memory.
    volatile const FloatType& min_dim =
      reinterpret_cast<volatile const FloatType&>(box.getMin()[dim]);
    volatile const FloatType& max_dim =
      reinterpret_cast<volatile const FloatType&>(box.getMax()[dim]);

      // NOTE: There is a possibility for a read-after-write hazard, where the
      // uncached store of an AABB on one thread isn't visible when another
      // thread calls this method to read the value. However, this doesn't seem to
      // be an issue on Volta; the atomicAdd used to terminate the first thread
      // seems to correctly synchronize the prior atomic store operations for the
      // bounding box data.
      //
      // Just in case this changes, we poll for a non-sentinel value to be read
      // out. Naturally, this assumes that reads of sizeof(FloatType) don't tear.
  #ifdef SPIN_BVH_DEBUG_MEMORY_HAZARD
    while((min_pt[dim] = min_dim) == BBoxType::InvalidMin)
    {
      nreads++;
    }
    while((max_pt[dim] = max_dim) == BBoxType::InvalidMax)
    {
      nreads++;
    }
  #else
    while((min_pt[dim] = min_dim) == BBoxType::InvalidMin);
    while((max_pt[dim] = max_dim) == BBoxType::InvalidMax);
  #endif
  }

  #ifdef SPIN_BVH_DEBUG_MEMORY_HAZARD
  if(nreads > 0)
  {
    printf("Warning: needed %d extra reads for address %p\n", nreads, &box);
  }
  #endif

  return BBoxType {min_pt, max_pt};

#else  // AXOM_DEVICE_CODE
  std::atomic_thread_fence(std::memory_order_acquire);
  return box;
#endif
}

//------------------------------------------------------------------------------
// Writes a bounding box to memory, synchronized with another thread's read.
// On the CPU, this is achieved with a release fence.
// On the GPU, this function uses atomicExch to write a value directly to the
// L2 cache, thus avoiding potential cache coherency issues between threads.
template <typename ExecSpace, typename BBoxType>
AXOM_HOST_DEVICE static inline void sync_store(BBoxType& box, const BBoxType& value)
{
#if defined(AXOM_DEVICE_CODE) && defined(AXOM_USE_RAJA)
  using PointType = typename BBoxType::PointType;

  constexpr int NDIMS = PointType::DIMENSION;

  // Cast away the underlying const so we can directly modify the box data.
  PointType& min_pt = const_cast<PointType&>(box.getMin());
  PointType& max_pt = const_cast<PointType&>(box.getMax());

  for(int dim = 0; dim < NDIMS; dim++)
  {
    axom::atomicExchange<ExecSpace>(&(min_pt[dim]), value.getMin()[dim]);
    axom::atomicExchange<ExecSpace>(&(max_pt[dim]), value.getMax()[dim]);
  }
#else  // __CUDA_ARCH__ || __HIP_DEVICE_COMPILE__
  box = value;
  std::atomic_thread_fence(std::memory_order_release);
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
AXOM_HOST_DEVICE static inline int atomic_increment(int* addr)
{
#ifdef AXOM_USE_RAJA
  return axom::atomicAdd<ExecSpace>(addr, 1);
#else
  static_assert(std::is_same<ExecSpace, SEQ_EXEC>::value, "Only SEQ_EXEC supported without RAJA");

  // TODO: use atomic_ref?
  int old = *addr;
  (*addr)++;
  return old;
#endif  // AXOM_USE_RAJA
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
void propagate_aabbs(RadixTree<FloatType, NDIMS>& data, int allocatorID)
{
  AXOM_ANNOTATE_SCOPE("propagate_abbs");

  using BoxType = primal::BoundingBox<FloatType, NDIMS>;

  const int inner_size = data.m_inner_size;
  const int leaf_size = data.m_inner_size + 1;
  SLIC_ASSERT(leaf_size == data.m_size);

  // Pointers and vars are redeclared because I have a faint memory
  // of a huge amount of pain and suffering due so cuda
  // labda captures of pointers indide a struct. Bad memories
  // of random segfaults ........ be warned
  const auto lchildren_ptr = data.m_left_children.view();
  const auto rchildren_ptr = data.m_right_children.view();
  const auto parent_ptr = data.m_parents.view();
  const auto leaf_aabb_ptr = data.m_leaf_aabbs.view();

  const auto inner_aabb_ptr = data.m_inner_aabbs.view();
  for_all<ExecSpace>(inner_size, AXOM_LAMBDA(IndexType idx) { inner_aabb_ptr[idx] = BoxType {}; });

  Array<std::int32_t> counters(inner_size, inner_size, allocatorID);
  const auto counters_ptr = counters.view();

  for_all<ExecSpace>(
    leaf_size,
    AXOM_LAMBDA(std::int32_t i) {
      BoxType aabb = leaf_aabb_ptr[i];
      std::int32_t last_node = inner_size + i;
      std::int32_t current_node = parent_ptr[inner_size + i];

      while(current_node != -1)
      {
        // TODO: If RAJA atomics get memory ordering policies in the future,
        // we should look at replacing the sync_load/sync_stores by changing
        // the below atomic to an acquire/release atomic.
        std::int32_t old = atomic_increment<ExecSpace>(&(counters_ptr[current_node]));

        if(old == 0)
        {
          // first thread to get here kills itself
          return;
        }

        std::int32_t lchild = lchildren_ptr[current_node];
        std::int32_t rchild = rchildren_ptr[current_node];

        std::int32_t other_child = (lchild == last_node) ? rchild : lchild;

        primal::BoundingBox<FloatType, NDIMS> other_aabb;
        if(other_child >= inner_size)
        {
          other_aabb = leaf_aabb_ptr[other_child - inner_size];
        }
        else
        {
          other_aabb = sync_load<ExecSpace>(inner_aabb_ptr[other_child]);
        }
        aabb.addBox(other_aabb);

        // Store the final AABB for this internal node coherently.
        sync_store<ExecSpace>(inner_aabb_ptr[current_node], aabb);

        last_node = current_node;
        current_node = parent_ptr[current_node];
      }
    });
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename BoxIndexable, typename FloatType, int NDIMS>
void build_radix_tree(const BoxIndexable boxes,
                      int size,
                      primal::BoundingBox<FloatType, NDIMS>& bounds,
                      RadixTree<FloatType, NDIMS>& radix_tree,
                      FloatType scale_factor,
                      int allocatorID)
{
  AXOM_ANNOTATE_SCOPE("build_radix_tree");

  // sanity checks
  SLIC_ASSERT(size > 0);

  radix_tree.allocate(size, allocatorID);

  // copy so we don't reorder the input
  transform_boxes<ExecSpace>(boxes, radix_tree.m_leaf_aabbs.view(), size, scale_factor);

  // evaluate global bounds
  bounds = reduce<ExecSpace, FloatType, NDIMS>(radix_tree.m_leaf_aabbs, size);

  // sort aabbs based on morton code
  // original positions of the sorted morton codes.
  // allows us to gather / sort other arrays.
  get_mcodes<ExecSpace, FloatType, NDIMS>(radix_tree.m_leaf_aabbs, size, bounds, radix_tree.m_mcodes);
  sort_mcodes<ExecSpace>(radix_tree.m_mcodes, size, radix_tree.m_leafs);

  reorder<ExecSpace>(radix_tree.m_leafs, radix_tree.m_leaf_aabbs, size, allocatorID);

  build_tree<ExecSpace>(radix_tree);

  propagate_aabbs<ExecSpace>(radix_tree, allocatorID);
}

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */
#endif
