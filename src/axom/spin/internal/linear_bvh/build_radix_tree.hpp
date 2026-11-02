// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

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
#include "axom/spin/BucketedBVHGrouping.hpp"

#include "axom/spin/internal/linear_bvh/RadixTree.hpp"
#include "axom/spin/MortonIndex.hpp"
#include "axom/spin/HilbertIndex.hpp"
#include "axom/spin/SpaceFillingCurve.hpp"

#include <atomic>
#include <type_traits>

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA)
  // NOTE: uses the cub installation that is  bundled with RAJA
  #include "cub/device/device_radix_sort.cuh"
#endif

namespace axom::spin::internal::linear_bvh
{
//------------------------------------------------------------------------------
// Encodes normalized coordinates in [0,1] using the selected spatial ordering.
// Uses 16 bits per axis in 2D and 10 bits per axis in 3D.
template <typename FloatType, int Dims>
static inline AXOM_HOST_DEVICE std::uint32_t encode_spatial_key(
  const primal::Vector<FloatType, Dims>& point,
  SpaceFillingCurve ordering)
{
  using PointType = primal::Point<std::int32_t, Dims>;
  using axom::utilities::clampVal;

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

  const PointType cell(int_coords);
  return ordering == SpaceFillingCurve::Hilbert
    ? convertPointToHilbert<std::uint32_t>(cell, NUM_BITS_PER_DIM)
    : convertPointToMorton<std::uint32_t>(cell);
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

  for_all<ExecSpace>(size, [=] AXOM_HOST_DEVICE(std::int32_t i) {
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

    for_all<ExecSpace>(size, [=] AXOM_HOST_DEVICE(std::int32_t i) {
      const primal::BoundingBox<FloatType, NDIMS>& aabb = aabbs[i];
      min_coord.min(aabb.getMin()[dim]);
      max_coord.max(aabb.getMax()[dim]);
    });

    min_pt[dim] = min_coord.get();
    max_pt[dim] = max_coord.get();
  }

  // The "false" parameter disables the BoundingBox bounds-fixing logic when min_pt[i] > max_pt[i].
  // This happens when the user provides no input boxes for the BVH, in that case we don't want to
  // return a valid bounding box for the BVH bounds.
  return primal::BoundingBox<FloatType, NDIMS>(min_pt, max_pt, false);
#else
  static_assert(std::is_same<ExecSpace, SEQ_EXEC>::value, "Only SEQ_EXEC supported without RAJA");

  primal::BoundingBox<FloatType, NDIMS> global_bounds;

  for_all<ExecSpace>(size, [&](std::int32_t i) { global_bounds.addBox(aabbs[i]); });

  return global_bounds;
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
void get_spatial_keys(ArrayView<const primal::BoundingBox<FloatType, NDIMS>> aabbs,
                      std::int32_t size,
                      const primal::BoundingBox<FloatType, NDIMS>& bounds,
                      const ArrayView<std::uint32_t> keys,
                      SpaceFillingCurve ordering)
{
  AXOM_ANNOTATE_SCOPE("get_spatial_keys");

  primal::Vector<FloatType, NDIMS> extent, inv_extent, min_coord;

  extent = primal::Vector<FloatType, NDIMS>(bounds.getMax());
  extent -= primal::Vector<FloatType, NDIMS>(bounds.getMin());
  min_coord = primal::Vector<FloatType, NDIMS>(bounds.getMin());

  for(int i = 0; i < NDIMS; ++i)
  {
    inv_extent[i] = utilities::isNearlyEqual<FloatType>(extent[i], .0f) ? 0.f : 1.f / extent[i];
  }

  for_all<ExecSpace>(size, [=] AXOM_HOST_DEVICE(std::int32_t i) {
    const primal::BoundingBox<FloatType, NDIMS>& aabb = aabbs[i];

    // get the center and normalize it
    primal::Vector<FloatType, NDIMS> centroid(aabb.getCentroid());
    centroid = primal::Vector<FloatType, NDIMS>((centroid - min_coord).array() * inv_extent.array());
    keys[i] = encode_spatial_key(centroid, ordering);
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

  for_all<ExecSpace>(size, [=] AXOM_HOST_DEVICE(std::int32_t i) { iterator[i] = start + i * step; });
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

  for_all<ExecSpace>(size, [=] AXOM_HOST_DEVICE(std::int32_t i) {
    std::int32_t in_idx = indices[i];
    temp_v[i] = array_v[in_idx];
  });

  array = std::move(temp);
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_RAJA) && \
  ((RAJA_VERSION_MAJOR > 0) || ((RAJA_VERSION_MAJOR == 0) && (RAJA_VERSION_MINOR >= 12)))

template <typename ExecSpace>
void sort_spatial_keys(ArrayView<std::uint32_t> keys, std::int32_t size, ArrayView<std::int32_t> iter)
{
  AXOM_ANNOTATE_SCOPE("sort_spatial_keys");

  array_counting<ExecSpace>(iter, size, 0, 1);

  {
    AXOM_ANNOTATE_SCOPE("raja_stable_sort");
    ArrayView<std::uint32_t> keys_view(keys.data(), size);
    ArrayView<std::int32_t> iter_view(iter.data(), size);
    axom::stable_sort_pairs<ExecSpace>(keys_view, iter_view);
  }
}

#else

// fall back to std::stable_sort
template <typename ExecSpace>
void sort_spatial_keys(Array<std::uint32_t>& keys, std::int32_t size, const ArrayView<std::int32_t> iter)
{
  AXOM_ANNOTATE_SCOPE("sort_spatial_keys");

  array_counting<ExecSpace>(iter, size, 0, 1);

  {
    AXOM_ANNOTATE_SCOPE("cpu_sort");

    std::stable_sort(iter.begin(), iter.begin() + size, [&](std::int32_t i1, std::int32_t i2) {
      return keys[i1] < keys[i2];
    });
  }

  const int allocID = axom::execution_space<ExecSpace>::allocatorID();
  reorder<ExecSpace>(iter, keys, size, allocID);
}

#endif /* RAJA version 0.12.0 and above */

//------------------------------------------------------------------------------
template <typename IntType, typename KeyType>
AXOM_HOST_DEVICE IntType
delta(const IntType& a, const IntType& b, const IntType& inner_size, axom::ArrayView<KeyType> keys)
{
  bool tie = false;
  bool out_of_range = (b < 0 || b > inner_size);
  //still make the call but with a valid adderss
  const std::int32_t bb = (out_of_range) ? 0 : b;
  const std::uint32_t a_key = keys[a];
  const std::uint32_t b_key = keys[bb];
  //use xor to find where they differ
  std::uint32_t exor = a_key ^ b_key;
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
  const auto keys_ptr = data.m_keys.view();

  for_all<ExecSpace>(inner_size, [=] AXOM_HOST_DEVICE(std::int32_t i) {
    //determine range direction
    std::int32_t d =
      0 > (delta(i, i + 1, inner_size, keys_ptr) - delta(i, i - 1, inner_size, keys_ptr)) ? -1 : 1;

    //find upper bound for the length of the range
    std::int32_t min_delta = delta(i, i - d, inner_size, keys_ptr);
    std::int32_t lmax = 2;
    while(delta(i, i + lmax * d, inner_size, keys_ptr) > min_delta)
    {
      lmax *= 2;
    }

    //binary search to find the lower bound
    std::int32_t l = 0;
    for(std::int32_t t = lmax / 2; t >= 1; t /= 2)
    {
      if(delta(i, i + (l + t) * d, inner_size, keys_ptr) > min_delta)
      {
        l += t;
      }
    }

    std::int32_t j = i + l * d;
    std::int32_t delta_node = delta(i, j, inner_size, keys_ptr);
    std::int32_t s = 0;
    FloatType div_factor = 2.f;
    //find the split postition using a binary search
    for(std::int32_t t = (std::int32_t)ceil(float32(l) / div_factor);;
        div_factor *= 2, t = (std::int32_t)ceil(float32(l) / div_factor))
    {
      if(delta(i, i + (s + t) * d, inner_size, keys_ptr) > delta_node)
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
  using FloatType = typename BBoxType::CoordType;
  using PointType = typename BBoxType::PointType;

  constexpr int NDIMS = PointType::DIMENSION;

  PointType min_pt {BBoxType::InvalidMin};
  PointType max_pt {BBoxType::InvalidMax};

  // Read each component atomically to avoid a racy read of the box struct.
  // This mirrors the device path's "poll until non-sentinel" approach.
  for(int dim = 0; dim < NDIMS; dim++)
  {
    while((min_pt[dim] = axom::atomicLoad<ExecSpace>(const_cast<FloatType*>(&box.getMin()[dim]))) ==
          BBoxType::InvalidMin);
    while((max_pt[dim] = axom::atomicLoad<ExecSpace>(const_cast<FloatType*>(&box.getMax()[dim]))) ==
          BBoxType::InvalidMax);
  }

  return BBoxType {min_pt, max_pt};
#endif
}

//------------------------------------------------------------------------------
// Writes a bounding box to memory, synchronized with another thread's read.
// On the CPU, we use atomic stores for each component to avoid data races.
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
  using PointType = typename BBoxType::PointType;
  constexpr int NDIMS = PointType::DIMENSION;

  // Cast away the underlying const so we can directly modify the box data.
  PointType& min_pt = const_cast<PointType&>(box.getMin());
  PointType& max_pt = const_cast<PointType&>(box.getMax());

  for(int dim = 0; dim < NDIMS; dim++)
  {
    axom::atomicStore<ExecSpace>(&(min_pt[dim]), value.getMin()[dim]);
    axom::atomicStore<ExecSpace>(&(max_pt[dim]), value.getMax()[dim]);
  }
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
  for_all<ExecSpace>(inner_size,
                     [=] AXOM_HOST_DEVICE(IndexType idx) { inner_aabb_ptr[idx] = BoxType {}; });

  Array<std::int32_t> counters(inner_size, inner_size, allocatorID);
  const auto counters_ptr = counters.view();

  for_all<ExecSpace>(leaf_size, [=] AXOM_HOST_DEVICE(std::int32_t i) {
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
/*!
 * \brief Allocate and sort the leaves of a radix tree.
 *
 * \note This does not construct the radix-tree topology. Call
 * `build_radix_tree_hierarchy()` after this function to complete the
 * tree.
 */
template <typename ExecSpace, typename BoxIndexable, typename FloatType, int NDIMS>
void initialize_spatially_sorted_leaves(const BoxIndexable boxes,
                                        int size,
                                        primal::BoundingBox<FloatType, NDIMS>& bounds,
                                        RadixTree<FloatType, NDIMS>& radix_tree,
                                        FloatType scale_factor,
                                        int allocatorID,
                                        SpaceFillingCurve ordering)
{
  radix_tree.allocate(size, allocatorID);

  // copy so we don't reorder the input
  transform_boxes<ExecSpace>(boxes, radix_tree.m_leaf_aabbs.view(), size, scale_factor);

  // evaluate global bounds
  bounds = reduce<ExecSpace, FloatType, NDIMS>(radix_tree.m_leaf_aabbs, size);

  // Keep original box IDs when sorting so queries can return input indices.
  get_spatial_keys<ExecSpace, FloatType, NDIMS>(radix_tree.m_leaf_aabbs,
                                                size,
                                                bounds,
                                                radix_tree.m_keys,
                                                ordering);
  sort_spatial_keys<ExecSpace>(radix_tree.m_keys, size, radix_tree.m_leafs);

  reorder<ExecSpace>(radix_tree.m_leafs, radix_tree.m_leaf_aabbs, size, allocatorID);
}

///@{
//! @name Functions which construct bucket ranges from a sorted list of keys,

using BucketRange = std::pair<std::int32_t, std::int32_t>;

//------------------------------------------------------------------------------
/*!
 * \brief Collect indices for primitives in each bucket based on range
 *
 * The resulting bucket ranges partition the sorted list of keys and retain 
 * spatial-key order, so their first primitive code can serve as a representative
 * key for the existing Karras hierarchy builder. 
 * 
 * \note Although the BucketedLinearBVH implementation is general enough to
 *   support non-contiugous bucket membership, all currently supported grouping
 *   strategies assume contiguous bucket membership. 
 */
template <typename ExecSpace, typename FloatType, int NDIMS>
void build_bucket_tree_from_ranges(const RadixTree<FloatType, NDIMS>& primitive_tree,
                                   const Array<BucketRange>& bucket_ranges,
                                   RadixTree<FloatType, NDIMS>& bucket_tree,
                                   Array<IndexType>& bucket_offsets,
                                   Array<std::int32_t>& bucket_primitive_ids,
                                   int allocatorID)
{
  AXOM_ANNOTATE_SCOPE("bucket_materialization");

  // Collect data from the original specialized BucketedLinearBVH
  const auto primitive_aabbs = primitive_tree.m_leaf_aabbs.view();
  const auto primitive_keys = primitive_tree.m_keys.view();
  const auto primitive_ids = primitive_tree.m_leafs.view();
  const std::int32_t num_primitives = primitive_tree.m_size;
  const auto num_buckets = static_cast<std::int32_t>(bucket_ranges.size());

  // Radix tree for buckets: each leaf stores a bucket union AABB, a
  // representative key (the first in the range), and its bucket ID.
  bucket_tree.allocate(num_buckets, allocatorID);
  auto bucket_tree_aabbs = bucket_tree.m_leaf_aabbs.view();
  auto bucket_tree_keys = bucket_tree.m_keys.view();
  auto bucket_tree_ids = bucket_tree.m_leafs.view();

  // This array maps each bucket ID to its first position in the flattened primitive-ID array.
  bucket_offsets = Array<IndexType>(num_buckets + 1, num_buckets + 1, allocatorID);
  auto bucket_offsets_view = bucket_offsets.view();

  // This array maps each position in the flattened bucket list to an original primitive ID.
  bucket_primitive_ids = Array<std::int32_t>(num_primitives, num_primitives, allocatorID);
  auto bucket_primitive_ids_view = bucket_primitive_ids.view();

  // Bucket ranges partition the sorted primitive list. Each iteration writes
  // a distinct bucket entry and a disjoint primitive-ID subrange.
  const auto bucket_ranges_view = bucket_ranges.view();
  for_all<ExecSpace>(num_buckets, [=] AXOM_HOST_DEVICE(std::int32_t bucket) {
    const std::int32_t begin = bucket_ranges_view[bucket].first;
    const std::int32_t end = bucket_ranges_view[bucket].second;

    // Combine all boxes for primitives in this Morton-prefix range.
    primal::BoundingBox<FloatType, NDIMS> bucket_aabb;
    for(std::int32_t primitive = begin; primitive < end; ++primitive)
    {
      bucket_aabb.addBox(primitive_aabbs[primitive]);

      // Preserve the original primitive ID at this spatially sorted position.
      bucket_primitive_ids_view[primitive] = primitive_ids[primitive];
    }

    // Add bucket parameters to radix tree
    bucket_tree_aabbs[bucket] = bucket_aabb;
    bucket_tree_keys[bucket] = primitive_keys[begin];
    bucket_tree_ids[bucket] = bucket;

    // Store the offset for this bucket
    bucket_offsets_view[bucket] = begin;
  });

  // Store the offset for the last bucket
  bucket_offsets_view[num_buckets] = num_primitives;
}

//------------------------------------------------------------------------------
/*!
 * \brief Append shared-spatial-key-prefix bucket ranges for a sorted key range.
 *
 * This grouping strategy defines a bucket as the largest range of sorted primitive keys
 * which share a bit prefix. If this range exceeds \a max_bucket_size, it is split
 * along its next key bit. Ties are broken by splitting the range sequentially.
 * 
 * \param[in] keys The sorted list of spatial keys
 * \param[in] begin, end The range of primitive indices to accept or split
 * \param[in] bit The bit position for a possible split
 * \param[in] max_bucket_size The maximum number of primitives allowed in a bucket
 * \param[out] bucket_ranges The output range for the primitives in each bucket 
 */
inline void append_shared_prefix_buckets(const ArrayView<const std::uint32_t> keys,
                                         std::int32_t begin,
                                         std::int32_t end,
                                         int bit,
                                         std::int32_t max_bucket_size,
                                         Array<BucketRange>& bucket_ranges)
{
  // Check if the current bucket is bigger than allowed
  const std::int32_t bucket_size = end - begin;
  if(bucket_size <= max_bucket_size)
  {
    bucket_ranges.push_back({begin, end});
    return;
  }

  if(bit < 0)
  {
    // If the current range has equal spatial keys, split it sequentially
    std::int32_t chunk_begin = begin;
    while(chunk_begin < end)
    {
      const std::int32_t chunk_size = utilities::min(max_bucket_size, end - chunk_begin);
      const std::int32_t chunk_end = chunk_begin + chunk_size;
      bucket_ranges.push_back({chunk_begin, chunk_end});
      chunk_begin = chunk_end;
    }
    return;
  }

  // Spatial-key sorting ensures the zero-valued side of this bit precedes the
  // one-valued side. Find their boundary without examining every primitive.
  const std::uint32_t mask = std::uint32_t {1} << bit;
  std::int32_t split = begin;
  std::int32_t search_end = end;
  while(split < search_end)
  {
    // Do a binary search to find where the bit flips, e.g.
    //  for bit == 1: 0000 0001 | 0010 0011
    const std::int32_t middle = split + (search_end - split) / 2;
    if((keys[middle] & mask) == 0)
    {
      split = middle + 1;
    }
    else
    {
      search_end = middle;
    }
  }

  // An empty child means every primitive has the same value at this bit.
  // Continue to the next bit and attempt to split the group there instead.
  if(split == begin || split == end)
  {
    append_shared_prefix_buckets(keys, begin, end, bit - 1, max_bucket_size, bucket_ranges);
    return;
  }

  // Recurse into the two halves of the input range, split at `bit`
  append_shared_prefix_buckets(keys, begin, split, bit - 1, max_bucket_size, bucket_ranges);
  append_shared_prefix_buckets(keys, split, end, bit - 1, max_bucket_size, bucket_ranges);
}

//------------------------------------------------------------------------------
/*!
 * \brief Form primitive buckets based on shared prefixes in the spatial keys of each primitive
 *
 * This policy will result in many buckets with fewer than \a max_bucket_size primitives,
 * but the buckets that result have very high spatial coherence, as primtiives *must* have
 * a shared prefix key to be bucketed together.  
 * 
 * \param[in] primitive_tree The radix tree with spatially sorted keys for each primitive
 * \param[in] max_bucket_size The maximum number of primitives allowed in a bucket
 * \param[out] bucket_tree The radix tree with spatially sorted keys for each bucket
 * \param[out] bucket_offsets Map from bucket ID to first primitive ID
 * \param[out] primitive_bucket_ids Map from positions in the bucket list to a primitive ID
 */
template <typename ExecSpace, typename FloatType, int NDIMS>
void build_shared_prefix_buckets(const RadixTree<FloatType, NDIMS>& primitive_tree,
                                 std::int32_t max_bucket_size,
                                 RadixTree<FloatType, NDIMS>& bucket_tree,
                                 Array<IndexType>& bucket_offsets,
                                 Array<std::int32_t>& bucket_primitive_ids,
                                 int allocatorID)
{
  AXOM_ANNOTATE_SCOPE("create_shared_prefix_buckets");
  SLIC_ASSERT(primitive_tree.m_size > 0);
  SLIC_ASSERT(max_bucket_size > 0);

  Array<BucketRange> bucket_ranges;
  append_shared_prefix_buckets(primitive_tree.m_keys.view(),
                               /* begin */ 0,
                               /* end */ primitive_tree.m_size,
                               /* bit */ 31,
                               max_bucket_size,
                               bucket_ranges);

  build_bucket_tree_from_ranges<ExecSpace>(primitive_tree,
                                           bucket_ranges,
                                           bucket_tree,
                                           bucket_offsets,
                                           bucket_primitive_ids,
                                           allocatorID);
}

//------------------------------------------------------------------------------
/*!
 * \brief Form primitive buckets sequentially by primitive spatial keys
 *
 * This policy will result in every bucket (except the last) having exactly \a max_bucket_size
 * primitives, This minimizes the number of buckets, but may result in low spatial coherence, 
 * as bucketed primtiives can have far-away spatial keys.
 * 
 * \param[in] primitive_tree The radix tree with spatially sorted keys for each primitive
 * \param[in] max_bucket_size The maximum number of primitives allowed in a bucket
 * \param[out] bucket_tree The radix tree with spatially sorted keys for each bucket
 * \param[out] bucket_offsets Maps each bucket ID to its range in \a bucket_primitive_ids.
 * \param[out] bucket_primitive_ids Flattened original primitive IDs for all buckets.
 */
template <typename ExecSpace, typename FloatType, int NDIMS>
void build_sequential_buckets(const RadixTree<FloatType, NDIMS>& primitive_tree,
                              std::int32_t max_bucket_size,
                              RadixTree<FloatType, NDIMS>& bucket_tree,
                              Array<IndexType>& bucket_offsets,
                              Array<std::int32_t>& bucket_primitive_ids,
                              int allocatorID)
{
  AXOM_ANNOTATE_SCOPE("create_sequential_buckets");
  SLIC_ASSERT(primitive_tree.m_size > 0);
  SLIC_ASSERT(max_bucket_size > 0);

  Array<BucketRange> bucket_ranges;
  bucket_ranges.reserve((primitive_tree.m_size + max_bucket_size - 1) / max_bucket_size);
  for(std::int32_t begin = 0; begin < primitive_tree.m_size; begin += max_bucket_size)
  {
    bucket_ranges.push_back({begin, utilities::min(begin + max_bucket_size, primitive_tree.m_size)});
  }

  build_bucket_tree_from_ranges<ExecSpace>(primitive_tree,
                                           bucket_ranges,
                                           bucket_tree,
                                           bucket_offsets,
                                           bucket_primitive_ids,
                                           allocatorID);
}
///@}

//------------------------------------------------------------------------------
/*!
 * \brief Construct a radix-tree topology from spatially sorted leaves.
 *
 * \pre \a radix_tree was populated by `initialize_spatially_sorted_leaves()` or has
 *      equivalent sorted leaf AABBs, spatial keys, and leaf IDs (such as from buckets)
 */
template <typename ExecSpace, typename FloatType, int NDIMS>
void build_radix_tree_hierarchy(RadixTree<FloatType, NDIMS>& radix_tree, int allocatorID)
{
  build_tree<ExecSpace>(radix_tree);

  propagate_aabbs<ExecSpace>(radix_tree, allocatorID);
}

//------------------------------------------------------------------------------
/*!
 * \brief Build a bucketed radix tree from primitive bounding boxes.
 *
 * Primitive leaves are scaled, spatially sorted according to \a ordering, 
 * then partitioned according to \a grouping. The resulting bucket leaves
 * are used to assemble the final radix-tree hierarchy.
 *
 * \param[in] boxes Indexed collection of input primitive bounding boxes.
 * \param[in] size Number of primitive bounding boxes in \a boxes.
 * \param[out] bounds Bounds of the scaled input primitive bounding boxes.
 * \param[out] bucket_tree Radix tree whose leaves represent primitive buckets.
 * \param[out] bucket_offsets Maps each bucket ID to its range in \a bucket_primitive_ids.
 * \param[out] bucket_primitive_ids Flattened original primitive IDs for all buckets.
 * \param[in] max_bucket_size Maximum number of primitives in a bucket.
 * \param[in] scale_factor Factor applied to each primitive bounding box.
 * \param[in] allocatorID AllocatFor for intermediate and output tree storage.
 * \param[in] ordering Space-filling curve used to spatially order primitives.
 * \param[in] grouping Policy used to group spatially sorted primitives into buckets.
 *
 * \pre Currently only supported for CPU execution spaces
 */
template <typename ExecSpace, typename BoxIndexable, typename FloatType, int NDIMS>
void build_bucketed_radix_tree(const BoxIndexable boxes,
                               int size,
                               primal::BoundingBox<FloatType, NDIMS>& bounds,
                               RadixTree<FloatType, NDIMS>& bucket_tree,
                               Array<IndexType>& bucket_offsets,
                               Array<std::int32_t>& bucket_primitive_ids,
                               std::int32_t max_bucket_size,
                               FloatType scale_factor,
                               int allocatorID,
                               SpaceFillingCurve ordering,
                               BucketedBVHGrouping grouping)
{
  AXOM_ANNOTATE_SCOPE("build_bucketed_radix_tree");
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  static_assert(
    std::is_same_v<ExecSpace, axom::SEQ_EXEC> || std::is_same_v<ExecSpace, axom::OMP_EXEC>,
    "Bucketed radix-tree construction currently supports only SEQ_EXEC and OMP_EXEC.");
#else
  static_assert(std::is_same_v<ExecSpace, axom::SEQ_EXEC>,
                "Bucketed radix-tree construction currently supports only SEQ_EXEC and OMP_EXEC.");
#endif
  SLIC_ASSERT(size > 0);

  RadixTree<FloatType, NDIMS> primitive_tree;
  initialize_spatially_sorted_leaves<ExecSpace>(boxes,
                                                size,
                                                bounds,
                                                primitive_tree,
                                                scale_factor,
                                                allocatorID,
                                                ordering);
  if(grouping == BucketedBVHGrouping::SharedPrefix)
  {
    build_shared_prefix_buckets<ExecSpace>(primitive_tree,
                                           max_bucket_size,
                                           bucket_tree,
                                           bucket_offsets,
                                           bucket_primitive_ids,
                                           allocatorID);
  }
  else
  {
    SLIC_ASSERT(grouping == BucketedBVHGrouping::Sequential);
    build_sequential_buckets<ExecSpace>(primitive_tree,
                                        max_bucket_size,
                                        bucket_tree,
                                        bucket_offsets,
                                        bucket_primitive_ids,
                                        allocatorID);
  }

  if(bucket_tree.m_size == 1)
  {
    // The binary hierarchy requires two leaves. Add an empty bucket whose
    // invalid AABB is skipped by traversal and whose offset range is empty.
    using BoxType = primal::BoundingBox<FloatType, NDIMS>;
    const BoxType valid_bucket_aabb = bucket_tree.m_leaf_aabbs[0];
    const auto bucket_key = bucket_tree.m_keys[0];

    RadixTree<FloatType, NDIMS> padded_bucket_tree;
    padded_bucket_tree.allocate(2, allocatorID);

    // Create an empty box
    BoxType empty_bucket_aabb;
    empty_bucket_aabb.clear();

    // Populate the padded tree
    padded_bucket_tree.m_leaf_aabbs[0] = valid_bucket_aabb;
    padded_bucket_tree.m_leaf_aabbs[1] = empty_bucket_aabb;
    padded_bucket_tree.m_keys[0] = bucket_key;
    padded_bucket_tree.m_keys[1] = bucket_key;
    padded_bucket_tree.m_leafs[0] = 0;
    padded_bucket_tree.m_leafs[1] = 1;
    bucket_tree = std::move(padded_bucket_tree);

    // Populate the indirection arrays for the padded tree
    Array<IndexType> padded_offsets(3, 3, allocatorID);
    padded_offsets[0] = 0;
    padded_offsets[1] = size;
    padded_offsets[2] = size;
    bucket_offsets = std::move(padded_offsets);
  }

  build_radix_tree_hierarchy<ExecSpace>(bucket_tree, allocatorID);
}

//------------------------------------------------------------------------------
/*!
 * \brief Build a radix tree from primitive bounding boxes.
 *
 * The primitive bounding boxes are scaled, spatially sorted using \a ordering,
 * and used to construct the radix-tree topology and internal bounding boxes.
 *
 * \param[in] boxes Indexed collection of input primitive bounding boxes.
 * \param[in] size Number of primitive bounding boxes in \a boxes.
 * \param[out] bounds Bounds of the scaled input primitive bounding boxes.
 * \param[out] radix_tree Radix tree populated with primitive leaves and hierarchy.
 * \param[in] scale_factor Factor applied to each primitive bounding box.
 * \param[in] allocatorID Allocator for radix-tree storage.
 * \param[in] ordering Space-filling curve used to spatially order primitives.
 *
 * \pre \a size is positive.
 */
template <typename ExecSpace, typename BoxIndexable, typename FloatType, int NDIMS>
void build_radix_tree(const BoxIndexable boxes,
                      int size,
                      primal::BoundingBox<FloatType, NDIMS>& bounds,
                      RadixTree<FloatType, NDIMS>& radix_tree,
                      FloatType scale_factor,
                      int allocatorID,
                      SpaceFillingCurve ordering = SpaceFillingCurve::Morton)
{
  AXOM_ANNOTATE_SCOPE("build_radix_tree");

  // sanity checks
  SLIC_ASSERT(size > 0);

  initialize_spatially_sorted_leaves<ExecSpace>(boxes,
                                                size,
                                                bounds,
                                                radix_tree,
                                                scale_factor,
                                                allocatorID,
                                                ordering);
  build_radix_tree_hierarchy<ExecSpace>(radix_tree, allocatorID);
}

}  // namespace axom::spin::internal::linear_bvh
