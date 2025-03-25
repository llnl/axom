// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_FlatMap_Util_HPP
#define Axom_Core_FlatMap_Util_HPP

#include "axom/config.hpp"
#include "axom/core/FlatMap.hpp"

namespace axom
{
namespace detail
{

struct SpinLock
{
  int value {0};

  AXOM_HOST_DEVICE bool tryLock()
  {
    int still_locked = 0;
#if defined(__HIP_DEVICE_COMPILE__)
    still_locked = __hip_atomic_exchange(&value, 1, __ATOMIC_ACQUIRE, __HIP_MEMORY_SCOPE_AGENT);
#elif defined(AXOM_USE_RAJA) && defined(__CUDA_ARCH__)
    still_locked = RAJA::atomicExchange<RAJA::cuda_atomic>(&value, 1);
    // We really want an acquire-fenced atomic here
    __threadfence();
#elif defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
    still_locked = RAJA::atomicExchange<RAJA::omp_atomic>(&value, 1);
    std::atomic_thread_fence(std::memory_order_acquire);
#endif
    return !still_locked;
  }

  AXOM_HOST_DEVICE void unlock()
  {
#if defined(__HIP_DEVICE_COMPILE__)
    __hip_atomic_exchange(&value, 0, __ATOMIC_RELEASE, __HIP_MEMORY_SCOPE_AGENT);
#elif defined(AXOM_USE_RAJA) && defined(__CUDA_ARCH__)
    // We really want a release-fenced atomic here
    __threadfence();
    RAJA::atomicExchange<RAJA::cuda_atomic>(&value, 0);
#elif defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
    std::atomic_thread_fence(std::memory_order_release);
    RAJA::atomicExchange<RAJA::omp_atomic>(&value, 0);
#else
    value = 0;
#endif
  }
};

}  // namespace detail

template <typename KeyType, typename ValueType, typename Hash>
template <typename ExecSpace>
auto FlatMap<KeyType, ValueType, Hash>::create(ArrayView<KeyType> keys,
                                               ArrayView<ValueType> values,
                                               Allocator allocator) -> FlatMap
{
  assert(keys.size() == values.size());

  IndexType num_elems = keys.size();

  FlatMap new_map(allocator);
  new_map.reserve(num_elems);

  using RajaReduce = typename axom::execution_space<ExecSpace>::reduce_policy;
  using HashResult = typename Hash::result_type;
  using GroupBucket = detail::flat_map::GroupBucket;

  // Grab some needed internal fields from the flat map.
  // We're going to be constructing metadata and the K-V pairs directly
  // in-place.
  int ngroups_pow_2 = new_map.m_numGroups2;
  const auto meta_group = new_map.m_metadata.view();
  const auto buckets = new_map.m_buckets.view();

  // Construct an array of locks per-group. This guards metadata updates for
  // each insertion.
  IndexType num_groups = 1 << ngroups_pow_2;
  Array<detail::SpinLock> lock_vec(num_groups, num_groups, allocator.get());
  const auto group_locks = lock_vec.view();

  // Add a counter for failed inserts.
#ifdef AXOM_USE_RAJA
  RAJA::ReduceSum<RajaReduce, int> num_failed_inserts(0);
#else
  int num_failed_inserts = 0;
#endif

  for_all<ExecSpace>(
    keys.size(),
    AXOM_LAMBDA(IndexType idx) {
      // Hash keys.
      auto hash = Hash {}(keys[idx]);

      // We use the k MSBs of the hash as the initial group probe point,
      // where ngroups = 2^k.
      int bitshift_right = ((CHAR_BIT * sizeof(HashResult)) - ngroups_pow_2);
      HashResult curr_group = hash >> bitshift_right;
      curr_group &= ((1 << ngroups_pow_2) - 1);

      std::uint8_t hash_8 = static_cast<std::uint8_t>(hash);

      int group_index = -1;
      int slot_index = -1;
      int iteration = 0;
      while(iteration < meta_group.size())
      {
        // Get the overflow state for the current group.
        bool overflow_group = meta_group[curr_group].template getMaybeOverflowed<true>(hash_8);
        bool group_locked = false;

        // Try to lock the group if we haven't overflowed.
        if(!overflow_group)
        {
          group_locked = group_locks[curr_group].tryLock();
        }

        if(group_locked)
        {
          // Re-check overflow state now that we've locked the group.
          // Other threads may have updated the overflow state before
          // actaully taking the lock.
          overflow_group = meta_group[curr_group].template getMaybeOverflowed<true>(hash_8);
          if(!overflow_group)
          {
            slot_index = meta_group[curr_group].getEmptyBucket();
            if(slot_index != GroupBucket::InvalidSlot)
            {
              // We have an empty slot in this group. Place an element.
              group_index = curr_group;
              meta_group[curr_group].template setBucket<true>(slot_index, hash_8);
            }
            else
            {
              // Group is full. Set overflow bit for the group.
              meta_group[curr_group].template setOverflow<true>(hash_8);
            }
          }
          // Unlock group once we're done.
          group_locks[curr_group].unlock();
        }

        // Either we got an empty slot to insert in, or we overflowed
        // on a group we waited on.
        if(group_index != -1)
        {
          // We have an empty slot. Exit.
          break;
        }
        else if(overflow_group)
        {
          // Move to next group.
          curr_group = (curr_group + LookupPolicy {}.getNext(iteration)) % meta_group.size();
          iteration++;
        }
      }

      if(group_index != -1)
      {
        IndexType bucket_index = group_index * GroupBucket::Size + slot_index;
        new(&buckets[bucket_index]) KeyValuePair(keys[idx], values[idx]);
      }
      else
      {
        // TODO: Is this even possible?
        num_failed_inserts += 1;
      }
    });

  new_map.m_size = keys.size();
  new_map.m_loadCount = keys.size();

#ifdef AXOM_DEBUG
  if(num_failed_inserts > 0)
  {
    axom::fmt::print("FlatMap::create: {} elements failed to insert\n");
  }
#endif

  return new_map;
}

}  // namespace axom

#endif
