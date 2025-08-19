// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_FlatMap_Util_HPP
#define Axom_Core_FlatMap_Util_HPP

#include "axom/config.hpp"
#include "axom/core/FlatMap.hpp"
#include "axom/core/execution/reductions.hpp"

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

/*!
 * \class KVPairIterator
 * \brief Implements a zip-iterator concept for a key-value pair.
 */
template <typename KeyType, typename ValueType>
class KVPairIterator : public IteratorBase<KVPairIterator<KeyType, ValueType>, IndexType>
{
private:
  using BaseType = IteratorBase<KVPairIterator<KeyType, ValueType>, IndexType>;
  using KeyIterator = const KeyType*;
  using ValueIterator = ValueType*;

public:
  // Iterator traits required to satisfy LegacyRandomAccessIterator concept
  // before C++20
  // See: https://en.cppreference.com/w/cpp/iterator/iterator_traits
  using difference_type = IndexType;
  using value_type = std::pair<const KeyType, ValueType>;
  using reference = ValueType&;
  using pointer = ValueType*;
  using iterator_category = std::random_access_iterator_tag;

  KVPairIterator() = default;

  AXOM_HOST_DEVICE KVPairIterator(KeyIterator keyIter, ValueIterator valueIter)
    : m_keyIter {keyIter}
    , m_valueIter {valueIter}
  { }

  /**
   * \brief Returns the current iterator value
   */
  AXOM_HOST_DEVICE value_type operator*() const { return {*m_keyIter, *m_valueIter}; }

  AXOM_HOST_DEVICE value_type operator[](IndexType idx) const { return *(*this + idx); }

protected:
  /** Implementation of advance() as required by IteratorBase */
  AXOM_HOST_DEVICE void advance(IndexType n)
  {
    BaseType::m_pos += n;
    m_keyIter += n;
    m_valueIter += n;
  }

private:
  KeyIterator m_keyIter {nullptr};
  ValueIterator m_valueIter {nullptr};
};

}  // namespace detail

template <typename KeyType, typename ValueType, typename Hash>
template <typename ExecSpace>
auto FlatMap<KeyType, ValueType, Hash>::create(ArrayView<KeyType> keys,
                                               ArrayView<ValueType> values,
                                               Allocator allocator) -> FlatMap
{
  assert(keys.size() == values.size());

  const IndexType num_elems = keys.size();
  detail::KVPairIterator zip_iterator {keys.data(), values.data()};

  FlatMap new_map(allocator);
  new_map.insert<ExecSpace>(zip_iterator, zip_iterator + num_elems);

  return new_map;
}

template <typename KeyType, typename ValueType, typename Hash>
template <typename ExecSpace, typename InputIt>
void FlatMap<KeyType, ValueType, Hash>::insert(InputIt kv_begin, InputIt kv_end)
{
  static_assert(std::is_base_of<std::random_access_iterator_tag,
                                typename std::iterator_traits<InputIt>::iterator_category>::value,
                "InputIt must be a random-access iterator for batched construction");

  using HashResult = typename Hash::result_type;
  using GroupBucket = detail::flat_map::GroupBucket;

  // Assume that all elements will be inserted into an empty slot.
  IndexType num_elems = std::distance(kv_begin, kv_end);
  this->reserve(this->size() + num_elems);

  // Grab some needed internal fields from the flat map.
  // We're going to be constructing metadata and the K-V pairs directly
  // in-place.
  const int ngroups_pow_2 = this->m_numGroups2;
  const auto meta_group = this->m_metadata.view();
  const auto buckets = this->m_buckets.view();

  // Construct an array of locks per-group. This guards metadata updates for
  // each insertion.
  const IndexType num_groups = 1 << ngroups_pow_2;
  Array<detail::SpinLock> lock_vec(num_groups, num_groups, this->m_allocator.getID());
  const auto group_locks = lock_vec.view();

  // Map bucket slots to k-v pair indices. This is used to deduplicate pairs
  // with the same key value.
  Array<IndexType> key_index_dedup_vec(0, 0, this->m_allocator.getID());
  key_index_dedup_vec.resize(num_groups * GroupBucket::Size, -1);
  const auto key_index_dedup = key_index_dedup_vec.view();

  // Map k-v pair indices to bucket slots. This is essentially the inverse of
  // the above mapping.
  Array<IndexType> key_index_to_bucket_vec(num_elems, num_elems, this->m_allocator.getID());
  const auto key_index_to_bucket = key_index_to_bucket_vec.view();

  for_all<ExecSpace>(
    num_elems,
    AXOM_LAMBDA(IndexType idx) {
      // Hash keys.
      auto hash = Hash {}(kv_begin[idx].first);

      // We use the k MSBs of the hash as the initial group probe point,
      // where ngroups = 2^k.
      int bitshift_right = ((CHAR_BIT * sizeof(HashResult)) - ngroups_pow_2);
      HashResult curr_group = hash >> bitshift_right;
      curr_group &= ((1 << ngroups_pow_2) - 1);

      std::uint8_t hash_8 = static_cast<std::uint8_t>(hash);

      IndexType duplicate_bucket_index = -1;
      IndexType empty_bucket_index = -1;
      int iteration = 0;
      while(iteration < meta_group.size())
      {
        // Try to lock the group. We do this in a non-blocking manner to avoid
        // intra-warp progress hazards.
        bool group_locked = group_locks[curr_group].tryLock();

        if(group_locked)
        {
          // Every bucket visit - check prior filled buckets for duplicate
          // keys.
          int empty_slot_index =
            meta_group[curr_group].visitHashOrEmptyBucket(hash_8, [&](int matching_slot) {
              IndexType bucket_index = curr_group * GroupBucket::Size + matching_slot;

              if(kv_begin[key_index_dedup[bucket_index]].first == kv_begin[idx].first)
              {
                // Highest-indexed kv pair wins.
                axom::atomicMax<ExecSpace>(&key_index_dedup[bucket_index], idx);
                key_index_to_bucket[idx] = bucket_index;
                duplicate_bucket_index = bucket_index;
              }
            });

          if(duplicate_bucket_index == -1)
          {
            if(empty_slot_index == GroupBucket::InvalidSlot)
            {
              // Group is full. Set overflow bit for the group.
              meta_group[curr_group].template setOverflow<true>(hash_8);
            }
            else
            {
              // Got to end of probe sequence without a duplicate.
              // Update empty bucket index.
              empty_bucket_index = curr_group * GroupBucket::Size + empty_slot_index;
              meta_group[curr_group].template setBucket<true>(empty_slot_index, hash_8);
              key_index_dedup[empty_bucket_index] = idx;
              key_index_to_bucket[idx] = empty_bucket_index;
            }
          }
          // Unlock group once we're done.
          group_locks[curr_group].unlock();

          if(duplicate_bucket_index != -1 || empty_bucket_index != -1)
          {
            // We've found an empty slot or a duplicate key to place the
            // value at. Empty slots should only occur at the end of the
            // probe sequence, since we're only inserting.
            break;
          }
          else
          {
            // Move to next group.
            curr_group = (curr_group + LookupPolicy {}.getNext(iteration)) % meta_group.size();
            iteration++;
          }
        }
      }
    });

  // Add a counter for duplicated inserts.
  axom::ReduceSum<ExecSpace, IndexType> total_inserts(0);

  // Using key-deduplication map, assign unique k-v pairs to buckets.
  for_all<ExecSpace>(
    num_elems,
    AXOM_LAMBDA(IndexType kv_idx) {
      IndexType bucket_idx = key_index_to_bucket[kv_idx];
      IndexType winning_idx = key_index_dedup[bucket_idx];
      // Place k-v pair at bucket_idx.
      if(kv_idx == winning_idx)
      {
#if defined(__CUDA_ARCH__)
        // HACK: std::pair constructor is not host-device annotated, but CUDA
        // requires passing in --expt-relaxed-constexpr for it to work.
        // Instead of requiring this flag, construct each member of the pair
        // individually.
        KeyType& key_dst = const_cast<KeyType&>(buckets[bucket_idx].get().first);
        ValueType& value_dst = buckets[bucket_idx].get().second;
        new(&key_dst) KeyType {kv_begin[kv_idx].first};
        new(&value_dst) ValueType {kv_begin[kv_idx].second};
#else
        new(&buckets[bucket_idx]) KeyValuePair(kv_begin[kv_idx]);
#endif
        total_inserts += 1;
      }
    });

  this->m_size += total_inserts.get();
  this->m_loadCount += total_inserts.get();
}

}  // namespace axom

#endif
