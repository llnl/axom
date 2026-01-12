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

#if defined(AXOM_USE_CUDA)
template <typename KeyType, typename ValueType>
AXOM_DEVICE void constructPairInPlace(std::pair<const KeyType, ValueType>& pair,
                                      KeyType key,
                                      ValueType value)
{
  // HACK: std::pair constructor is not host-device annotated, but CUDA
  // requires passing in --expt-relaxed-constexpr for it to work.
  // Instead of requiring this flag, construct each member of the pair
  // individually.
  KeyType& key_dst = const_cast<KeyType&>(pair.first);
  ValueType& value_dst = pair.second;
  new(&key_dst) KeyType {key};
  new(&value_dst) ValueType {value};
}
#endif

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
  AXOM_HOST_DEVICE value_type operator*() const
  {
#if defined(__CUDA_ARCH__)
    value_type kv_pair;
    constructPairInPlace(kv_pair, *m_keyIter, *m_valueIter);
    return kv_pair;
#else
    return {*m_keyIter, *m_valueIter};
#endif
  }

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

/*!
 * \class FlatMapOffsetIterator
 * \brief Iterator helper for iterating over filled buckets given an array of
 *  bucket indices.
 */
template <typename KeyType, typename ValueType>
class FlatMapOffsetIterator
  : public IteratorBase<FlatMapOffsetIterator<KeyType, ValueType>, IndexType>
{
private:
  using BaseType = IteratorBase<FlatMapOffsetIterator<KeyType, ValueType>, IndexType>;
  using KeyValuePair = std::pair<const KeyType, ValueType>;
  using PairStorage = detail::flat_map::TypeErasedStorage<KeyValuePair>;

public:
  // Iterator traits required to satisfy LegacyRandomAccessIterator concept
  // before C++20
  // See: https://en.cppreference.com/w/cpp/iterator/iterator_traits
  using difference_type = IndexType;
  using value_type = std::pair<const KeyType, ValueType>;
  using reference = value_type&;
  using pointer = value_type*;
  using iterator_category = std::random_access_iterator_tag;

  FlatMapOffsetIterator() = default;

  AXOM_HOST_DEVICE FlatMapOffsetIterator(PairStorage* buckets, IndexType* offsets)
    : m_buckets {buckets}
    , m_offsets {offsets}
  { }

  AXOM_HOST_DEVICE reference operator*() const { return m_buckets[m_offsets[this->m_pos]].get(); }

protected:
  /** Implementation of advance() as required by IteratorBase */
  AXOM_HOST_DEVICE void advance(IndexType n) { BaseType::m_pos += n; }

private:
  PairStorage* m_buckets {nullptr};
  IndexType* m_offsets {nullptr};
};

/**
 * \brief Helper function to gather filled buckets within a FlatMap.
 *
 *  Workaround for a limitation within CUDA where a lambda cannot be defined in
 *  a protected or private member function.
 */
template <typename ExecSpace>
void gatherFilledBuckets(ArrayView<flat_map::GroupBucket> group_metadata,
                         ArrayView<IndexType> filled_bucket_indexes,
                         IndexType num_buckets,
                         int allocator_id)
{
  using flat_map::GroupBucket;

  IndexType* counter = axom::allocate<IndexType>(1, allocator_id);
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA)
  if(detail::getAllocatorSpace(allocator_id) == MemorySpace::Device)
  {
    for_all<ExecSpace>(1, AXOM_LAMBDA(IndexType) { *counter = 0; });
  }
  else
  {
    *counter = 0;
  }
#else
  *counter = 0;
#endif

  for_all<ExecSpace>(
    num_buckets,
    AXOM_LAMBDA(IndexType bucket_idx) {
      IndexType group_idx = bucket_idx / GroupBucket::Size;
      int slot_idx = bucket_idx % GroupBucket::Size;
      if(group_metadata[group_idx].metadata.buckets[slot_idx] > GroupBucket::Sentinel)
      {
        // Bucket contains an element.
        IndexType dest = axom::atomicAdd<ExecSpace>(counter, IndexType {1});
        filled_bucket_indexes[dest] = bucket_idx;
      }
    });

  axom::deallocate(counter);
}

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
template <typename ExecSpace>
void FlatMap<KeyType, ValueType, Hash>::parallelRehash(IndexType count)
{
  using detail::FlatMapOffsetIterator;

  // If the FlatMap is constructed in device-only memory, construct in
  // parallel on the device.
  axom::Array<IndexType> filled_bucket_idx_vec(m_size, m_size, m_allocator.getID());

  detail::gatherFilledBuckets<ExecSpace>(m_metadata.view(),
                                         filled_bucket_idx_vec.view(),
                                         m_buckets.size(),
                                         m_allocator.getID());

  FlatMapOffsetIterator bucket_begin {m_buckets.data(), filled_bucket_idx_vec.data()};

  FlatMap new_map(count, m_allocator);
  new_map.template insert<ExecSpace>(bucket_begin, bucket_begin + m_size);

  this->swap(new_map);
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

  IndexType num_elems = std::distance(kv_begin, kv_end);

  // Batched insertion assumes probing sequences are gap-free
  // (i.e., there are no tombstones from prior erase() operations).
  // When tombstones exist, the parallel insertion logic can mishandle duplicates
  // under contention (e.g. OpenMP) and produce incorrect size/value results.
  //
  // If tombstones exist, rehash to compact the table and restore the invariants required by this algorithm.
  if(this->m_loadCount != static_cast<std::uint64_t>(this->m_size))
  {
    this->rehash(this->m_size + num_elems);
  }

  const bool is_gap_free = (this->m_loadCount == static_cast<std::uint64_t>(this->m_size));

  // Assume that all elements will be inserted into an empty slot.
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

  axom::ReduceSum<ExecSpace, IndexType> total_overwrites(0);

  for_all<ExecSpace>(
    num_elems,
    AXOM_LAMBDA(IndexType idx) {
      // Construct key.
      KeyType key = (*(kv_begin + idx)).first;

      // Hash keys.
      auto hash = Hash {}(key);

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
          meta_group[curr_group].visitHashBucket(hash_8, [&](int matching_slot) -> bool {
            IndexType bucket_index = curr_group * GroupBucket::Size + matching_slot;

            if(buckets[bucket_index].get().first == key)
            {
              duplicate_bucket_index = bucket_index;
              return false;  // Don't need to search other buckets.
            }
            return true;
          });
          int empty_slot_index = meta_group[curr_group].getEmptyBucket();

          if(duplicate_bucket_index == -1 && empty_bucket_index == -1)
          {
            // Default probing behavior: no duplicate found yet, and no empty
            // bucket found prior.
            if(empty_slot_index == GroupBucket::InvalidSlot)
            {
              // Group is full. Set overflow bit for the group.
              meta_group[curr_group].template setOverflow<true>(hash_8);
            }
            else
            {
              // Update empty bucket index with first empty slot we encounter.
              empty_bucket_index = curr_group * GroupBucket::Size + empty_slot_index;
              key_index_dedup[empty_bucket_index] = idx;
              key_index_to_bucket[idx] = empty_bucket_index;

              // Insert initial element, this will be updated with the value of
              // the "winning" key-value pair.
              meta_group[curr_group].template setBucket<true>(empty_slot_index, hash_8);
#if defined(__CUDA_ARCH__)
              detail::constructPairInPlace(buckets[empty_bucket_index].get(),
                                           key,
                                           (*(kv_begin + idx)).second);
#else
              new(&buckets[empty_bucket_index]) KeyValuePair(*(kv_begin + idx));
#endif
            }
          }
          else if(duplicate_bucket_index != -1)
          {
            // Found a duplicate bucket.
            if(!is_gap_free && empty_bucket_index != -1)
            {
              // We've already encountered an empty bucket earlier to place a
              // k-v pair. This may occur if a probing sequence contains gaps
              // (insertions followed by erasures).
              //
              // Just erase this element.
              total_overwrites += 1;

              int slot_index = duplicate_bucket_index - curr_group * GroupBucket::Size;
              buckets[duplicate_bucket_index].get().~KeyValuePair();
              meta_group[curr_group].clearBucket(slot_index);
            }
            else
            {
              if(key_index_dedup[duplicate_bucket_index] == -1)
              {
                // The k-v pair matches an already-existing pair in the map.
                // Keep track of the number of overwrites so that we don't
                // double-count them when incrementing the size.
                total_overwrites += 1;
              }
              // Highest-indexed kv pair wins.
              axom::atomicMax<ExecSpace>(&key_index_dedup[duplicate_bucket_index], idx);
              key_index_to_bucket[idx] = duplicate_bucket_index;
            }
          }
          // Unlock group once we're done.
          group_locks[curr_group].unlock();

          if(duplicate_bucket_index != -1)
          {
            // We've found a duplicate key to overwrite.
            break;
          }
          else if(empty_bucket_index != -1 &&
                  (is_gap_free || !meta_group[curr_group].getMaybeOverflowed(hash_8)))
          {
            // If we're inserting into a gap-free map, empty bucket signals the
            // end of the probing sequence.
            // Otherwise, we need to check the overflow mask to continue probing.
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
        detail::constructPairInPlace(buckets[bucket_idx].get(),
                                     (*(kv_begin + kv_idx)).first,
                                     (*(kv_begin + kv_idx)).second);
#else
        new(&buckets[bucket_idx]) KeyValuePair(*(kv_begin + kv_idx));
#endif
        total_inserts += 1;
      }
    });

  this->m_size += total_inserts.get() - total_overwrites.get();
  this->m_loadCount += total_inserts.get() - total_overwrites.get();
}

}  // namespace axom

#endif
