// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_Detail_FlatTable_Hpp
#define Axom_Core_Detail_FlatTable_Hpp

#include <climits>

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/utilities/BitUtilities.hpp"

#if defined(_MSC_VER)
  // MSVC does *not* define __SSE2__ to indicate SSE2 intrinsic support.
  #if defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
    // Either:
    //  - We're building for x86_64 (SSE2 is always available), or...
    //  - We're building for a 32-bit platform with SSE2 support.
    #define _AXOM_CORE_HAVE_SSE2
    #include <intrin.h>
  #endif
#elif defined(__SSE2__)
  // GCC, Clang, ICC all define the __SSE2__ macro.
  #define _AXOM_CORE_HAVE_SSE2
  #include <emmintrin.h>
#endif

namespace axom
{
namespace detail
{
namespace flat_map
{
struct QuadraticProbing
{
  /*!
   * \brief Returns the next offset to jump given an iteration count.
   *
   *  Each probe point for a given iteration is given by the formula
   *  H(i) = H_0 + i/2 + i^2/2 mod m. For m = 2^n, the sequence H(i) for
   *  i in [0, m-1) is a permutation on [0, m-1).
   *
   *  We return the offset from H(i) to H(i+1), which is i+1.
   */
  AXOM_HOST_DEVICE int getNext(int iter) const { return iter + 1; }
};

/*!
 * \brief A bit mixer used to ensure avalanching of a hash function.
 *
 *  Uses the "mxmxm" function described here:
 *  https://jonkagstrom.com/bit-mixer-construction/index.html
 */
template <typename KeyType, template <typename> class HashFunc>
struct HashMixer64
{
  using argument_type = typename HashFunc<KeyType>::argument_type;
  using result_type = typename HashFunc<KeyType>::result_type;

  AXOM_HOST_DEVICE uint64_t operator()(const KeyType& key) const
  {
    uint64_t hash = HashFunc<KeyType> {}(key);
    hash *= 0xbf58476d1ce4e5b9ULL;
    hash ^= hash >> 32;
    hash *= 0x94d049bb133111ebULL;
    hash ^= hash >> 32;
    hash *= 0x94d049bb133111ebULL;
    return hash;
  }
};

// We follow the design of boost::unordered_flat_map, which uses a 128-bit chunk
// of metadata for each group of 15 buckets.
// This is split up into an "overflow bit", and 15 bytes representing the
// state of each bucket.
//
// Each bucket byte has the value:
// - 0 if the bucket is "deleted"
// - 1 to signal the end of the bucket array.
// - if the bucket has an element, a reduced hash in the range [2, 255].
//
// The overflow bit acts as a Bloom filter to determine if probing should
// terminate.
struct GroupBucket
{
  constexpr static std::uint8_t Empty = 0;
  constexpr static std::uint8_t Sentinel = 1;
  constexpr static int InvalidSlot = -1;

  constexpr static int Size = 15;

  // We load the metadata into an SSE register as a 128-bit integer. Since
  // x86_64 is little endian, the expected order of bytes are reversed, and
  // the overflow byte is the LSB of the register.
  constexpr static int SSEBucketMask = ~(0x1);

  AXOM_HOST_DEVICE GroupBucket() : data {0ULL, 0ULL} { }

  AXOM_HOST_DEVICE int getEmptyBucket() const
  {
#if !defined(AXOM_DEVICE_CODE) && defined(_AXOM_CORE_HAVE_SSE2)
    auto match_simd = _mm_set1_epi8(Empty);
    auto metadata_simd = _mm_load_si128(reinterpret_cast<const __m128i*>(&metadata));
    match_simd = _mm_cmpeq_epi8(match_simd, metadata_simd);

    int mask_match = _mm_movemask_epi8(match_simd) & SSEBucketMask;
    // Little endian: count zeroes from right, bytes reversed from perspective
    // of SSE register
    int empty_bucket = axom::utilities::countr_zero(mask_match) - 1;

    return mask_match == 0 ? InvalidSlot : empty_bucket;
#else
    for(int i = 0; i < Size; i++)
    {
      if(metadata.buckets[i] == GroupBucket::Empty)
      {
        return i;
      }
    }
#endif
    // Bucket not found.
    return InvalidSlot;
  }

  AXOM_HOST_DEVICE int nextFilledBucket(int start_index) const
  {
    for(int i = start_index + 1; i < Size; i++)
    {
      // We intentionally don't check for the sentinel here. This gives us the
      // index of the sentinel bucket at the end, which is needed as a "stop"
      // for forward iteration.
      if(metadata.buckets[i] != GroupBucket::Empty)
      {
        return i;
      }
    }
    return InvalidSlot;
  }

  template <typename Func>
  AXOM_HOST_DEVICE int visitHashBucket(std::uint8_t hash, Func&& visitor) const
  {
    std::uint8_t reducedHash = reduceHash(hash);
#if !defined(AXOM_DEVICE_CODE) && defined(_AXOM_CORE_HAVE_SSE2)
    // Broadcast reduced hash across SIMD register
    auto hash_simd = _mm_set1_epi8(reducedHash);
    auto metadata_simd = _mm_load_si128(reinterpret_cast<const __m128i*>(&metadata));
    auto match_simd = _mm_cmpeq_epi8(hash_simd, metadata_simd);

    int mask_match = _mm_movemask_epi8(match_simd) & SSEBucketMask;
    while(mask_match != 0)
    {
      // Little endian: count zeroes from right, bytes reversed from perspective
      // of SSE register
      int bucket = axom::utilities::countr_zero(mask_match) - 1;
      if(!visitor(bucket))
      {
        // Found a match - break
        break;
      }
      mask_match ^= 1 << (bucket + 1);
    }

#else
    for(int i = 0; i < Size; i++)
    {
      if(metadata.buckets[i] == reducedHash)
      {
        visitor(i);
      }
    }
#endif
    return InvalidSlot;
  }

  /*!
   * \brief Visits matching hash buckets until an empty bucket is encountered.
   *
   *  This is used when performing batched insertion: since elements are only
   *  inserted, not deleted, an empty bucket will always be encountered only at
   *  the very end of a given probe sequence.
   *  The visitor function is used to allow for detecting duplicate keys.
   *
   * \param [in] hash reduced hash to search for
   * \param [in] visitor functor to call for each matching bucket slot
   *
   * \return the first empty slot found, or InvalidSlot
   */
  template <typename Func>
  AXOM_HOST_DEVICE int visitHashOrEmptyBucket(std::uint8_t hash, Func&& visitor) const
  {
    std::uint8_t reducedHash = reduceHash(hash);
    for(int i = 0; i < Size; i++)
    {
      if(metadata.buckets[i] == reducedHash)
      {
        visitor(i);
      }
      else if(metadata.buckets[i] == GroupBucket::Empty)
      {
        return i;
      }
    }
    return InvalidSlot;
  }

  template <bool Atomic = false>
  AXOM_HOST_DEVICE void setBucket(int index, std::uint8_t hash)
  {
    std::uint8_t reduced_hash = reduceHash(hash);
    if(Atomic)  // TODO: should be constexpr
    {
#if defined(AXOM_USE_CUDA)
      // CUDA workaround for lack of 8-bit atomicStore.
      volatile std::uint8_t* bucket = &(metadata.buckets[index]);
      *bucket = reduced_hash;
#else
      axom::atomicStore<axom::auto_atomic>(&(metadata.buckets[index]), reduced_hash);
#endif
      return;
    }
    else
    {
      metadata.buckets[index] = reduced_hash;
    }
  }

  AXOM_HOST_DEVICE void clearBucket(int index) { metadata.buckets[index] = Empty; }

  template <bool Atomic = false>
  AXOM_HOST_DEVICE void setOverflow(std::uint8_t hash)
  {
    std::uint8_t hashOfwBit = 1 << (hash % 8);
    if(Atomic)  // TODO: should be constexpr
    {
#if defined(AXOM_USE_HIP) || defined(AXOM_USE_CUDA)
      // Workaround for a lack of an atomicOr builtin for uint8_t.
      GroupBucket pack_data {};
      pack_data.metadata.ofw = hashOfwBit;
      axom::atomicOr<axom::auto_atomic>(&data[0], pack_data.data[0]);
#else
      axom::atomicOr<axom::auto_atomic>(&(metadata.ofw), hashOfwBit);
#endif
    }
    else
    {
      metadata.ofw |= hashOfwBit;
    }
  }

  template <bool Atomic = false>
  AXOM_HOST_DEVICE bool getMaybeOverflowed(std::uint8_t hash) const
  {
    std::uint8_t hashOfwBit = 1 << (hash % 8);
    std::uint8_t curr_ofw;
    if(Atomic)  // TODO: should be constexpr
    {
#if defined(AXOM_USE_CUDA)
      // CUDA workaround for lack of 8-bit atomicLoad.
      curr_ofw = *const_cast<const volatile std::uint8_t*>(&(metadata.ofw));
#else
      // TODO: why is the const_cast required? (RAJA issue?)
      curr_ofw = axom::atomicLoad<axom::auto_atomic>(const_cast<std::uint8_t*>(&(metadata.ofw)));
#endif
    }
    else
    {
      curr_ofw = metadata.ofw;
    }
    return (curr_ofw & hashOfwBit);
  }

  bool hasSentinel() const { return metadata.buckets[Size - 1] == Sentinel; }

  AXOM_HOST_DEVICE void setSentinel() { metadata.buckets[Size - 1] = Sentinel; }

  // We need to map hashes in the range [0, 255] to [2, 255], since 0 and 1
  // are taken by the "empty" and "sentinel" values respectively.
  AXOM_HOST_DEVICE static std::uint8_t reduceHash(std::uint8_t hash)
  {
    return (hash < 2) ? (hash + 8) : hash;
  }

  union alignas(16)
  {
    struct
    {
      std::uint8_t ofw;
      std::uint8_t buckets[Size];
    } metadata;
    std::uint64_t data[2];
    static_assert(sizeof(metadata) == sizeof(data),
                  "flat_map::GroupBucket: sizeof(data_bytes) != sizeof(data)");
  };
};

static_assert(sizeof(GroupBucket) == 16, "flat_map::GroupBucket: size != 16 bytes");
static_assert(std::alignment_of<GroupBucket>::value == 16,
              "flat_map::GroupBucket: alignment != 16 bytes");
static_assert(std::is_standard_layout<GroupBucket>::value,
              "flat_map::GroupBucket: not standard layout");

template <typename HashType, typename ProbePolicy = QuadraticProbing>
struct SequentialLookupPolicy : ProbePolicy
{
  constexpr static int NO_MATCH = -1;

  /*!
   * \brief Inserts a hash into the first empty bucket in an array of groups
   *  for an open-addressing hash map.
   *
   * \param [in] ngroups_pow_2 the number of groups, expressed as a power of 2
   * \param [in] metadata the array of metadata for the groups in the hash map
   * \param [in] hash the hash to insert
   */
  IndexType probeEmptyIndex(int ngroups_pow_2, ArrayView<GroupBucket> metadata, HashType hash) const
  {
    // We use the k MSBs of the hash as the initial group probe point,
    // where ngroups = 2^k.
    int bitshift_right = ((CHAR_BIT * sizeof(HashType)) - ngroups_pow_2);
    HashType curr_group = hash >> bitshift_right;
    curr_group &= ((1 << ngroups_pow_2) - 1);
    int empty_group = NO_MATCH;
    int empty_bucket = NO_MATCH;

    std::uint8_t hash_8 = static_cast<std::uint8_t>(hash);
    for(int iteration = 0; iteration < metadata.size(); iteration++)
    {
      int tentative_empty_bucket = metadata[curr_group].getEmptyBucket();
      if(tentative_empty_bucket != GroupBucket::InvalidSlot && empty_group == NO_MATCH)
      {
        empty_group = curr_group;
        empty_bucket = tentative_empty_bucket;
      }

      if((!metadata[curr_group].getMaybeOverflowed(hash_8) && empty_group != NO_MATCH))
      {
        // We've reached the last group that might contain the hash.
        // Stop probing.
        break;
      }
      else if(empty_group == NO_MATCH)
      {
        // Set the overflow bit and continue probing.
        metadata[curr_group].setOverflow(hash_8);
      }
      curr_group = (curr_group + this->getNext(iteration)) % metadata.size();
    }
    if(empty_group != NO_MATCH)
    {
      return empty_group * GroupBucket::Size + empty_bucket;
    }
    return NO_MATCH;
  }

  /*!
   * \brief Finds the next potential bucket index for a given hash in a group
   *  array for an open-addressing hash map.
   *
   * \param [in] ngroups_pow_2 the number of groups, expressed as a power of 2
   * \param [in] metadata the array of metadata for the groups in the hash map
   * \param [in] hash the hash to insert
   * \param [in] on_hash_found functor to call for a bucket index with a
   *  matching hash
   */
  template <typename FoundIndex>
  AXOM_HOST_DEVICE void probeIndex(int ngroups_pow_2,
                                   ArrayView<const GroupBucket> metadata,
                                   HashType hash,
                                   FoundIndex&& on_hash_found) const
  {
    // We use the k MSBs of the hash as the initial group probe point,
    // where ngroups = 2^k.
    int bitshift_right = ((CHAR_BIT * sizeof(HashType)) - ngroups_pow_2);
    HashType curr_group = hash >> bitshift_right;
    curr_group &= ((1 << ngroups_pow_2) - 1);

    std::uint8_t hash_8 = static_cast<std::uint8_t>(hash);
    bool keep_going = true;
    for(int iteration = 0; iteration < metadata.size(); iteration++)
    {
      metadata[curr_group].visitHashBucket(hash_8, [&](IndexType bucket_index) -> bool {
        keep_going = on_hash_found(curr_group * GroupBucket::Size + bucket_index);
        return keep_going;
      });

      if(!metadata[curr_group].getMaybeOverflowed(hash_8))
      {
        // Stop probing if the "overflow" bit is not set.
        keep_going = false;
      }

      if(!keep_going)
      {
        break;
      }
      // Probe the next bucket.
      curr_group = (curr_group + this->getNext(iteration)) % metadata.size();
    }
  }

  void setBucketHash(ArrayView<GroupBucket> metadata, IndexType bucket, HashType hash)
  {
    int group_index = bucket / GroupBucket::Size;
    int slot_index = bucket % GroupBucket::Size;

    metadata[group_index].setBucket(slot_index, hash);
  }

  bool clearBucket(ArrayView<GroupBucket> metadata, IndexType bucket, HashType hash)
  {
    int group_index = bucket / GroupBucket::Size;
    int slot_index = bucket % GroupBucket::Size;

    metadata[group_index].clearBucket(slot_index);

    // Return if the overflow bit is set on the bucket. That indicates whether
    // we are deleting an element in the middle of a probing sequence.
    return metadata[group_index].getMaybeOverflowed(hash);
  }

  AXOM_HOST_DEVICE IndexType nextValidIndex(ArrayView<const GroupBucket> metadata, int last_bucket) const
  {
    if(last_bucket >= metadata.size() * GroupBucket::Size - 1)
    {
      return last_bucket;
    }
    int group_index = last_bucket / GroupBucket::Size;
    int slot_index = last_bucket % GroupBucket::Size;

    do
    {
      slot_index = metadata[group_index].nextFilledBucket(slot_index);
      if(slot_index == GroupBucket::InvalidSlot)
      {
        group_index++;
        slot_index = -1;
      }
    } while(slot_index == GroupBucket::InvalidSlot && group_index < metadata.size());

    return group_index * GroupBucket::Size + slot_index;
  }
};

template <typename T>
struct alignas(T) TypeErasedStorage
{
  unsigned char data[sizeof(T)];

  AXOM_HOST_DEVICE const T& get() const { return *(reinterpret_cast<const T*>(&data)); }

  AXOM_HOST_DEVICE T& get() { return *(reinterpret_cast<T*>(&data)); }
};

}  // namespace flat_map
}  // namespace detail
}  // namespace axom

#undef _AXOM_CORE_HAVE_SSE2

#endif  // Axom_Core_Detail_FlatTable_Hpp
