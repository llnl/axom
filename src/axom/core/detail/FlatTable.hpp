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
  int getNext(int iter) const { return iter + 1; }
};

/*!
 * \brief A bit mixer used to ensure avalanching of a hash function.
 *
 *  Uses the "mxmxm" function described here:
 *  https://jonkagstrom.com/bit-mixer-construction/index.html
 */
template <typename KeyType, typename HashFunc>
struct HashMixer64
{
  uint64_t operator()(const KeyType& key) const
  {
    uint64_t hash = HashFunc {}(key);
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

  GroupBucket() : data {0ULL, 0ULL} { }

  int getEmptyBucket() const
  {
    for(int i = 0; i < Size; i++)
    {
      if(metadata.buckets[i] == GroupBucket::Empty)
      {
        return i;
      }
    }
    // Bucket not found.
    return InvalidSlot;
  }

  int nextFilledBucket(int start_index) const
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
  int visitHashBucket(std::uint8_t hash, Func&& visitor) const
  {
    std::uint8_t reducedHash = reduceHash(hash);
    for(int i = 0; i < Size; i++)
    {
      if(metadata.buckets[i] == reducedHash)
      {
        visitor(i);
      }
    }
    return InvalidSlot;
  }

  void setBucket(int index, std::uint8_t hash) { metadata.buckets[index] = reduceHash(hash); }

  void clearBucket(int index) { metadata.buckets[index] = Empty; }

  void setOverflow(std::uint8_t hash)
  {
    std::uint8_t hashOfwBit = 1 << (hash % 8);
    metadata.ofw |= hashOfwBit;
  }

  bool getMaybeOverflowed(std::uint8_t hash) const
  {
    std::uint8_t hashOfwBit = 1 << (hash % 8);
    return (metadata.ofw & hashOfwBit);
  }

  bool hasSentinel() const { return metadata.buckets[Size - 1] == Sentinel; }

  AXOM_HOST_DEVICE void setSentinel() { metadata.buckets[Size - 1] = Sentinel; }

  // We need to map hashes in the range [0, 255] to [2, 255], since 0 and 1
  // are taken by the "empty" and "sentinel" values respectively.
  static std::uint8_t reduceHash(std::uint8_t hash) { return (hash < 2) ? (hash + 8) : hash; }

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
  void probeIndex(int ngroups_pow_2,
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
      metadata[curr_group].visitHashBucket(hash_8, [&](IndexType bucket_index) {
        keep_going = on_hash_found(curr_group * GroupBucket::Size + bucket_index);
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

  IndexType nextValidIndex(ArrayView<const GroupBucket> metadata, int last_bucket) const
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

  const T& get() const { return *(reinterpret_cast<const T*>(&data)); }

  T& get() { return *(reinterpret_cast<T*>(&data)); }
};

}  // namespace flat_map
}  // namespace detail
}  // namespace axom

#endif  // Axom_Core_Detail_FlatTable_Hpp
