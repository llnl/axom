// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_FlatMap_View_HPP
#define Axom_Core_FlatMap_View_HPP

#include "axom/core/FlatMap.hpp"

namespace axom
{

/*!
 * \class FlatMap
 *
 * \brief Provides a read-only view of a key-value container.
 *
 * \see FlatMapView
 */
template <typename KeyType, typename ValueType, typename Hash = detail::flat_map::HashMixer64<KeyType, DeviceHash>>
class FlatMapView : detail::flat_map::SequentialLookupPolicy<typename Hash::result_type>
{
private:
  using LookupPolicy = detail::flat_map::SequentialLookupPolicy<typename Hash::result_type>;
  using LookupPolicy::NO_MATCH;

  class IteratorImpl;
  friend class IteratorImpl;

  using KeyValuePair = std::pair<const KeyType, ValueType>;

public:
  using key_type = KeyType;
  using mapped_type = ValueType;
  using size_type = IndexType;
  using value_type = KeyValuePair;
  using iterator = IteratorImpl;
  using const_iterator = IteratorImpl;

  FlatMapView(const FlatMap<KeyType, ValueType, Hash>& other)
    : m_numGroups2(other.m_numGroups2)
    , m_size(other.m_size)
    , m_metadata(other.m_metadata.view())
    , m_buckets(other.m_buckets.view())
  { }

  /*!
   * \brief Returns an iterator to the first valid object in the bucket array.
   */
  /// @{
  const_iterator begin() const
  {
    IndexType firstBucketIndex = this->nextValidIndex(m_metadata, NO_MATCH);
    return const_iterator(this, firstBucketIndex);
  }
  const_iterator cbegin() const
  {
    IndexType firstBucketIndex = this->nextValidIndex(m_metadata, NO_MATCH);
    return const_iterator(this, firstBucketIndex);
  }
  /// @}

  /*!
   * \brief Returns an iterator to "one past" the last valid object in the
   *  bucket array.
   */
  /// @{
  const_iterator end() const { return const_iterator(*this, bucket_count()); }
  const_iterator cend() const { return const_iterator(*this, bucket_count()); }
  /// @}

  /*!
   * \brief Try to find an entry with a given key.
   *
   * \param [in] key the key to search for
   *
   * \return An iterator pointing to the corresponding key-value pair, or end()
   *  if the key wasn't found.
   */
  /// @{
  const_iterator find(const KeyType& key) const;
  /// @}

  /*!
   * \brief Returns true if there are no entries in the FlatMap, false
   *  otherwise.
   */
  bool empty() const { return m_size == 0; }

  /*!
   * \brief Returns the number of entries stored in the FlatMap.
   */
  IndexType size() const { return m_size; }

  /*!
   * \brief Return the number of entries matching a given key.
   *
   *  This method will always return 0 or 1.
   *
   * \param [in] key the key to search for
   */
  IndexType count(const KeyType& key) const
  {
    return contains(key) ? IndexType {1} : IndexType {0};
  }

  /*!
   * \brief Return true if the FlatMap contains a key, false otherwise.
   *
   * \param [in] key the key to search for
   */
  bool contains(const KeyType& key) const { return (find(key) != end()); }

  /*!
   * \brief Returns the number of buckets allocated in the FlatMap.
   *
   *  The maximum number of elements that can be stored in the FlatMap without
   *  resizing and rehashing is bucket_count() * max_load_factor().
   */
  IndexType bucket_count() const { return m_buckets.size(); }

private:
  IndexType m_numGroups2;
  IndexType m_size;
  axom::ArrayView<const detail::flat_map::GroupBucket> m_metadata;

  // Storage details:
  using PairStorage = detail::flat_map::TypeErasedStorage<KeyValuePair>;
  axom::ArrayView<const PairStorage> m_buckets;
};

template <typename KeyType, typename ValueType, typename Hash>
class FlatMapView<KeyType, ValueType, Hash>::IteratorImpl
{
private:
  using MapType = FlatMapView<KeyType, ValueType, Hash>;

  friend class FlatMapView<KeyType, ValueType, Hash>;

public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename MapType::value_type;
  using difference_type = IndexType;

  using DataType = const value_type;
  using pointer = DataType*;
  using reference = DataType&;

public:
  IteratorImpl(const MapType& map, IndexType internalIdx) : m_map(map), m_internalIdx(internalIdx)
  {
    assert(m_internalIdx >= 0 && m_internalIdx <= m_map.bucket_count());
  }

  friend bool operator==(const IteratorImpl& lhs, const IteratorImpl& rhs)
  {
    return lhs.m_internalIdx == rhs.m_internalIdx;
  }

  friend bool operator!=(const IteratorImpl& lhs, const IteratorImpl& rhs)
  {
    return lhs.m_internalIdx != rhs.m_internalIdx;
  }

  IteratorImpl& operator++()
  {
    m_internalIdx = m_map.nextValidIndex(m_map.m_metadata, m_internalIdx);
    return *this;
  }

  IteratorImpl operator++(int)
  {
    IteratorImpl next = *this;
    ++(*this);
    return next;
  }

  reference operator*() const { return m_map.m_buckets[m_internalIdx].get(); }

  pointer operator->() const { return &(m_map.m_buckets[m_internalIdx].get()); }

private:
  MapType m_map;
  IndexType m_internalIdx;
};

template <typename KeyType, typename ValueType, typename Hash>
auto FlatMapView<KeyType, ValueType, Hash>::find(const KeyType& key) const -> const_iterator
{
  auto hash = Hash {}(key);
  iterator found_iter = end();
  this->probeIndex(m_numGroups2, m_metadata, hash, [&](IndexType bucket_index) -> bool {
    if(this->m_buckets[bucket_index].get().first == key)
    {
      found_iter = iterator(*this, bucket_index);
      // Stop tracking.
      return false;
    }
    return true;
  });
  return found_iter;
}

template <typename KeyType, typename ValueType, typename Hash>
auto FlatMap<KeyType, ValueType, Hash>::view() const -> View
{
  return View(*this);
}

}  // namespace axom

#endif
