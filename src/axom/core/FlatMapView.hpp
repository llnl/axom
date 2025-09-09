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
 * \class FlatMapView
 *
 * \brief Provides a read-only view of a key-value container.
 *
 * \see FlatMap
 */
template <typename KeyType,
          typename ValueType,
          bool IsConst,
          typename Hash = detail::flat_map::HashMixer64<KeyType, DeviceHash>>
class FlatMapView : detail::flat_map::SequentialLookupPolicy<typename Hash::result_type>
{
private:
  using LookupPolicy = detail::flat_map::SequentialLookupPolicy<typename Hash::result_type>;
  using LookupPolicy::NO_MATCH;

  class IteratorImpl;
  friend class IteratorImpl;

  using BaseKVPair = std::pair<const KeyType, ValueType>;
  using KeyValuePair = std::conditional_t<IsConst, const BaseKVPair, BaseKVPair>;

public:
  using key_type = KeyType;
  using mapped_type = ValueType;
  using size_type = IndexType;
  using value_type = KeyValuePair;
  using iterator = IteratorImpl;
  using const_iterator = IteratorImpl;

  using FlatMapType =
    std::conditional_t<IsConst, const FlatMap<KeyType, ValueType, Hash>, FlatMap<KeyType, ValueType, Hash>>;

  AXOM_HOST_DEVICE FlatMapView() = default;

  FlatMapView(FlatMapType& other)
    : m_numGroups2(other.m_numGroups2)
    , m_size(other.m_size)
    , m_metadata(other.m_metadata.view())
    , m_buckets(other.m_buckets.view())
    , m_defaultValue()
  { }

  /*!
   * \brief Returns an iterator to the first valid object in the bucket array.
   */
  /// @{
  AXOM_HOST_DEVICE const_iterator begin() const
  {
    IndexType firstBucketIndex = this->nextValidIndex(m_metadata, NO_MATCH);
    return const_iterator(this, firstBucketIndex);
  }
  AXOM_HOST_DEVICE const_iterator cbegin() const
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
  AXOM_HOST_DEVICE const_iterator end() const { return const_iterator(*this, bucket_count()); }
  AXOM_HOST_DEVICE const_iterator cend() const { return const_iterator(*this, bucket_count()); }
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
  AXOM_HOST_DEVICE const_iterator find(const KeyType& key) const;
  /// @}

  /*!
   * \brief Find an entry with a given key.
   *
   *  If a corresponding value does not exist, a default value for the value
   *  type will be returned (but not inserted into the map).
   *
   * \param [in] key the key to search for
   *
   * \return The corresponding value, or a default value if the key does not exist
   *
   * \pre ValueType is default-constructible
   */
  /// @{
  AXOM_HOST_DEVICE const ValueType& operator[](const KeyType& key) const
  {
    static_assert(std::is_default_constructible<ValueType>::value,
                  "Cannot use axom::FlatMapView::operator[] when value type is not "
                  "default-constructible.");
    auto it = this->find(key);
    if(it != this->end())
    {
      return it->second;
    }
    return m_defaultValue;
  }
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
  AXOM_HOST_DEVICE IndexType bucket_count() const { return m_buckets.size(); }

private:
  /*!
   * \brief Returns whether another FlatMapView points to the same base FlatMap
   *  instance as this FlatMapView.
   *
   *  This is intended as a helper function for iterator comparisons.
   */
  AXOM_HOST_DEVICE bool isViewOfSameMap(const FlatMapView& other) const
  {
    return (m_metadata.data() == other.m_metadata.data() &&
            m_buckets.data() == other.m_buckets.data());
  }

  IndexType m_numGroups2 {-1};
  IndexType m_size {0};
  axom::ArrayView<const detail::flat_map::GroupBucket> m_metadata;

  // Storage details:
  using PairStorage = std::conditional_t<IsConst,
                                         const detail::flat_map::TypeErasedStorage<BaseKVPair>,
                                         detail::flat_map::TypeErasedStorage<BaseKVPair>>;
  axom::ArrayView<PairStorage> m_buckets;

  // Supporting functionality for FlatMapView::operator[]
  struct Dummy
  { };
  using DefaultValueType =
    std::conditional_t<std::is_default_constructible<ValueType>::value, ValueType, Dummy>;
  DefaultValueType m_defaultValue;
};

template <typename KeyType, typename ValueType, bool IsConst, typename Hash>
class FlatMapView<KeyType, ValueType, IsConst, Hash>::IteratorImpl
{
private:
  using MapType = FlatMapView<KeyType, ValueType, IsConst, Hash>;

  friend class FlatMapView<KeyType, ValueType, IsConst, Hash>;

public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename MapType::value_type;
  using difference_type = IndexType;

  using DataType = value_type;
  using pointer = DataType*;
  using reference = DataType&;

public:
  AXOM_HOST_DEVICE IteratorImpl() = default;

  AXOM_HOST_DEVICE IteratorImpl(const MapType& map, IndexType internalIdx)
    : m_map(map)
    , m_internalIdx(internalIdx)
  {
    assert(m_internalIdx >= 0 && m_internalIdx <= m_map.bucket_count());
  }

  AXOM_HOST_DEVICE friend bool operator==(const IteratorImpl& lhs, const IteratorImpl& rhs)
  {
    return (lhs.isViewOfSameMap(rhs) && lhs.m_internalIdx == rhs.m_internalIdx);
  }

  AXOM_HOST_DEVICE friend bool operator!=(const IteratorImpl& lhs, const IteratorImpl& rhs)
  {
    return !(lhs == rhs);
  }

  AXOM_HOST_DEVICE IteratorImpl& operator++()
  {
    m_internalIdx = m_map.nextValidIndex(m_map.m_metadata, m_internalIdx);
    return *this;
  }

  AXOM_HOST_DEVICE IteratorImpl operator++(int)
  {
    IteratorImpl next = *this;
    ++(*this);
    return next;
  }

  AXOM_HOST_DEVICE reference operator*() const { return m_map.m_buckets[m_internalIdx].get(); }

  AXOM_HOST_DEVICE pointer operator->() const { return &(m_map.m_buckets[m_internalIdx].get()); }

private:
  /*!
   * Wrapper function to call FlatMapView::isViewOfSameMap().
   *
   * This is a workaround for the following issue with friending:
   * https://www.open-std.org/jtc1/sc22/wg21/docs/cwg_active.html#1699
   */
  AXOM_HOST_DEVICE bool isViewOfSameMap(const IteratorImpl& other) const
  {
    return this->m_map.isViewOfSameMap(other.m_map);
  }

  MapType m_map;
  IndexType m_internalIdx {0};
};

template <typename KeyType, typename ValueType, bool IsConst, typename Hash>
AXOM_HOST_DEVICE auto FlatMapView<KeyType, ValueType, IsConst, Hash>::find(const KeyType& key) const
  -> const_iterator
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
auto FlatMap<KeyType, ValueType, Hash>::view() -> View
{
  return View(*this);
}

template <typename KeyType, typename ValueType, typename Hash>
auto FlatMap<KeyType, ValueType, Hash>::view() const -> ConstView
{
  return ConstView(*this);
}

}  // namespace axom

#endif
