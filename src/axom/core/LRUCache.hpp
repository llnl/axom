// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_LRUCACHE_HPP_
#define AXOM_CORE_LRUCACHE_HPP_

#include <cassert>
#include <cstddef>
#include <list>
#include <unordered_map>
#include <utility>

namespace axom
{

/*!
 * \brief Fixed-capacity least-recently-used cache.
 *
 * \tparam Key Cache key type. The default implementation uses
 *  `std::unordered_map`, so `Key` must be hashable and equality comparable.
 * \tparam Value Cached value type.
 *
 * `LRUCache` stores at most `capacity()` entries. When insertion of a new key
 * would exceed the configured capacity, the least recently used entry is
 * evicted. Successful calls to `find()` and `insert()` promote the affected
 * entry to most-recently-used position.
 *
 * Lookup, insertion, promotion, and eviction are average constant time. Pointers
 * and references to cached values remain valid until the corresponding entry is
 * replaced, erased, or evicted. The class is not internally synchronized;
 * callers that share a cache across threads must provide external locking.
 */
template <typename Key, typename Value>
class LRUCache
{
public:
  /*!
   * \brief Construct a cache with a fixed positive capacity.
   *
   * \param [in] capacity Maximum number of entries retained by this cache.
   *
   * \pre `capacity > 0`.
   */
  explicit LRUCache(std::size_t capacity) : m_capacity(capacity) { assert(m_capacity > 0); }

  /*!
   * \brief Find an entry by key and promote it to most-recently-used position.
   *
   * \param [in] key Key to find.
   * \return Pointer to the cached value, or `nullptr` when no entry exists.
   */
  Value* find(const Key& key)
  {
    auto it = m_entries.find(key);
    if(it == m_entries.end())
    {
      return nullptr;
    }

    touch(it);
    return &it->second.value;
  }

  /*!
   * \brief Find an entry by key without changing its recency.
   *
   * \param [in] key Key to find.
   * \return Pointer to the cached value, or `nullptr` when no entry exists.
   */
  const Value* peek(const Key& key) const
  {
    auto it = m_entries.find(key);
    return it == m_entries.end() ? nullptr : &it->second.value;
  }

  /*!
   * \brief Insert or replace an entry and promote it to most-recently-used position.
   *
   * If `key` already exists, this replaces the value and promotes the entry.
   * If `key` is new and the cache is full, the least recently used entry is
   * evicted before insertion.
   *
   * \param [in] key Key for the cache entry.
   * \param [in] value Value to store.
   * \return Reference to the stored value.
   */
  Value& insert(Key key, Value value)
  {
    auto it = m_entries.find(key);
    if(it != m_entries.end())
    {
      it->second.value = std::move(value);
      touch(it);
      return it->second.value;
    }

    if(m_entries.size() >= m_capacity)
    {
      m_entries.erase(m_lru_order.back());
      m_lru_order.pop_back();
    }

    m_lru_order.push_front(key);
    auto inserted = m_entries.emplace(std::move(key), Entry {std::move(value), m_lru_order.begin()});
    return inserted.first->second.value;
  }

  /*!
   * \brief Check whether a key is present without changing recency.
   */
  bool contains(const Key& key) const { return m_entries.find(key) != m_entries.end(); }

  /*!
   * \brief Remove all entries from the cache.
   */
  void clear()
  {
    m_entries.clear();
    m_lru_order.clear();
  }

  /*!
   * \brief Return the number of entries currently stored.
   */
  std::size_t size() const { return m_entries.size(); }

  /*!
   * \brief Return the maximum number of entries retained.
   */
  std::size_t capacity() const { return m_capacity; }

  /*!
   * \brief Return true if the cache stores no entries.
   */
  bool empty() const { return m_entries.empty(); }

private:
  struct Entry
  {
    Value value;
    typename std::list<Key>::iterator lru_position;
  };

  using EntryIterator = typename std::unordered_map<Key, Entry>::iterator;

  void touch(EntryIterator it)
  {
    m_lru_order.splice(m_lru_order.begin(), m_lru_order, it->second.lru_position);
    it->second.lru_position = m_lru_order.begin();
  }

  std::size_t m_capacity;
  std::unordered_map<Key, Entry> m_entries;
  std::list<Key> m_lru_order;
};

}  // namespace axom

#endif  // AXOM_CORE_LRUCACHE_HPP_
