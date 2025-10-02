// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MAP_HPP_
#define AXOM_MAP_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/atomics.hpp"
#include "axom/core/Types.hpp"

// C/C++ includes
#include <functional>
#include <iostream>

namespace axom
{
namespace experimental
{
namespace axom_map
{
/// \name Map Supporting Data Structures
/// @{
/*!
   * \class Node
   *
   * \brief Node is the foundational data type of the Map, being the lowest-level data storage unit,
   *  with the key and value data for a given item, along with a next value for the linked list.
   *
   * \tparam Key the type of key for the key-value pair stored within the Node instance.
   * \tparam T the type of value for the key-value pair stored within the Node instance.
   */
template <typename Key, typename T>
struct Node
{
  IndexType next;
  Key key;
  T value;

  friend bool operator==(const Node& lhs, const Node& rhs)
  {
    return (lhs.next == rhs.next) && (lhs.key == rhs.key);
  }
};

/*!
   * \class Pair
   *
   * \brief Pair is strictly for returning from insert, to match the pair-returning setup of STL
   *  unordered_map.
   *
   * \tparam Key the type of key stored in the associated Node.
   * \tparam T the type of value stored in the associated Node.
   */
template <typename Key, typename T>
struct Pair
{
  axom_map::Node<Key, T>* first;
  bool second;

  Pair(axom_map::Node<Key, T>* node, bool status) : first(node), second(status) { }
};

/*!
   * \class Bucket
   *
   * \brief Bucket is the implementation of a bucket for this Map with chaining. It's a singly-linked
   * list, and supports search, insertion, and removal.
   *
   * \tparam Key the type of the keys stored in each Node in the list.
   * \tparapm T the type of the values stored in each Node in the list.
   */
template <typename Key, typename T>
class Bucket
{
public:
  /// \name Bucket Construction
  ///@{
  /*!
     * \brief Default constructor so new Buckets just construct into nothing.
     *  Can likely just be removed.
     */
  Bucket() { }

  /*!
     * \brief Basic constructor, creates singly-linked list with provided size.
     *
     * \param [in] len the number of items this Bucket instance needs to be able to store.
     *
     * \pre len > 0
     */
  Bucket(IndexType len) { init(len); }

  /*!
     * \brief If a Bucket instance is created with the default contructor, this
     *  is used to initialize and allocate all the member variables.
     *
     * \param [in] len the number of items this Bucket instance needs to be able to store.
     *
     * \pre len > 0
     */
  void init(IndexType len)
  {
    m_list = axom::allocate<axom_map::Node<Key, T>>(len);
    for(int i = 0; i < len - 1; i++)
    {
      m_list[i].next = i + 1;
    }
    m_list[len - 1].next = -1;
    m_free = 0;
    m_head = -1;
    m_capacity = len;
    m_size = 0;
    m_end.key = Key {0};
    m_end.value = T {0};
    m_end.next = -2;
  }
  /// @}
  /// \name Bucket Methods
  /// @{
  /*!
     * \brief Inserts a given key-value pair into the linked list if an item with
     *  given key does not already exist. Fails if list is full or item with given key
     *  already exists.
     *
     * \param [in] key the key of the key-value pair to be inserted.
     * \param [in] value the value of the key-value pair to be inserted.
     *
     * \return Returns a Pair object, with the first element being a pointer to the inserted node
     *  and the second element being boolean value true if successful. If insertion fails,
     *  the bool is set to False, and the node pointed to is the end node if the list was full,
     *  or the item occupying the queried slot if an item with the given key already existed.
     */
  axom_map::Pair<Key, T> insert_no_update(const Key& key, const T& value)
  {
    if(m_free != -1)
    {
      IndexType ind = m_head;
      if(m_head == -1)
      {
        m_head = m_free;
        ind = m_free;
        m_free = m_list[m_free].next;
        m_list[m_head].next = -1;
        m_list[m_head].key = key;
        m_list[m_head].value = value;
        m_size++;
      }
      else
      {
        while(m_list[ind].next != -1)
        {
          if(m_list[ind].key == key)
          {
            axom_map::Pair<Key, T> ret(&(m_list[ind]), false);
            return ret;
          }
          ind = m_list[ind].next;
        }
        m_list[ind].next = m_free;
        ind = m_free;
        m_free = m_list[m_free].next;
        m_list[ind].next = -1;
        m_list[ind].key = key;
        m_list[ind].value = value;
        m_size++;
      }
      axom_map::Pair<Key, T> ret(&(m_list[ind]), true);
      return ret;
    }
    axom_map::Pair<Key, T> ret(&(m_end), false);
    return ret;
  }

  /*!
     * \brief Inserts a given key-value pair into the linked list if an item with
     *  given key does not already exist, and updates the item's current value with the supplied
     *  one if item already exists. Fails if list is full.
     *
     * \param [in] key the key of the key-value pair to be inserted.
     * \param [in] value the value of the key-value pair to be inserted.
     *
     * \return Returns a Pair object, with the first element being a pointer to the inserted node
     *  and the second element being boolean value true if successful. If insertion fails, or assignment
     *  was performed, the bool is set to False. If assignment occured, the node pointed to is the item
     *  corresponding to the supplied key. Otherwise, if insertion and assignment failed, the node pointed to is the end node.
     */
  axom_map::Pair<Key, T> insert_update(const Key& key, const T& value)
  {
    IndexType ind = m_head;
    if(m_head == -1)
    {
      m_head = m_free;
      ind = m_free;
      m_free = m_list[m_free].next;
      m_list[m_head].next = -1;
      m_list[m_head].key = key;
      m_list[m_head].value = value;
      m_size++;
      axom_map::Pair<Key, T> ret(&(m_list[m_head]), true);
      return ret;
    }
    else
    {
      while(m_list[ind].next != -1)
      {
        if(m_list[ind].key == key)
        {
          m_list[ind].value = value;
          axom_map::Pair<Key, T> ret(&(m_list[ind]), false);
          return ret;
        }
        ind = m_list[ind].next;
      }

      if(m_list[ind].key == key)
      {
        m_list[ind].value = value;
        axom_map::Pair<Key, T> ret(&(m_list[ind]), false);
        return ret;
      }
      else if(m_free != -1)
      {
        m_list[ind].next = m_free;
        ind = m_free;
        m_free = m_list[m_free].next;
        m_list[ind].next = -1;
        m_list[ind].key = key;
        m_list[ind].value = value;
        m_size++;

        axom_map::Pair<Key, T> ret(&(m_list[ind]), true);
        return ret;
      }
    }
    axom_map::Pair<Key, T> ret(&(m_end), false);
    return ret;
  }

  /*!
     * \brief Removes item with the given key from the linked list, if item exists. Relinks accordingly.
     *
     * \param [in] key the key of the item to be found and removed.
     *
     * \return Returns true if item was found and removed, false otherwise.
     *
     */
  bool remove(const Key& key)
  {
    if(m_head == -1)
    {
      return false;
    }
    IndexType ind = m_head;
    IndexType prev = -1;
    do
    {
      if(m_list[ind].key == key)
      {
        if(prev != -1)
        {
          m_list[prev].next = m_list[ind].next;
        }
        if(ind == m_head)
        {
          m_head = m_list[ind].next;
        }
        m_list[ind].next = m_free;
        m_free = ind;
        m_size--;
        return true;
      }
      prev = ind;
      ind = m_list[ind].next;
    } while(ind != -1);
    return false;
  }

  /*!
    * \brief Returns pointer to node containing key-value pair matching provided key, if possible.
    *
    * \param [in] key the key of the item to be searched for.
    *
    * \return Returns pointer to the requested node if possible, and to special end node otherwise.
    *
    */
  axom_map::Node<Key, T>& find(const Key& key)
  {
    if(m_head == -1)
    {
      return m_end;
    }
    IndexType ind = m_head;
    do
    {
      if(m_list[ind].key == key)
      {
        return m_list[ind];
      }

      ind = m_list[ind].next;
    } while(ind != -1);
    return m_end;
  }

  /*!
   * \brief Returns maximum number of items that can be stored in this Bucket instance.
   *
   * \return Returns maximum number of items that can be stored in this Bucket instance.
   *
   */
  int get_capacity() const { return m_capacity; }

  /*!
   * \brief Returns current number of items in this Bucket instance.
   *
   * \return Returns current number of items in this Bucket instance.
   */
  int get_size() const { return m_size; }
  ///@}

  axom_map::Node<Key, T>* m_list;
  axom_map::Node<Key, T> m_end;
  IndexType m_head;
  IndexType m_free;
  IndexType m_capacity, m_size;
};

/// @}
}  // namespace axom_map

/*!
 * \class Map
 *
 * \brief Provides a hashmap implementation, relying upon chaining to
 *  resolve hash collisions. Simple hashmap for now, future work will use
 *  RAJA and Umpire to allow use of map within kernels. Resembles
 *  unordered_map, which we can't use due to a lack of host-device decoration.
 *
 * \tparam Key the type of keys. Must allow equality comparison, must be copyable.
 * \tparam T the type of values to hold, must be copyable.
 * \tparam Hash functor that takes an object of Key type and returns
 *  a hashed value of size_t (for now)
 */

template <typename Key, typename T, typename Hash = std::hash<Key>, typename Policy = axom::SEQ_EXEC>
class Map
{
public:
  using key_type = Key;
  using mapped_type = T;

public:
  /// \name Map Constructors
  /// @{

  /*!
   *\brief Constructs a Map instance with the given number of buckets.
   *
   * \param [in] num_buckets the number of buckets used in the Map.
   * \param [in] bucket_len the maximum number of items that can co-exist in a bucket.
   *
   * \pre num_buckets > 0
   * \pre bucket_len > 0
   */
  Map(IndexType num_buckets, IndexType bucket_len = 10)
  {
    init(num_buckets, bucket_len);
    m_end.key = Key {0};
    m_end.value = T {0};
    m_end.next = -2;
  }

  /*!
   *\brief Move constructor for Map
   *
   * \param [in] other The Map to move from
   */
  Map(Map&& other) noexcept
    : m_buckets(nullptr)
    , m_bucket_count(0)
    , m_bucket_len(0)
    , m_size(0)
    , m_load_factor(0)
    , m_bucket_fill(false)
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
    , m_locks(nullptr)
#endif
  {
    *this = std::move(other);
  }

  /// @}

  /*!
   *\brief Move assignment for Map
   *
   * \param [in] other The Map to move from
   */
  Map& operator=(Map&& other) noexcept
  {
    if(this != &other)
    {
      clear();

      m_buckets = other.m_buckets;
      m_bucket_count = other.m_bucket_count;
      m_bucket_len = other.m_bucket_len;
      m_size = other.m_size;
      m_load_factor = other.m_load_factor;
      m_end = other.m_end;
      pol = other.pol;
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
      m_locks = other.m_locks;
#endif

      other.m_buckets = nullptr;
      other.m_bucket_count = 0;
      other.m_bucket_len = 0;
      other.m_size = 0;
      other.m_load_factor = 0;
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
      other.m_locks = nullptr;
#endif
    }
    return *this;
  }

  /// @}

  /*!
   * Destructor. Frees associated buffers.
   */
  ~Map() { clear(); }

  /// \name Map Methods
  ///@{
  /*!
   * \brief Resizes Map by creating new Map instance and inserting all values
   * from original Map.
   *
   * \param [in] buckets user-specified number of buckets in new Map. -1 to fall back to multiplicative factor.
   * \param [in] factor user-specified multiplicative factor to determine size of new Map based upon current size.
   *  -1 to fall back to default -- would be more efficient to specify neither input.
   *
   * \return true if the hash table is now in a safe state, false otherwise
   *
   * \note if neither input is specified, the new Map will be of size 2*old Map size.
   * \note both map instances exist in memory until rehash returns. Risks running out of memory.
   * \note Rehash should generally ensure a previously filled bucket is no longer filled, and the load_factor is correct,
   *  if used properly. However, since this highly manual process has room for error this implementation does
   *  perform a check for if any bucket is now full.
   */
  bool rehash(IndexType buckets = -1, int factor = -1)
  {
    IndexType newlen = 0;
    IndexType ind;
    std::size_t hashed;
    if(buckets != -1 && buckets > m_bucket_count)
    {
      newlen = buckets;
    }
    else if(factor != -1 && factor > 1)
    {
      newlen = factor * m_bucket_count;
    }
    else
    {
      assert(factor == -1);
      assert(buckets == -1);
      newlen = 2 * m_bucket_count;
    }
    m_bucket_fill = false;
    axom_map::Bucket<Key, T>* new_list = alloc_map(newlen, m_bucket_len);

    for(IndexType i = 0; i < m_bucket_count; i++)
    {
      ind = m_buckets[i].m_head;
      while(ind != -1)
      {
        hashed = (get_hash(m_buckets[i].m_list[ind].key) % newlen);
        new_list[hashed].insert_no_update(m_buckets[i].m_list[ind].key,
                                          m_buckets[i].m_list[ind].value);
        ind = m_buckets[i].m_list[ind].next;
      }
      axom::deallocate(m_buckets[i].m_list);
    }
    axom::deallocate(m_buckets);
    destroy_locks(pol);
    m_buckets = new_list;
    m_bucket_count = newlen;
    for(IndexType i = 0; i < m_bucket_count; i++)
    {
      if(m_buckets[i].get_size() == m_buckets[i].get_capacity())
      {
        m_bucket_fill = true;
      }
    }
    init_locks(pol);
    return check_rehash();
  }
  /*!
   * \brief Deallocates memory resources for this Map instance, resets state as if constructed with 0 elements.
   *  Also used in destructor. Once this is used, insert will not function properly until init() is used, because there
   *  is no memory available for elements to be inserted.
   *
   */
  void clear()
  {
    for(IndexType i = 0; i < m_bucket_count; i++)
    {
      axom::deallocate(m_buckets[i].m_list);
    }
    axom::deallocate(m_buckets);
    m_bucket_count = 0;
    m_bucket_len = 0;
    m_size = 0;
    m_load_factor = 0;
    destroy_locks(pol);
  }

  /*!
   * \brief Initializes this Map instance, mostly in case clear() was used and there's a desire to continue use.
   *  Required if use is to continue after a clear().
   *
   */
  void init(IndexType num_buckets, IndexType bucket_len = 10)
  {
    m_bucket_count = num_buckets;
    m_bucket_len = bucket_len;
    m_size = 0;
    m_load_factor = 0;
    m_bucket_fill = false;
    m_buckets = alloc_map(m_bucket_count, m_bucket_len);
    init_locks(pol);
  }
  /*!
   * \brief Inserts given key-value pair into the Map instance.
   *
   * \param [in] key the key used to index into the Map to store this item
   * \param [in] val the data value of the pair to insert into the map
   *
   * \note Deviation from STL unordered_map is that insertion fails when map full, rather than automatically rehashing.
   *  Manual call of rehash() required.
   *
   * \note Can set "full bucket" flag for Map, which is only unset by rehashing.
   *
   * \return A pair whose first element is a pointer to the inserted item and bool set to true
   *  if successful. Otherwise, the second element is set to false, and the first to one of two values: a sentinel node
   *  if the map is overfilled, and the actual node with the given key if an item with given key already exists.
   *
   */
  axom_map::Pair<Key, T> insert(const Key& key, const T& val)
  {
    //If branching becomes an issue, this protective statement will be
    //removed in lieu of expecting the user to not try to insert after
    //clear().
    if(m_bucket_count * m_bucket_len == 0)
    {
      return axom_map::Pair<Key, T>(&m_end, false);
    }
    std::size_t index = bucket(key);
    bucket_lock(index, pol);
    axom_map::Bucket<Key, T>* target = &(m_buckets[index]);
    axom_map::Pair<Key, T> ret = target->insert_no_update(key, val);
    //Candidate to get cut out if branching becomes too much of an issue.
    if(ret.second == true)
    {
      axom::atomicAdd<axom::auto_atomic>(&m_size, IndexType {1});

      if(target->get_size() == target->get_capacity())
      {
        m_bucket_fill = true;
      }
    }
    bucket_unlock(index, pol);
    return ret;
  }

  /*!
   * \brief Inserts given key-value pair into the Map instance, or updates value of existing key-value pair if item
   *  with supplied key already exists.
   *
   * \param [in] key the key used to index into the Map to store this item
   * \param [in] val the data value of the pair to insert into the map
   *
   * \note Deviation from STL unordered_map is that insertion fails when map full, rather than automatically rehashing.
   *  Manual call of rehash() required.
   *
   * \note Can set "full bucket" flag for Map, which is only unset by rehashing.
   *
   * \return A pair whose first element is a pointer to the inserted item and bool set to true
   *  if insertion was performed. Otherwise, the second element is set to false, and the first to one of two values: a sentinel node
   *  if the map is overfilled, and the actual node with the given key if an assignment occured.
   *
   */
  axom_map::Pair<Key, T> insert_or_assign(const Key& key, const T& val)
  {
    //If branching becomes an isssue, this protective statement will be
    //removed in lieu of expecting the user to not try to insert after
    //clear().
    if(m_bucket_count * m_bucket_len == 0)
    {
      return axom_map::Pair<Key, T>(&m_end, false);
    }

    std::size_t index = bucket(key);
    bucket_lock(index, pol);
    axom_map::Bucket<Key, T>* target = &(m_buckets[index]);
    axom_map::Pair<Key, T> ret = target->insert_update(key, val);
    //Candidate to get cut out if branching becomes too much of an issue.
    if(ret.second == true)
    {
      axom::atomicAdd<axom::auto_atomic>(&m_size, IndexType {1});

      if(target->get_size() == target->get_capacity())
      {
        m_bucket_fill = true;
      }
    }
    bucket_unlock(index, pol);
    return ret;
  }

  /*!
   * \brief Finds reference to value of key-value pair mapped to supplied pair, returns value from sentinel node
   *  if fails, meaning value at reference will be data type T initialized with 0.
   *
   * \param [in] key the key used to index into the Map to find item
   *
   * \note Major deviation from STL unordered_map in that this function will not perform insertions, and cannot
   *  stand in for insert(), nor can accesses be guaranteed to point to valid data. Also returns a const reference,
   *  thus meaning it cannot update existing values. Use insert_or_assign() for such functionality.
   *
   * \return A const reference to the value of key-value pair in the Map instance that is associated with supplied Key.
   *
   */
  const T& operator[](const Key& key) const
  {
    //Since we can't throw an exception, the safest solution is to be read-only, and for the user to be careful with their
    //accesses.
    axom_map::Node<Key, T>& ins_result = find(key);

    return ins_result.value;
  }

  /*!
   * \brief Removes key-value pair with given key from the Map instance.
   *
   * \param [in] key the key of the item to be removed from the Map instance.
   *
   * \note Deviation from STL unordered_map in returning a boolean instead of an iterator to next item in map.
   *  Mostly useful on sparing logic when this moves to paralellism, since the bucket boundaries would become inconvenient.
   *
   * \return A bool value of true if the item was successfully removed, false otherwise.
   */
  bool erase(const Key& key)
  {
    std::size_t index = bucket(key);
    bucket_lock(index, pol);
    axom_map::Bucket<Key, T>* target = &(m_buckets[index]);
    //Candidate to get cut out if branching becomes too much of an issue.
    bool ret = target->remove(key);
    if(ret == true)
    {
      axom::atomicSub<axom::auto_atomic>(&m_size, IndexType {1});
    }
    bucket_unlock(index, pol);
    return ret;
  }

  /*!
   * \brief Returns key-value pair with given key from the Map instance.
   *
   * \param [in] key the key of the item to be searched for in the Map instance.
   *
   * \return A reference to the requested item if found, sentinel node end otherwise.
   */
  axom_map::Node<Key, T>& find(const Key& key) const
  {
    std::size_t index = bucket(key);
    bucket_lock(index, pol);
    axom_map::Bucket<Key, T>* target = &(m_buckets[index]);
    axom_map::Node<Key, T>& out = target->find(key);
    bucket_unlock(index, pol);
    return out;
  }

  /*!
   * \brief Returns id of bucket associated with a 64-bit integer, intended to be the hash
   *  of a key.
   *
   * \param [in] hash the id of the bucket to be queried.
   *
   * \note Likely better off inlined, included here in case of later changes to hash collision resolution method.
   *
   * \return A pointer to the queried bucket.
   */
  std::size_t bucket(Key key) const { return get_hash(key) % m_bucket_count; }

  /*!
   * \brief Returns the maximum number of items per bucket.
   * \return bucket_len the maximum number of items per bucket.
   */
  IndexType bucket_size() const { return m_bucket_len; }
  /*!
   * \brief Returns the number of buckets in the Map.
   * \return bucket_count
   */
  IndexType bucket_count() const { return m_bucket_count; }
  /*!
   * \brief Returns the amount of items in the Map instance.
   * \return size the amount of items in the Map instance.
   */
  IndexType size() const { return m_size; }
  /*!
   * \brief Returns the overall capacity of the Map instance.
   * \return capacity the overall capacity of the Map instance.
   */
  IndexType max_size() { return m_bucket_len * m_bucket_count; }
  /*!
   * \brief Checks if the container has no elements.
   * \return true if the container is empty, false otherwise
   */
  bool empty() const { return (m_size == 0); }
  /*!
   * \brief Returns const reference to "end" node, for comparison to check
   *  return results of other Map functions.
   * \return m_end the values of the sentinel node for failure checks
   */
  const axom_map::Node<Key, T>& end() const { return m_end; }
  /*!
   * \brief Returns maximum load factor (ratio between size and bucket count)
   * hash table will reach before check_rehash() returns true.
   *
   * \return max_load_factor maximum load factor desired for this hash table. Default 1.0
   */
  float max_load_factor() const { return m_load_factor; }
  /*!
   * \brief Sets maximum load factor (ratio between size and bucket count)
   * hash table will reach before check_rehash returns true. Default value is
   *  1.0.
   *
   * \param [in] load_factor desired maximum load factor
   */
  void max_load_factor(float load_factor) { m_load_factor = load_factor; }
  /*!
   * \brief Returns current load factor (ratio between size and bucket count).
   *
   * \return load_factor the ratio between the amount of items in the Map and the amount of buckets.
   */
  float load_factor() const { return m_size / m_bucket_count; }

  /*!
   * \brief Returns whether a rehash is necessary for this Map instance. Does so based on two metrics. The
   *  first is whether the load factor of the map (ratio between size and bucket count) has exceeded its
   *  maximum load factor (default 1.0). The second is whether any of the buckets are currently full.
   *
   * \note This function does not exist in the STL unordered_map, because it performs automatic rehashes as
   *  necessary. To support the portability constraints put on this data structure, the user must regularly
   *  call this function to verify whether they need to rehash to maintain performance. Another metric is
   *  when an insertion fails and returns end(), but this will yield an answer before such an event.
   *
   * \return true if Map should be rehashed now, false otherwise
   */
  bool check_rehash() const
  {
    if(m_size / m_bucket_count >= m_load_factor || m_bucket_fill == true)
    {
      return true;
    }
    return false;
  }
  ///@}
private:
  /// \name Private Map Methods
  ///@{
  /*!
   * \brief Memory allocation handler, for Map construction, both at initialization and rehash. Allocates and
   *  initializes linked lists used for data storage, and returns pointer to said array, to be used elsewhere.
   *  Kept separate from constructor and rehasher for sake of modularity in case this needs to be swapped
   *  with a different method in the future.
   *
   * \param [in] bucount the number of buckets to be in side the allocated Map instance.
   * \param [in] bucklen the amount of items that can fit in a bucket in the allocated Map instance.
   *
   * \return A pointer to the now-allocated array of linked lists.
   */
  axom_map::Bucket<Key, T>* alloc_map(IndexType bucount, IndexType bucklen)
  {
    axom_map::Bucket<Key, T>* tmp = axom::allocate<axom_map::Bucket<Key, T>>(bucount);
    for(IndexType i = 0; i < bucount; i++)
    {
      tmp[i].init(bucklen);
    }
    //update when we're testing our returns
    return tmp;
  }
  /*!
   * \brief Returns hash value for a given input.
   *
   * \note Currently uses std::hash, will switch to a custom hasher in future.
   *
   * \return A 64-bit integer containing the hashed value of the key.
   */
  std::size_t get_hash(const Key& input) const
  {
    std::size_t hashed = Hash {}(input);

    return hashed;
  }

  //NOTE: Every single locking function will use the OpenMP locks if Axom was compiled with OpenMP support,
  //regardless of whether the user code is actually using it.

  /*!
   * \brief Empty function for sequential execution environment.
   *
   * \param [in] index the index needed for locking, if this weren't sequential.
   * \param [in] overload execution space object for the sake of function overloading.
   */
  void bucket_lock(std::size_t index, axom::SEQ_EXEC overload) const
  {
    AXOM_UNUSED_VAR(index);
    AXOM_UNUSED_VAR(overload);
  }

  /*!
   * \brief Empty function for sequential execution environment.
   *
   * \param [in] index the index needed for locking, if this weren't sequential.
   * \param [in] overload execution space object for the sake of function overloading.
   */
  void bucket_unlock(std::size_t index, axom::SEQ_EXEC overload) const
  {
    AXOM_UNUSED_VAR(index);
    AXOM_UNUSED_VAR(overload);
  }

  /*!
   * \brief Empty function for sequential execution environment.
   *
   * \param [in] overload execution space object for the sake of function overloading.
   */
  void init_locks(axom::SEQ_EXEC overload) const { AXOM_UNUSED_VAR(overload); }

  /*!
   * \brief Empty function for sequential execution environment.
   *
   * \param [in] overload execution space object for the sake of function overloading.
   */
  void destroy_locks(axom::SEQ_EXEC overload) const { AXOM_UNUSED_VAR(overload); }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)

  /*!
   * \brief Acquires lock for a bucket.
   *
   * \param [in] index the index of the bucket which the thread wants the lock for.
   */
  void bucket_lock(std::size_t index, axom::OMP_EXEC overload) const
  {
    AXOM_UNUSED_VAR(overload);
    omp_set_lock(m_locks + index);
  }

  /*!
   * \brief Releases lock for a bucket.
   *
   * \param [in] index the index of the bucket which the thread wants to release the lock for.
   */
  void bucket_unlock(std::size_t index, axom::OMP_EXEC overload) const
  {
    AXOM_UNUSED_VAR(overload);
    omp_unset_lock(m_locks + index);
  }

  /*!
   * \brief Allocates and initializes a lock for every bucket making up the Map instance.
   */
  void init_locks(axom::OMP_EXEC overload) const
  {
    AXOM_UNUSED_VAR(overload);
    m_locks = axom::allocate<omp_lock_t>(m_bucket_count);
    for(IndexType i = 0; i < m_bucket_count; i++)
    {
      omp_init_lock(m_locks + i);
    }
  }

  /*!
   * \brief Deallocates and destroys every lock used for this Map instance.
   */
  void destroy_locks(axom::OMP_EXEC overload) const
  {
    AXOM_UNUSED_VAR(overload);
    for(IndexType i = 0; i < m_bucket_count; i++)
    {
      omp_destroy_lock(m_locks + i);
    }
    axom::deallocate<omp_lock_t>(m_locks);
  }
#endif
  /// @}

  /// \name Private Data Members
  /// @{

  axom_map::Bucket<Key, T>* m_buckets; /*!< array of pointers to linked lists containing data */
  IndexType m_bucket_count;            /*!< the number of buckets in the Map instance */
  IndexType m_bucket_len; /*!< the number of items that can be contained in a bucket in this Map instance */
  IndexType m_size;    /*!< the number of items currently stored in this Map instance */
  float m_load_factor; /*!< currently unused value, used in STL unordered_map to determine when to resize, which we don't do internally at the moment */
  axom_map::Node<Key, T> m_end; /*!< the sentinel node enabling verification of operation success or failure */
  bool m_bucket_fill; /*!<  status of buckets in general -- if at least one is full, this is set to true, false otherwise*/
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  mutable omp_lock_t* m_locks;
#endif
  Policy pol;
  /// @}
};
} /* namespace experimental */
} /* namespace axom */

#endif /* AXOM_MAP_HPP_ */
