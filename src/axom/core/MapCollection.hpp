// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file MapCollection.hpp
 *
 * \brief   Header file for MapCollection.
 *
 *          MapCollection is an implemenation of ItemCollection to
 *          hold a collection of items of a fixed type that can be accessed
 *          accessed by string name or axom::IndexType.
 *
 *          The primary intent is to decouple the implementation of the
 *          collection of times from the Group class which owns collections of
 *          View and child Group objects. This may have other uses,
 *          so it is not dependent on the Group class. This class is
 *          templated on the item type so that it can be used
 *          to hold either View or Group object pointers without
 *          having to code a separate class for each.
 *
 *          \attention This class should be robust against any potential
 *                     user interaction. It doesn't report errors and leaves
 *                     checking of return values to calling code.
 *
 *          \attention Template parameter type must provide a method
 *                     "getName()" that returns a reference to a string object.
 *
 *          \attention The interface is as follows:
 *
 *          \verbatim
 *
 *          - // Return number of items in collection.
 *
 *               size_t getNumItems() const;
 *
 *          - // Return first valid item index (i.e., smallest index over
 *            // all items).
 *            // axom::InvalidIndex returned if no items in collection
 *
 *               IndexType getFirstValidIndex() const;
 *
 *          - // Return next valid item index after given index (i.e., smallest
 *            // index over all indices larger than given one).
 *            // axom::InvalidIndex returned
 *
 *               IndexType getNextValidIndex(IndexType idx) const;
 *
 *          - // Return true if item with given name in collection; else false.
 *
 *               bool hasItem(const std::string& name) const;
 *
 *          - // Return true if item with given index in collection; else false.
 *
 *               bool hasItem(IndexType idx) const;
 *
 *          - // Return pointer to item with given name (nullptr if none).
 *
 *               T* getItem(const std::string& name);
 *               T const* getItem(const std::string& name) const ;
 *
 *          - // Return pointer to item with given index (nullptr if none).
 *
 *               T* getItem(IndexType idx);
 *               T const* getItem(IndexType idx) const;
 *
 *          - // Return name of object with given index
 *            // (axom::utilities::string::InvalidName if none).
 *
 *               const std::string& getItemName(IndexType idx) const;
 *
 *          - // Return index of object with given name
 *            // (axom::InvalidIndex if none).
 *
 *               IndexType getItemIndex(const std::string& name) const;
 *
 *          - // Insert item with given name; return index if insertion
 *            // succeeded, and InvalidIndex otherwise.
 *
 *               IndexType insertItem(T* item, const std::string& name);
 *
 *          - // Remove item with given name if it exists and return a
 *            // pointer to it. If it doesn't exist, return nullptr.
 *
 *               T* removeItem(const std::string& name);
 *
 *          - // Remove item with given index if it exists and return a
 *            // pointer to it. If it doesn't exist, return nullptr.
 *
 *               T* removeItem(IndexType idx);
 *
 *          - // Remove all items (items not destroyed).
 *
 *               void removeAllItems();
 *
 *          - // Clear all items and destroy them.
 *
 *               void deleteAllItems();
 *
 *          \endverbatim
 *
 ******************************************************************************
 */

#ifndef AXOM_MAP_COLLECTIONS_HPP_
#define AXOM_MAP_COLLECTIONS_HPP_

// Standard C++ headers
#include <map>
#include <stack>
#include <string>
#include <vector>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

// Sidre project headers
#include "ItemCollection.hpp"

#if defined(AXOM_USE_SPARSEHASH)
  #include "axom/sparsehash/dense_hash_map"
#else
  #include <unordered_map>
#endif

namespace axom
{
////////////////////////////////////////////////////////////////////////
//
// MapCollection keeps an index constant for each item
// as long as it remains in the collection; i.e., don't shift indices
// around.  It has the additional benefit that users can hold on to
// item indices without them being changed without notice.
//
////////////////////////////////////////////////////////////////////////

/*!
 *************************************************************************
 *
 * \class MapCollection
 *
 * \brief MapCollection is a container class template for holding
 *        a collection of items of template parameter type T
 *
 * \warning Only std::map and std::unordered_map have been tried so far.
 *          These classes have identical APIs for the functionality we
 *          are using.
 *
 *************************************************************************
 */
template <typename T>
class MapCollection : public ItemCollection<T>
{
public:
  using value_type = T;
  using iterator = typename ItemCollection<T>::iterator;
  using const_iterator = typename ItemCollection<T>::const_iterator;

public:
  //
  // Default compiler-generated ctor, dtor, copy ctor, and copy assignment
  // operator suffice for this class.
  //

  ///
  size_t getNumItems() const { return m_items.size() - m_free_ids.size(); }

  ///
  IndexType getFirstValidIndex() const;

  ///
  IndexType getNextValidIndex(IndexType idx) const;

  ///
  bool hasItem(const std::string& name) const
  {
    // Implementation note: Sparsehash requires sentinel values for 'empty' and 'deleted'
    // entries. Our usage in this class only enforces this when we insert or delete
    // items from the map, so we add explicit checks for empty maps in functions
    // that do not mutate the map to avoid undefined behavior.
    if(m_name2idx_map.empty())
    {
      return false;
    }
    else
    {
      const auto mit = m_name2idx_map.find(name);
      return (mit != m_name2idx_map.end() ? true : false);
    }
  }

  ///
  bool hasItem(IndexType idx) const
  {
    return (idx >= 0 && static_cast<std::size_t>(idx) < m_items.size() &&
            m_items[static_cast<std::size_t>(idx)]);
  }

  ///
  T* getItem(const std::string& name)
  {
    // Implementation note: This explicit check avoids undefined behavior. See note above.
    if(m_name2idx_map.empty())
    {
      return nullptr;
    }
    else
    {
      const auto mit = m_name2idx_map.find(name);
      return (mit != m_name2idx_map.end() ? m_items[mit->second] : nullptr);
    }
  }

  ///
  T const* getItem(const std::string& name) const
  {
    // Implementation note: This explicit check avoids undefined behavior. See note above.
    if(m_name2idx_map.empty())
    {
      return nullptr;
    }
    else
    {
      const auto mit = m_name2idx_map.find(name);
      return (mit != m_name2idx_map.end() ? m_items[mit->second] : nullptr);
    }
  }

  ///
  T* getItem(IndexType idx)
  {
    return (hasItem(idx) ? m_items[static_cast<std::size_t>(idx)] : nullptr);
  }

  ///
  T const* getItem(IndexType idx) const
  {
    return (hasItem(idx) ? m_items[static_cast<std::size_t>(idx)] : nullptr);
  }

  ///
  const std::string& getItemName(IndexType idx) const
  {
    return (hasItem(idx) ? m_items[static_cast<std::size_t>(idx)]->getName()
                         : axom::utilities::string::InvalidName);
  }

  ///
  IndexType getItemIndex(const std::string& name) const
  {
    // Implementation note: This explicit check avoids undefined behavior. See note above.
    if(m_name2idx_map.empty())
    {
      return InvalidIndex;
    }
    else
    {
      const auto mit = m_name2idx_map.find(name);
      return (mit != m_name2idx_map.end() ? mit->second : InvalidIndex);
    }
  }

  ///
  IndexType insertItem(T* item, const std::string& name);

  ///
  T* removeItem(const std::string& name);

  ///
  T* removeItem(IndexType idx);

  ///
  void removeAllItems()
  {
    m_items.clear();
    while(!m_free_ids.empty())
    {
      m_free_ids.pop();
    }
#if defined(AXOM_USE_SPARSEHASH)
    if(m_name2idx_map.empty() && m_empty_key != "DENSE_MAP_EMPTY_KEY")
    {
      m_empty_key = "DENSE_MAP_EMPTY_KEY";
      m_name2idx_map.set_empty_key(m_empty_key);
      m_name2idx_map.set_deleted_key("DENSE_MAP_DELETED_KEY");
    }
#endif

    m_name2idx_map.clear();
  }

  iterator begin() { return iterator(this, true); }
  iterator end() { return iterator(this, false); }

  const_iterator cbegin() const { return const_iterator(this, true); }
  const_iterator cend() const { return const_iterator(this, false); }

  const_iterator begin() const { return const_iterator(this, true); }
  const_iterator end() const { return const_iterator(this, false); }

private:
  std::vector<T*> m_items;
  std::stack<IndexType> m_free_ids;

#if defined(AXOM_USE_SPARSEHASH)
  using MapType = axom::google::dense_hash_map<std::string, IndexType>;
#else
  using MapType = std::unordered_map<std::string, IndexType>;
#endif

  MapType m_name2idx_map;
#if defined(AXOM_USE_SPARSEHASH)
  std::string m_empty_key;
#endif
};

template <typename T>
IndexType MapCollection<T>::getFirstValidIndex() const
{
  IndexType idx = 0;
  while(static_cast<std::size_t>(idx) < m_items.size() &&
        m_items[static_cast<std::size_t>(idx)] == nullptr)
  {
    ++idx;
  }
  return ((static_cast<std::size_t>(idx) < m_items.size()) ? idx : InvalidIndex);
}

template <typename T>
IndexType MapCollection<T>::getNextValidIndex(IndexType idx) const
{
  if(idx == InvalidIndex)
  {
    return InvalidIndex;
  }

  idx++;
  while(static_cast<std::size_t>(idx) < m_items.size() &&
        m_items[static_cast<std::size_t>(idx)] == nullptr)
  {
    idx++;
  }
  return ((static_cast<std::size_t>(idx) < m_items.size()) ? idx : InvalidIndex);
}

template <typename T>
IndexType MapCollection<T>::insertItem(T* item, const std::string& name)
{
  bool use_recycled_index = false;
  IndexType idx = m_items.size();
  if(!m_free_ids.empty())
  {
    idx = m_free_ids.top();
    m_free_ids.pop();
    use_recycled_index = true;
  }

#if defined(AXOM_USE_SPARSEHASH)
  if(m_name2idx_map.empty() && m_empty_key != "DENSE_MAP_EMPTY_KEY")
  {
    m_empty_key = "DENSE_MAP_EMPTY_KEY";
    m_name2idx_map.set_empty_key(m_empty_key);
    m_name2idx_map.set_deleted_key("DENSE_MAP_DELETED_KEY");
  }
#endif

  if(m_name2idx_map.insert(std::make_pair(name, idx)).second)
  {
    // name was inserted into map
    if(use_recycled_index)
    {
      m_items[idx] = item;
    }
    else
    {
      m_items.push_back(item);
    }
    return idx;
  }
  else
  {
    // name was NOT inserted into map, return free index if necessary
    if(use_recycled_index)
    {
      m_free_ids.push(idx);
    }
    return InvalidIndex;
  }
}

template <typename T>
T* MapCollection<T>::removeItem(const std::string& name)
{
  T* ret_val = nullptr;

  auto mit = m_name2idx_map.find(name);
  if(mit != m_name2idx_map.end())
  {
    IndexType idx = mit->second;

    ret_val = m_items[idx];

    m_name2idx_map.erase(mit);
    m_items[idx] = nullptr;
    m_free_ids.push(idx);
  }

  return ret_val;
}

template <typename T>
T* MapCollection<T>::removeItem(IndexType idx)
{
  if(hasItem(idx))
  {
    T* item = removeItem(m_items[idx]->getName());
    return item;
  }
  else
  {
    return nullptr;
  }
}

}  // namespace axom

#endif  // AXOM_MAP_COLLECTIONS_HPP_
