// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SIDRE_INDEXED_COLLECTION_HPP_
#define SIDRE_INDEXED_COLLECTION_HPP_

// Standard C++ headers
#include <iostream>
#include <map>
#include <stack>
#include <string>
#include <vector>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/ItemCollection.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

namespace axom
{
/*!
 *************************************************************************
 *
 * \class IndexedCollection
 *
 * \brief IndexedCollection is a container for a collection of pointers
 * to items of template parameter type T, each with a corresponding index
 *
 * Each item has an associated index which will always be in the range
 * between 0 and \a getLastAvailableEmptyIndex()
 *************************************************************************
 */
template <typename T>
class IndexedCollection : public ItemCollection<T>
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

  /// Gets the number of items stored in the collection
  size_t getNumItems() const { return m_num_items; }

  /// Returns the index of the first valid item or InvalidIndex if there are none
  IndexType getFirstValidIndex() const;

  /// Returns the next valid index after \idx or InvalidIndex if there are none
  IndexType getNextValidIndex(IndexType idx) const;

  /// Return true if there is an item at index \a idx
  bool hasItem(IndexType idx) const
  {
    return isInHalfOpenRange(idx) && m_items[static_cast<unsigned>(idx)] != nullptr;
  }

  /// Return the \a item at index \idx or nullptr if that index is empty
  T* getItem(IndexType idx)
  {
    return (hasItem(idx) ? m_items[static_cast<unsigned>(idx)] : nullptr);
  }

  /// Return the \a item at index \idx or nullptr if that index is empty
  T const* getItem(IndexType idx) const
  {
    return (hasItem(idx) ? m_items[static_cast<unsigned>(idx)] : nullptr);
  }

  /// Insert \a item into the next available free index
  IndexType insertItem(T* item) { return insertItem(item, getValidEmptyIndex()); }

  /*!
   * \brief Insert \a item into the next available free index
   *
   * \note The second parameter is unused and only present to conform to
   * the ItemCollection interface
   */
  IndexType insertItem(T* item, const std::string& AXOM_UNUSED_PARAM(name))
  {
    return insertItem(item);
  }

  /*!
   * \brief  Insert \a item at index \a idx if that index is not already occupied
   *
   * \return Index at which \a item was inserted, if successful; axom::InvalidIndex otherwise
   */
  IndexType insertItem(T* item, IndexType idx)
  {
    if(hasItem(idx))
    {
      return axom::InvalidIndex;
    }

    if(idx < 0)
    {
      return axom::InvalidIndex;
    }

    // grow capacity to support insertion at index
    if(!isInHalfOpenRange(idx))
    {
      m_items.reserve(idx);
      for(auto i = getLastAvailableEmptyIndex(); i < idx; ++i)
      {
        m_free_ids.push(i);
        m_items.push_back(nullptr);
      }
    }

    if(idx == getLastAvailableEmptyIndex())
    {
      m_items.push_back(item);
      ++m_num_items;
    }
    else
    {
      // Remove this index from the free_ids stack if it's at the top
      if(!m_free_ids.empty() && m_free_ids.top() == idx)
      {
        m_free_ids.pop();
      }
      if(m_items[idx] == nullptr)
      {
        ++m_num_items;
      }
      m_items[idx] = item;
    }

    return idx;
  }

  /// Removes the item from index \a idx but does not destroy it
  T* removeItem(IndexType idx);

  /*!
   * \brief Removes all items from the collection, but does not destroy them
   *
   * \warning This function can leak memory if the collection stores
   * a pointer to the only copy of the items
   */
  void removeAllItems()
  {
    m_items.clear();
    while(!m_free_ids.empty())
    {
      m_free_ids.pop();
    }
    m_num_items = 0;
  }

  /*!
   * \brief Return the index of a valid empty slot in the collection
   *
   * Finds an empty (unused) index at which an item can be inserted
   */
  IndexType getValidEmptyIndex()
  {
    IndexType newIndex = axom::InvalidIndex;
    bool found_empty_index = false;

    // try to find an empty index from the stack
    while(!m_free_ids.empty() && !found_empty_index)
    {
      newIndex = m_free_ids.top();
      if(hasItem(newIndex))
      {
        m_free_ids.pop();
      }
      else
      {
        found_empty_index = true;
      }
    }

    // if empty index not found, extend the array
    if(!found_empty_index)
    {
      newIndex = m_items.size();
    }

#ifdef AXOM_DEBUG
    if(!isInClosedRange(newIndex) || hasItem(newIndex))
    {
      std::cerr << "Index " << newIndex << " in IndexedCollection is not a valid empty index"
                << std::endl;
    }
    assert(isInClosedRange(newIndex) && !hasItem(newIndex));
#endif

    return newIndex;
  }

  /*!
   * \brief Gets the empty index at the end of the range of available indices
   *
   * \note This index will always be empty
   */
  IndexType getLastAvailableEmptyIndex() const { return m_items.size(); }

  iterator begin() { return iterator(this, true); }
  iterator end() { return iterator(this, false); }

  const_iterator cbegin() const { return const_iterator(this, true); }
  const_iterator cend() const { return const_iterator(this, false); }

  const_iterator begin() const { return const_iterator(this, true); }
  const_iterator end() const { return const_iterator(this, false); }

private:
  /// Predicate to check if index \a idx is in the half-open range, 0 <= idx < getLastAvailableEmptyIndex()
  bool isInHalfOpenRange(IndexType idx) const
  {
    return idx >= 0 && idx < getLastAvailableEmptyIndex();
  }

  /// Predicate to check if index \a idx is in the closed range, 0 <= idx <= getLastAvailableEmptyIndex()
  bool isInClosedRange(IndexType idx) const
  {
    return idx >= 0 && idx <= getLastAvailableEmptyIndex();
  }

private:
  std::vector<T*> m_items;
  std::stack<IndexType> m_free_ids;
  int m_num_items {0};
};

// -----------------------------------------------------------------------------

template <typename T>
IndexType IndexedCollection<T>::getFirstValidIndex() const
{
  IndexType idx = 0;
  while(static_cast<unsigned>(idx) < m_items.size() && m_items[static_cast<unsigned>(idx)] == nullptr)
  {
    ++idx;
  }
  return ((static_cast<unsigned>(idx) < m_items.size()) ? idx : InvalidIndex);
}

template <typename T>
IndexType IndexedCollection<T>::getNextValidIndex(IndexType idx) const
{
  if(idx == InvalidIndex)
  {
    return InvalidIndex;
  }

  idx++;
  while(static_cast<unsigned>(idx) < m_items.size() && m_items[static_cast<unsigned>(idx)] == nullptr)
  {
    idx++;
  }
  return ((static_cast<unsigned>(idx) < m_items.size()) ? idx : InvalidIndex);
}

template <typename T>
T* IndexedCollection<T>::removeItem(IndexType idx)
{
  if(hasItem(idx))
  {
    T* item = m_items[idx];
    m_items[idx] = nullptr;
    m_free_ids.push(idx);
    --m_num_items;

    return item;
  }
  return nullptr;
}

}  // end namespace axom

#endif  // SIDRE_INDEXED_COLLECTION_HPP_
