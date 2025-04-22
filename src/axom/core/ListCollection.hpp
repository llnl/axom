// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file ListCollection.hpp
 *
 * \brief   Header file for ListCollection.
 *
 *          This is an implementation of ItemCollection to hold a
 *          collection of items of a fixed type. This implementation
 *          is intended to hold items that may have no name. If they do
 *          have names, those names are ignored.
 *
 *          This class is templated on the item type so that the same
 *          class can be used to hold either View or Group object pointers
 *          without having to code a separate class for each.
 *
 *          \attention This class should be robust against any potential
 *                     user interaction. It doesn't report errors and leaves
 *                     checking of return values to calling code.
 *
 *          \attention The interface is as follows:
 *
 *          \verbatim
 *
 *          - // Return number of items in collection.
 *
 *               size_t getNumItems() const;
 *
 *          - // Return first item index for iteration.
 *            // axom::InvalidIndex returned if no items in collection
 *
 *               IndexType getFirstValidIndex() const;
 *
 *          - // Return next valid item index for iteration.
 *            // axom::InvalidIndex returned if there are no further items
 *
 *               IndexType getNextValidIndex(IndexType idx) const;
 * *
 *          - // Return true if item with given index in collection; else false.
 *
 *               bool hasItem(IndexType idx) const;
 *
 *          - // Return pointer to item with given index (nullptr if none).
 *
 *               T* getItem(IndexType idx);
 *               T const* getItem(IndexType idx) const;
 *
 *
 *          - // Insert item; the name argument will be ignored.
 *            // Return index if insertion succeeded, and InvalidIndex
 *            // otherwise.
 *
 *               IndexType insertItem(T* item, const std::string& name);
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

#ifndef AXOM_LISTCOLLECTIONS_HPP_
#define AXOM_LISTCOLLECTIONS_HPP_

// Standard C++ headers
#include <iostream>
#include <map>
#include <list>
#include <stack>
#include <string>
#include <vector>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Types.hpp"

// Sidre project headers
#include "axom/core/ItemCollection.hpp"

namespace axom
{
////////////////////////////////////////////////////////////////////////
//
// ListCollection keeps an index constant for each item
// as long as it remains in the collection; i.e., don't shift indices
// around.  It has the additional benefit that users can hold on to
// item indices without them being changed without notice.
//
////////////////////////////////////////////////////////////////////////

/*!
 *************************************************************************
 *
 * \class ListCollection
 *
 * \brief ListCollection is a container class template for holding
 *        a collection of items of template parameter type T, using
 *        a list container.
 *
 *************************************************************************
 */
template <typename T>
class ListCollection : public ItemCollection<T>
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
  bool hasItem(IndexType idx) const
  {
    return (idx >= 0 && static_cast<unsigned>(idx) < m_items.size() &&
            m_items[static_cast<unsigned>(idx)]);
  }

  ///
  T* getItem(IndexType idx)
  {
    return (hasItem(idx) ? m_items[static_cast<unsigned>(idx)] : nullptr);
  }

  ///
  T const* getItem(IndexType idx) const
  {
    return (hasItem(idx) ? m_items[static_cast<unsigned>(idx)] : nullptr);
  }

  ///
  IndexType insertItem(T* item, const std::string& name = "");

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
    m_index_list.clear();
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

  std::list<IndexType> m_index_list;
};

template <typename T>
IndexType ListCollection<T>::getFirstValidIndex() const
{
  IndexType idx = 0;
  while(static_cast<unsigned>(idx) < m_items.size() && m_items[static_cast<unsigned>(idx)] == nullptr)
  {
    idx++;
  }
  return ((static_cast<unsigned>(idx) < m_items.size()) ? idx : InvalidIndex);
}

template <typename T>
IndexType ListCollection<T>::getNextValidIndex(IndexType idx) const
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
IndexType ListCollection<T>::insertItem(T* item, const std::string& name)
{
#ifdef AXOM_DEBUG
  if(!name.empty())
  {
    std::cerr << "Item " << name << " added to Group "
              << "which holds items in list format. "
              << "The name of this item will be ignored." << std::endl;
  }
#else
  AXOM_UNUSED_VAR(name);
#endif

  bool use_recycled_index = false;
  IndexType idx = m_items.size();
  if(!m_free_ids.empty())
  {
    idx = m_free_ids.top();
    m_free_ids.pop();
    use_recycled_index = true;
  }

  m_index_list.push_back(idx);

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

template <typename T>
T* ListCollection<T>::removeItem(IndexType idx)
{
  T* ret_val = nullptr;
  if(hasItem(idx))
  {
    for(auto itr = m_index_list.begin(); itr != m_index_list.end(); ++itr)
    {
      if(*itr == idx)
      {
        ret_val = m_items[idx];
        m_index_list.erase(itr);
        m_items[idx] = nullptr;
        m_free_ids.push(idx);
        break;
      }
    }
  }

  return ret_val;
}

} /* end namespace axom */

#endif /* AXOM_LIST_COLLECTIONS_HPP_ */
