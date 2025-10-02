// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file ItemCollection.hpp
 *
 * \brief   Header file for ItemCollection.
 *
 *          This is a templated abstract base class defining an interface for
 *          classes holding a collection of items of a fixed
 *          type that can be accessed by string name or axom::IndexType.
 *
 *          The primary intent is to decouple the implementation of the
 *          collections from the Group class which owns collections of
 *          View and child Group objects. They may have other uses,
 *          so they are not dependent on the Group class. This class is
 *          templated on the item type so that derived classes can be used
 *          to hold either View or Group object pointers without
 *          having to code a separate class for each.
 *
 *          Derived implemenations of this class can be used to explore
 *          alternative collection implementations for performance
 *          (insertion, lookup, etc.) and memory overhead.
 *
 *          \attention These classes should be robust against any potential
 *                     user interaction. They don't report errors and leave
 *                     checking of return values to calling code.
 *
 *          \attention The interface defined by this class is as follows:
 *
 *          \verbatim
 *
 *          - // Return number of items in collection.
 *
 *               size_t getNumItems() const;
 *
 *          - // Return first valid item index for iteration.
 *            // axom::InvalidIndex returned if no items in collection
 *
 *               IndexType getFirstValidIndex() const;
 *
 *          - // Return next valid item index for iteration.
 *            // axom::InvalidIndex returned if there are no more items
 *            // to be iterated over.
 *
 *               IndexType getNextValidIndex(IndexType idx) const;
 *
 *          - // Return true if item with given index in collection; else false.
 *
 *               bool hasItem(IndexType idx) const;
 *
 *          - // Return pointer to item with given index (nullptr if none).
 *
 *               T* getItem(IndexType idx);
 *               T const* getItem(IndexType idx) const;
 *
 *          - // Insert item with given name; return index if insertion
 *            // succeeded, and InvalidIndex otherwise.
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
 *          \endverbatim
 *
 ******************************************************************************
 */

#ifndef AXOM_ITEMCOLLECTIONS_HPP_
#define AXOM_ITEMCOLLECTIONS_HPP_

#include <string>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/IteratorBase.hpp"

namespace axom
{
/*!
 *************************************************************************
 *
 * \class ItemCollection
 *
 * \brief ItemCollection is an abstract base class template for holding
 *        a collection of items of template parameter type T.  Derived
 *        child classes can determine how to specifically store the items.
 *
 *************************************************************************
 */
template <typename T>
class ItemCollection
{
public:
  using value_type = T;

  // Forward declare iterator classes and helpers
  class iterator;
  class const_iterator;
  class iterator_adaptor;
  class const_iterator_adaptor;

public:
  virtual ~ItemCollection() { }

  //
  // Default compiler-generated ctor, dtor, copy ctor, and copy assignment
  // operator suffice for this class.
  //

  ///
  virtual size_t getNumItems() const = 0;

  ///
  virtual IndexType getFirstValidIndex() const = 0;

  ///
  virtual IndexType getNextValidIndex(IndexType idx) const = 0;

  ///
  virtual bool hasItem(IndexType idx) const = 0;

  ///
  virtual T* getItem(IndexType idx) = 0;

  ///
  virtual T const* getItem(IndexType idx) const = 0;

  ///
  virtual IndexType insertItem(T* item, const std::string& name = "") = 0;

  ///
  virtual T* removeItem(IndexType idx) = 0;

  ///
  virtual void removeAllItems() = 0;

public:
  virtual iterator begin() = 0;
  virtual iterator end() = 0;

  virtual const_iterator cbegin() const = 0;
  virtual const_iterator cend() const = 0;

  virtual const_iterator begin() const = 0;
  virtual const_iterator end() const = 0;

  /// Returns an adaptor wrapping this collection in support of iteration
  iterator_adaptor getIteratorAdaptor() { return iterator_adaptor(this); }

  /// Returns a const adaptor wrapping this collection in support of iteration
  const_iterator_adaptor getIteratorAdaptor() const { return const_iterator_adaptor(this); }
};

/*!
 * \brief An std-compliant forward iterator for an ItemCollection
 */
template <typename T>
class ItemCollection<T>::iterator : public IteratorBase<iterator, IndexType>
{
private:
  using BaseType = IteratorBase<iterator, IndexType>;
  using CollectionType = ItemCollection<T>;

public:
  // Iterator traits required to satisfy LegacyRandomAccessIterator concept
  // before C++20
  // See: https://en.cppreference.com/w/cpp/iterator/iterator_traits
  using difference_type = IndexType;
  using value_type = typename std::remove_cv<T>::type;
  using reference = T&;
  using pointer = T*;
  using iterator_category = std::forward_iterator_tag;

public:
  iterator(CollectionType* coll, bool is_first) : m_collection(coll)
  {
    assert(coll != nullptr);

    BaseType::m_pos = is_first ? coll->getFirstValidIndex() : axom::InvalidIndex;
  }

  IndexType index() const { return BaseType::m_pos; }

  pointer operator->() { return m_collection->getItem(BaseType::m_pos); }

  reference operator*() { return *m_collection->getItem(BaseType::m_pos); }

private:
  // Remove backwards iteration functions
  using BaseType::operator--;
  using BaseType::operator-=;

protected:
  /// Implementation of advance() as required by IteratorBase
  void advance(IndexType n)
  {
    for(int i = 0; i < n; ++i)
    {
      BaseType::m_pos = m_collection->getNextValidIndex(BaseType::m_pos);
    }
  }

private:
  CollectionType* m_collection;
};

/*!
 * \brief An std-compliant forward iterator for a const ItemCollection
 */
template <typename T>
class ItemCollection<T>::const_iterator : public IteratorBase<const_iterator, IndexType>
{
private:
  using BaseType = IteratorBase<const_iterator, IndexType>;
  using CollectionType = ItemCollection<T>;

public:
  // Iterator traits required to satisfy LegacyRandomAccessIterator concept
  // before C++20
  // See: https://en.cppreference.com/w/cpp/iterator/iterator_traits
  using difference_type = IndexType;
  using value_type = typename std::remove_cv<T>::type;
  using reference = const T&;
  using pointer = const T*;
  using iterator_category = std::forward_iterator_tag;

public:
  const_iterator(const CollectionType* coll, bool is_first) : m_collection(coll)
  {
    assert(coll != nullptr);

    BaseType::m_pos = is_first ? coll->getFirstValidIndex() : axom::InvalidIndex;
  }

  IndexType index() const { return BaseType::m_pos; }

  pointer operator->() { return m_collection->getItem(BaseType::m_pos); }

  reference operator*() { return *m_collection->getItem(BaseType::m_pos); }

private:
  // Remove backwards iteration functions
  using BaseType::operator--;
  using BaseType::operator-=;

protected:
  /// Implementation of advance() as required by IteratorBase
  void advance(IndexType n)
  {
    for(int i = 0; i < n; ++i)
    {
      BaseType::m_pos = m_collection->getNextValidIndex(BaseType::m_pos);
    }
  }

private:
  const CollectionType* m_collection;
};

/*!
 * \brief Utility class to wrap an ItemCollection in support of iteration
 */
template <typename T>
class ItemCollection<T>::iterator_adaptor
{
public:
  using CollectionType = ItemCollection<T>;

public:
  iterator_adaptor(CollectionType* coll) : m_collection(coll) { }

  std::size_t size() const { return m_collection ? m_collection->getNumItems() : 0; }

  iterator begin() { return iterator(m_collection, true); }
  iterator end() { return iterator(m_collection, false); }

  const_iterator cbegin() const { return const_iterator(m_collection, true); }
  const_iterator cend() const { return const_iterator(m_collection, false); }

  operator const_iterator_adaptor() const { return const_iterator_adaptor(m_collection); }

private:
  CollectionType* m_collection {nullptr};
};

/*!
 * \brief Utility class to wrap a const ItemCollection in support of iteration
 */
template <typename T>
class ItemCollection<T>::const_iterator_adaptor
{
public:
  using CollectionType = ItemCollection<T>;

public:
  const_iterator_adaptor(const CollectionType* coll) : m_collection(coll) { }

  std::size_t size() const { return m_collection ? m_collection->getNumItems() : 0; }

  const_iterator begin() { return const_iterator(m_collection, true); }
  const_iterator end() { return const_iterator(m_collection, false); }

  const_iterator cbegin() const { return const_iterator(m_collection, true); }
  const_iterator cend() const { return const_iterator(m_collection, true); }

private:
  const CollectionType* m_collection {nullptr};
};

} /* end namespace axom */

#endif /* AXOM_ITEMCOLLECTIONS_HPP_ */
