// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file DynamicSet.hpp
 *
 * \brief Contains a DynamicSet class, whose size can change dynamically
 * at runtime
 */

#ifndef SLAM_DYNAMIC_SET_H_
#define SLAM_DYNAMIC_SET_H_

#include "axom/config.hpp"
#include "axom/core/IteratorBase.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/RangeSet.hpp"

namespace axom
{
namespace slam
{
/**
 * \class DynamicSet
 * \brief A Set class that supports dynamically adding and removing set items
 *
 * \detail An entry in the set is valid if it is not equal to INVALID_ENTRY.
 *
 * An example to traverse the elements
 * \code
 * DynamicSet<> some_set;
 * ... // initialize some_set
 *
 * const IndexType N = some_set.size()
 * for(IndexType i=0; i< N; ++i)
 * {
 *   if( some_set.isValidEntry(i) )
 *   {
 *     ElementType el = some_set[i];
 *
 *     ... // do something with el
 *
 *   } // END if the entry is valid
 * } //END for all set entries
 *
 * \endcode
 */

template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType,
          typename SizePolicy = policies::DynamicRuntimeSize<PosType>>
class DynamicSet : public Set<PosType, ElemType>, SizePolicy
{
public:
  using PositionType = PosType;
  using ElementType = ElemType;
  using SetVectorType = std::vector<ElementType>;
  using SizePolicyType = SizePolicy;

  /// value to mark indices of deleted elements
  static constexpr ElementType INVALID_ENTRY = ~0;

  // predeclare SetBuilder  struct
  struct SetBuilder;

  // types for DynamicSet iterator
  template <typename T, bool C>
  class DynamicSetIterator;

  using const_iterator = DynamicSetIterator<ElementType, true>;
  using const_iterator_pair = std::pair<const_iterator, const_iterator>;

  using iterator = DynamicSetIterator<ElementType, false>;
  using iterator_pair = std::pair<iterator, iterator>;

public:
  /**
   * \brief Constructor for a DynamicSet
   * \param size The initial size of the set
   *
   * \note The set entries will be initialized such that set[i] = i
   */
  DynamicSet(PositionType size = SizePolicyType::DEFAULT_VALUE) : SizePolicy(size)
  {
    fill_array_default(size);
  };

  /**
   * \brief Constructor for a DynamicSet from a SetBuilder
   *
   * \note The set entries will be initialized such that set[i] = i
   */
  DynamicSet(const SetBuilder& builder) : SizePolicy(builder.m_size)
  {
    fill_array_default(builder.m_size.size());
  }

  //~DynamicSet();

public:
  /**
   * \class SetBuilder
   * \brief Helper class for constructing a DynamicSet.
   */
  struct SetBuilder
  {
    friend class DynamicSet;

    /** \brief Set the size of the DynamicSet using SizePolicy */
    SetBuilder& size(PositionType sz)
    {
      m_size = SizePolicy(sz);
      return *this;
    }

  private:
    SizePolicy m_size;
  };

  /**
   * \class DynamicSetIterator
   * \brief An stl-compliant random iterator type for a DynamicSet
   *
   * Uses the set's policies for efficient iteration
   * \tparam T The result type of the iteration
   * \tparam Const Boolean to indicate if this is a const iterator
   *
   * \note Most operators are implemented via the \a IteratorBase class
   *
   * \note Use of a const template parameter with conditional member and pointer
   * operations based on ideas from https://stackoverflow.com/a/49425072
   */
  template <typename T, bool Const>
  class DynamicSetIterator : public IteratorBase<DynamicSetIterator<T, Const>, PositionType>
  {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = PositionType;
    using reference = typename std::conditional<Const, const T&, T&>::type;
    using pointer = typename std::conditional<Const, const T*, T*>::type;
    using DynamicSetType = typename std::conditional<Const, const DynamicSet, DynamicSet>::type;
    using IterBase = IteratorBase<DynamicSetIterator<T, Const>, PositionType>;
    using IterBase::m_pos;

  public:
    /// \name Constructors, copying and assignment
    /// \{
    DynamicSetIterator() = default;

    DynamicSetIterator(PositionType pos, DynamicSetType& dSet) : IterBase(pos), m_dynamicSet(&dSet)
    { }

    /// \}

    /// \name Member and pointer operators
    /// \note We use the \a enable_if construct to implement both
    /// const and non-const iterators in the same implementation.
    /// \{

    /// Indirection operator for non-const iterator
    template <bool _Const = Const>
    typename std::enable_if<!_Const, reference>::type operator*()
    {
      return (*m_dynamicSet)[m_pos];
    }

    /// Indirection operator for const iterator
    template <bool _Const = Const>
    typename std::enable_if<_Const, reference>::type operator*() const
    {
      return (*m_dynamicSet)[m_pos];
    }

    /// Structure dereference operator for non-const iterator
    template <bool _Const = Const>
    typename std::enable_if<!_Const, pointer>::type operator->()
    {
      return &((*m_dynamicSet)[m_pos]);
    }

    /// Structure dereference operator for const iterator
    template <bool _Const = Const>
    typename std::enable_if<_Const, pointer>::type operator->() const
    {
      return &((*m_dynamicSet)[m_pos]);
    }

    /// Subscript operator for non-const iterator
    template <bool _Const = Const>
    typename std::enable_if<!_Const, reference>::type operator[](PositionType n)
    {
      return *(*this + n);
    }

    /// Subscript operator for const iterator
    template <bool _Const = Const>
    typename std::enable_if<_Const, reference>::type operator[](PositionType n) const
    {
      return *(*this + n);
    }

    /// \}

    /// \name Conversion operators
    /// \{

    /// Convert from iterator type to const_iterator type
    template <typename U>
    operator DynamicSetIterator<U, true>() const
    {
      return DynamicSetIterator<U, true>(this->m_pos, *this->m_dynamicSet);
    }
    /// \}

    PositionType index() const { return m_pos; }

    bool isValidEntry() const { return m_dynamicSet->isValidEntry(m_pos); };

  protected:
    /** Implementation of advance() as required by IteratorBase */
    void advance(PositionType n) { m_pos += n; }

  private:
    DynamicSetType* m_dynamicSet {nullptr};
  };

public:  // Functions related to iteration
  iterator begin() { return iterator(0, *this); }

  const_iterator begin() const { return const_iterator(0, *this); }

  const_iterator cbegin() const { return const_iterator(0, *this); }

  iterator end() { return iterator(size(), *this); }

  const_iterator end() const { return const_iterator(size(), *this); }

  const_iterator cend() const { return const_iterator(size(), *this); }

public:
  /// \name DynamicSet element access functions
  /// @{

  /**
   * \brief Access the element at position \a pos
   *
   * \pre pos must be between 0 and size()
   */
  ElementType at(PositionType pos) const { return operator[](pos); };

  /**
   * \brief Access the element at position \a pos
   *
   * \pre pos must be between 0 and size()
   */
  const ElementType& operator[](IndexType pos) const
  {
    verifyPositionImpl(pos);
    return m_data[pos];
  };

  /**
   * \brief Access the element at position \a pos
   *
   * \pre pos must be between 0 and size()
   */
  ElementType& operator[](IndexType pos)
  {
    verifyPositionImpl(pos);
    return m_data[pos];
  };

  /** \brief Returns a reference to the underlying set data */
  SetVectorType& data() { return m_data; }

  /** \brief Returns a const reference to the underlying set data */
  const SetVectorType& data() const { return m_data; }

  /**
   * \brief Given a value, find the index of the first entry containing it
   *
   * \return The index of the first element with value \a e, or INVALID_ENTRY
   * if none can be found.
   * \note This is an O(n) operation in the size of the set
   */
  IndexType findIndex(ElementType e) const
  {
    const int sz = size();
    for(int i = 0; i < sz; ++i)
    {
      if(m_data[i] == e)
      {
        return i;
      }
    }
    return INVALID_ENTRY;
  };

  /**
   * \brief Checks whether an element exists within the DynamicSet
   *
   * \return \a true if the set contains element with value \a e, false otherwise
   * \note This is an O(n) operation in the size of the set
   */
  bool contains(ElementType e) const
  {
    const int sz = size();
    for(int i = 0; i < sz; ++i)
    {
      if(m_data[i] == e)
      {
        return true;
      }
    }
    return false;
  };

  /// @}

public:
  /// \name Functions that deal with the set cardinality and indexing
  /// @{

  /**
   * \brief Returns the number of possible elements in the set
   *
   * \note Not all elements are necessarily valid since some elements
   * could have been deleted
   * \sa numberOfValidEntries(), isValidEntry()
   */
  AXOM_HOST_DEVICE inline PositionType size() const { return SizePolicy::size(); };

  /// \brief Uses \a SizePolicy::empty() to determine if the set is empty
  AXOM_HOST_DEVICE bool empty() const { return SizePolicy::empty(); };

  /// \brief Returns a positionset over the set elements
  PositionSet<PositionType> positions() const { return PositionSet<PositionType>(size()); }

  /**
   * \brief Return the number of valid entries in the set.
   *
   * \details This is an O(n) operation, because the class makes no assumption
   * that data was not changed by the user
   */
  PositionType numberOfValidEntries() const
  {
    PositionType nvalid = 0;

    const int sz = size();
    for(int i = 0; i < sz; ++i)
    {
      nvalid += (m_data[i] != INVALID_ENTRY) ? 1 : 0;
    }
    return nvalid;
  }

  /// \brief Returns true if this set is a subset of another set
  bool isSubset() const { return false; };

  /// @}

public:
  /// \name Functions that deal with validity checks
  /// @{

  /**
   * \brief Predicate to check if the entry at index \a i is valid
   *
   * The entry is valid when 0 <= i < size() and the value at index
   * \a i is not marked as \a INVALID_ENTRY
   */
  inline bool isValidEntry(IndexType i) const
  {
    return i >= 0 && i < size() && m_data[i] != INVALID_ENTRY;
  };

  /**
   * \brief Returns true if the DynamicSet instance is valid
   *
   * A DynamicSet is valid if each of its policies claim it to be valid.
   * This includes its \a SizePolicy, \a OffsetPolicy and \a StridePolicy
   */
  bool isValid(bool verboseOutput = false) const
  {
    bool bValid = SizePolicy::isValid(verboseOutput);
    return bValid;
  };

  /// @}

public:
  /// \name Functions that modify the set cardinality
  /// @{

  /**
   * \brief Insert an entry at the end of the set with value = ( size()-1 )
   */
  IndexType insert() { return insert(size()); }

  /**
   * \brief Insert an entry at the end of the set with the given value.
   * \param val the value of the inserted entry
   */
  IndexType insert(ElementType val)
  {
    m_data.push_back(val);
    SizePolicy::m_sz = m_data.size();
    return size() - 1;
  };

  /**
   * \brief Mark the corresponding entry as invalid
   *
   * \note It is not a problem to mark an INVALID_ENTRY as INVALID_ENTRY
   */
  void remove(IndexType idx)
  {
    verifyPositionImpl(idx);
    m_data[idx] = INVALID_ENTRY;
  };

  /**
   * \brief Resets to a default DynamicSet of size \a sz
   * 
   * \details The entry at index i will have the value of i for 0 <= i < sz
   */
  void reset(PositionType sz)
  {
    SizePolicy::m_sz = sz;
    fill_array_default(sz);
  }

  /// @}

private:
  /// \brief Debug check that the index \a pos is not out-of-range
  inline void verifyPosition(PositionType pos) const { verifyPositionImpl(pos); };

  /** \brief Debug check that the index \a pos is not out-of-range */
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(pos)) const
  {
    SLIC_ASSERT_MSG((pos >= 0) && (pos < size()),
                    "SLAM::DynamicSet -- requested out-of-range element at position "
                      << pos << ", but set only has " << size() << " elements.");
  };

  /** Fill each entry of the set such that its value is equal to its index. */
  void fill_array_default(PositionType sz)
  {
    if(sz < 0)
    {
      return;
    }

    m_data.resize(sz);
    for(int i = 0; i < sz; ++i)
    {
      m_data[i] = i;
    }
  }

private:
  SetVectorType m_data;
};

template <typename P, typename E, typename S>
constexpr typename DynamicSet<P, E, S>::ElementType DynamicSet<P, E, S>::INVALID_ENTRY;

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_DYNAMIC_SET_H_
