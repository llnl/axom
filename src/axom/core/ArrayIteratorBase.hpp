// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAYITERATORBASE_HPP_
#define AXOM_ARRAYITERATORBASE_HPP_

#include "axom/core/IteratorBase.hpp"  // for Iterator

namespace axom
{
/// \name ArrayIteratorBase to iterate through Array-like types
/// @{

/**
 * \class   ArrayIteratorBase
 * \brief   An iterator type for Array-like types.
 *          Each increment operation advances the iterator to the next
 *          element in the Array-like.
 * \tparam ArrayType The type of the array to iterate over
 * \tparam ValueType The type of the array's elements
 */
template <typename ArrayType, typename ValueType>
class ArrayIteratorBase : public IteratorBase<ArrayIteratorBase<ArrayType, ValueType>, IndexType>
{
private:
  using BaseType = IteratorBase<ArrayIteratorBase<ArrayType, ValueType>, IndexType>;

public:
  // Iterator traits required to satisfy LegacyRandomAccessIterator concept
  // before C++20
  // See: https://en.cppreference.com/w/cpp/iterator/iterator_traits
  using difference_type = IndexType;
  using value_type = typename std::remove_cv<ValueType>::type;
  using reference = ValueType&;
  using pointer = ValueType*;
  using iterator_category = std::random_access_iterator_tag;

public:
  using ArrayPointerType = ArrayType*;

  ArrayIteratorBase() : BaseType(0) { }

  AXOM_HOST_DEVICE
  ArrayIteratorBase(IndexType pos, ArrayPointerType arr) : BaseType(pos), m_arrayPtr(arr) { }

  /**
   * \brief Returns the current iterator value
   */
  AXOM_HOST_DEVICE
  ValueType& operator*() const { return m_arrayPtr->flatIndex(BaseType::m_pos); }

protected:
  /** Implementation of advance() as required by IteratorBase */
  AXOM_HOST_DEVICE
  void advance(IndexType n) { BaseType::m_pos += n; }

protected:
  ArrayPointerType m_arrayPtr {nullptr};
};  // end of ArrayIteratorBase class

/// @}

} /* namespace axom */

#endif /* AXOM_ARRAYITERATORBASE_HPP_ */
