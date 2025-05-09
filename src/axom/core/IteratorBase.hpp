// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IteratorBase.hpp
 *
 * \brief Contains iterator base classes
 */

#ifndef AXOM_ITERBASE_HPP_
#define AXOM_ITERBASE_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

namespace axom
{
/**
 * \class IteratorBase
 *
 * \brief Base class for a random access iterator over positions in a set
 *
 * This class is for keeping track of the position value in an iterator class.
 * It uses the Curiously Recurring Template Pattern (CRTP): The first template
 * parameter should be the actual iterator class. \n
 * The derived iterator class must implement `void advance(PositionType)`
 * to update the \a m_pos variable.
 * e.g.
 * \code{.cpp}
 * class IteratorDerived : public IteratorBase < IteratorDerived, PositionType >
 * {
 *   using IteratorBase<IteratorDerived, DataType>::m_pos;

 *   //...implementation of the accessing functions using m_pos...
 *
 *   //derive iterator class must implement this function, ideally protected.
 * protected:
 *   void advance( PositionType n) {...}
 * };
 * \endcode
 *
 * \tparam  IterType The iterator class to be the derived class
 * \tparam  PosType The type used to index into the iterator
 *          Must be a signed integral type
 */

template <typename IterType, typename PosType>
class IteratorBase
{
  static_assert(std::is_signed<PosType>::value, "PosType must be signed");
  static_assert(std::is_integral<PosType>::value, "PosType must be integral");

protected:
  IteratorBase() : m_pos(PosType()) { }

  AXOM_HOST_DEVICE
  explicit IteratorBase(PosType pos) : m_pos(pos) { }

private:
  /**
   *  Utility class to access protected function IterType::advance().
   *  Idea borrowed from: https://accu.org/index.php/journals/296
   */
  struct accessor : IterType
  {
    AXOM_HOST_DEVICE accessor(const IterType& base) : IterType(base) { }

    AXOM_SUPPRESS_HD_WARN
    AXOM_HOST_DEVICE
    static void adv(IterType& instance, PosType n)
    {
      // Protected member functions may only be accessed from a derived class
      // via an instance of the derived class.
      // As a workaround, we construct an instance of accessor, call
      // accessor::advance(), then slice it back to the base instance.
      accessor derived(instance);
      derived.advance(n);
      instance = derived;
    }
  };

private:
  /// Call the derived iterator's advance() function
  AXOM_HOST_DEVICE
  void adv(IterType& derived, PosType n) const { accessor::adv(derived, n); }

public:
  /// \name Equality and relational operators
  /// \{

  using iterator = IteratorBase<IterType, PosType>;

  /// Equality operator
  friend bool operator==(const iterator& lhs, const iterator& rhs)
  {
    return lhs.m_pos == rhs.m_pos;
  }
  /// Inequality operator
  AXOM_HOST_DEVICE
  friend bool operator!=(const iterator& lhs, const iterator& rhs)
  {
    return lhs.m_pos != rhs.m_pos;
  }

  /// Less than operator
  friend bool operator<(const iterator& lhs, const iterator& rhs) { return lhs.m_pos < rhs.m_pos; }
  /// Less than or equal operator
  friend bool operator<=(const iterator& lhs, const iterator& rhs)
  {
    return lhs.m_pos <= rhs.m_pos;
  }
  /// Greater than operator
  friend bool operator>(const iterator& lhs, const iterator& rhs) { return lhs.m_pos > rhs.m_pos; }
  /// Greater than or equal operator
  friend bool operator>=(const iterator& lhs, const iterator& rhs)
  {
    return lhs.m_pos >= rhs.m_pos;
  }
  /// \}

  /// \name Iterator advance and distance operators
  /// \{

  /// Pre-increment operator
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  IterType& operator++()
  {
    adv(getIter(), 1);
    return getIter();
  }
  /// Post-increment operator
  IterType operator++(int)
  {
    IterType ret = getIter();
    adv(getIter(), 1);
    return ret;
  }
  /// Pre-decrement operator
  IterType& operator--()
  {
    adv(getIter(), -1);
    return getIter();
  }
  /// Post-decrement operator
  IterType operator--(int)
  {
    IterType ret = getIter();
    adv(getIter(), -1);
    return ret;
  }

  /// Addition-assignment operator
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  IterType& operator+=(PosType n)
  {
    adv(getIter(), n);
    return getIter();
  }
  /// Subtraction-assignment operator
  IterType& operator-=(PosType n)
  {
    adv(getIter(), -n);
    return getIter();
  }

  /// Addition operator with iterator on left and position on right
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  friend IterType operator+(const IterType& it, PosType n)
  {
    IterType ret(it);
    ret.adv(ret, n);
    return ret;
  }
  /// Addition operator with position on left and iterator on right
  friend IterType operator+(PosType n, const IterType& it) { return operator+(it, n); }

  /// Subtraction operator with iterator on left and position on right
  friend IterType operator-(const IterType& it, PosType n) { return operator+(it, -n); }
  /// Subtraction operator with position on left and iterator on right
  friend IterType operator-(PosType n, const IterType& it) { return operator+(it, -n); }

  /// Difference operator
  friend PosType operator-(const IterType& a, const IterType& b) { return (a.m_pos - b.m_pos); }
  /// \}

private:
  /// Accessor to derived class
  AXOM_HOST_DEVICE
  IterType& getIter() { return *static_cast<IterType*>(this); }

  /// Const accessor to derived class
  AXOM_HOST_DEVICE
  const IterType& getIter() const { return *static_cast<const IterType*>(this); }

protected:
  PosType m_pos;
};

}  // end namespace axom

#endif  //  AXOM_ITERBASE_HPP_
