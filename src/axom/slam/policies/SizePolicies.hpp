// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file SizePolicies.hpp
 *
 * \brief Size policies for SLAM
 *
 * Size policies are meant to represent the size of a SLAM entity (e.g. the size
 * of a set). A valid size policy must support the following interface:
 *   * [required]
 *   * DEFAULT_VALUE is a public static constant of type IntType
 *   * size() : IntType  -- returns the underlying integer size
 *   * empty() : bool -- returns whether the size is zero
 *   * isValid() : bool -- indicates whether the Size policy of the set is
 *      valid
 *   * [optional]
 *     * operator(): IntType -- alternate accessor for the size value
 */

#ifndef SLAM_POLICIES_SIZE_H_
#define SLAM_POLICIES_SIZE_H_

#include "axom/core/Macros.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
/**
 * \name OrderedSet_Size_Policies
 * \brief A few default policies for the size of an OrderedSet
 */

/// \{

/**
 * \brief A policy class for the size of a set whose value can be set at
 * runtime.
 */
template <typename IntType>
struct RuntimeSize
{
public:
  static const IntType DEFAULT_VALUE;

  AXOM_HOST_DEVICE RuntimeSize(IntType sz = DEFAULT_VALUE) : m_sz(sz) { }

  AXOM_HOST_DEVICE inline IntType size() const { return m_sz; }
  AXOM_HOST_DEVICE inline IntType& size() { return m_sz; }

  inline IntType operator()() const { return size(); }
  inline IntType& operator()() { return size(); }

  AXOM_HOST_DEVICE inline bool empty() const { return m_sz == IntType(); }
  inline bool isValid(bool) const
  {
    // We do not (currently) allow negatively sized sets
    return m_sz >= IntType();
  }

protected:
  IntType m_sz;
};

template <typename IntType>
const IntType RuntimeSize<IntType>::DEFAULT_VALUE = IntType {};

/**
 * \brief A policy class for the size of a set that can be modified at runtime
 */
template <typename IntType>
struct DynamicRuntimeSize : public RuntimeSize<IntType>
{
public:
  AXOM_HOST_DEVICE DynamicRuntimeSize(IntType sz = RuntimeSize<IntType>::DEFAULT_VALUE)
    : RuntimeSize<IntType>(sz)
  { }

  inline void setSize(IntType s)
  {
    if(s >= 0)
    {
      RuntimeSize<IntType>::m_sz = s;
    }
  }

  inline void add(IntType s)
  {
    if(s >= 0)
    {
      RuntimeSize<IntType>::m_sz += s;
    }
  }

  inline void subtract(IntType s)
  {
    if(s >= 0)
    {
      RuntimeSize<IntType>::m_sz -= s;
    }
  }
};

/**
 * \brief A policy class for a compile-time known set size
 */
template <typename IntType, IntType INT_VAL>
struct CompileTimeSize
{
  static const IntType DEFAULT_VALUE = INT_VAL;

  AXOM_HOST_DEVICE CompileTimeSize(IntType val = INT_VAL)
  {
    AXOM_UNUSED_VAR(val);
    SLIC_ASSERT_MSG(val == INT_VAL,
                    "slam::CompileTimeSize -- tried to initialize a compile time size "
                      << "policy with value (" << val << " ) that differs from the "
                      << "template parameter of " << INT_VAL << ".");
  }

  AXOM_HOST_DEVICE inline IntType size() const { return INT_VAL; }

  inline IntType operator()() const { return size(); }

  AXOM_HOST_DEVICE inline bool empty() const { return INT_VAL == IntType {}; }
  inline bool isValid(bool) const
  {
    // We do not (currently) allow negatively sized sets
    return INT_VAL >= IntType {};
  }
};

/**
 * \brief A policy class for an empty set (no size)
 */
template <typename IntType>
struct ZeroSize
{
  static const IntType DEFAULT_VALUE;

  AXOM_HOST_DEVICE ZeroSize(IntType val = DEFAULT_VALUE)
  {
    AXOM_UNUSED_VAR(val);
    SLIC_ASSERT_MSG(val == DEFAULT_VALUE,
                    "slam::ZeroSize policy-- tried to initialize a NoSize set with "
                      << "value with value (" << val << " ) but should always be zero.");
  }

  AXOM_HOST_DEVICE inline IntType size() const { return DEFAULT_VALUE; }
  inline IntType operator()() const { return size(); }
  AXOM_HOST_DEVICE inline bool empty() const { return true; }
  inline bool isValid(bool) const { return true; }
};

template <typename IntType>
const IntType ZeroSize<IntType>::DEFAULT_VALUE = IntType {};

/// \}

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_SIZE_H_
