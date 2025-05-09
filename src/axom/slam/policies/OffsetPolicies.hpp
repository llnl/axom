// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file OffsetPolicies.hpp
 *
 * \brief Offset policies for SLAM
 *
 * Offset policies are meant to represent the offset to the first element of
 * SLAM ordered set.
 *
 * A valid offset policy must support the following interface:
 *   * [required]
 *     * DEFAULT_VALUE is a public static constant of type IntType
 *     * offset() : IntType  -- returns the offset
 *     * isValid() : bool -- indicates whether the Offset policy of the set is
 *       valid  [optional]
 *     * operator(): IntType -- alternate accessor for the offset value
 */

#ifndef SLAM_POLICIES_OFFSET_H_
#define SLAM_POLICIES_OFFSET_H_

#include "axom/core/Macros.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
/**
 * \name OrderedSet_Offset_Policies
 * \brief A few default policies for the offset of an OrderedSet
 */

/// \{

/**
 * \brief A policy class for the offset in a set.  The offset can be set at
 * runtime.
 */
template <typename IntType>
struct RuntimeOffset
{
public:
  static const IntType DEFAULT_VALUE;

  AXOM_HOST_DEVICE RuntimeOffset(IntType off = DEFAULT_VALUE) : m_off(off) { }

  AXOM_HOST_DEVICE inline IntType offset() const { return m_off; }
  AXOM_HOST_DEVICE inline IntType& offset() { return m_off; }

  inline IntType operator()() const { return offset(); }
  inline IntType& operator()() { return offset(); }

  inline bool isValid(bool) const { return true; }

private:
  IntType m_off;
};

template <typename IntType>
const IntType RuntimeOffset<IntType>::DEFAULT_VALUE = IntType {};

/**
 * \brief A policy class for a compile-time known set offset
 */
template <typename IntType, IntType INT_VAL>
struct CompileTimeOffset
{
  static constexpr IntType DEFAULT_VALUE = INT_VAL;

  AXOM_HOST_DEVICE CompileTimeOffset(IntType val = DEFAULT_VALUE)
  {
    AXOM_UNUSED_VAR(val);
    SLIC_ASSERT_MSG(val == INT_VAL,
                    "slam::CompileTimeOffset -- tried to initialize a compile time "
                      << "offset with value (" << val << " ) that differs from "
                      << "the template parameter of " << INT_VAL << ".");
  }

  AXOM_HOST_DEVICE inline IntType offset() const { return INT_VAL; }

  inline IntType operator()() const { return offset(); }

  inline bool isValid(bool) const { return true; }
};

/**
 * \brief A policy class for when we have no offset
 */
template <typename IntType>
using ZeroOffset = CompileTimeOffset<IntType, 0>;

/// \}

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_OFFSET_H_
