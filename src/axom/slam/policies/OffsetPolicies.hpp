// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
 *
 * \note The Runtime/CompileTime storage, constructors and validity checking
 *  are provided by the unified RuntimeValue/CompileTimeValue core in ValuePolicies.hpp.
 *  The policies below add only the named `offset()` accessor and the DEFAULT_VALUE member.
 */

#ifndef SLAM_POLICIES_OFFSET_H_
#define SLAM_POLICIES_OFFSET_H_

#include "axom/core/Macros.hpp"
#include "axom/slam/policies/ValuePolicies.hpp"

namespace axom::slam::policies
{
/**
 * \name OrderedSet_Offset_Policies
 * \brief A few default policies for the offset of an OrderedSet
 */

/// \{

/// \brief A policy class for the offset in a set.  The offset can be set at runtime.
template <typename IntType>
struct RuntimeOffset : RuntimeValue<OffsetTag<IntType>>
{
private:
  using BaseType = RuntimeValue<OffsetTag<IntType>>;

public:
  static const IntType DEFAULT_VALUE;

  using BaseType::BaseType;

  AXOM_HOST_DEVICE inline IntType offset() const { return this->value(); }
  AXOM_HOST_DEVICE inline IntType& offset() { return this->value(); }
};

template <typename IntType>
const IntType RuntimeOffset<IntType>::DEFAULT_VALUE = OffsetTag<IntType>::defaultValue();

/// \brief A policy class for a compile-time known set offset
template <typename IntType, IntType INT_VAL>
struct CompileTimeOffset : CompileTimeValue<INT_VAL, OffsetTag<IntType>>
{
private:
  using BaseType = CompileTimeValue<INT_VAL, OffsetTag<IntType>>;

public:
  static constexpr IntType DEFAULT_VALUE = INT_VAL;

  using BaseType::BaseType;

  AXOM_HOST_DEVICE inline IntType offset() const { return INT_VAL; }
};

/// \brief A policy class for when we have no offset
template <typename IntType>
using ZeroOffset = CompileTimeOffset<IntType, 0>;

/// \}

}  // end namespace axom::slam::policies

#endif  // SLAM_POLICIES_OFFSET_H_
