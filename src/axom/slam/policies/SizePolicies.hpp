// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
 *   * isValid() : bool -- indicates whether the Size policy of the set is valid
 *   * [optional]
 *     * operator(): IntType -- alternate accessor for the size value
 *
 * \note The Runtime/CompileTime storage, constructors and validity checking
 *  are provided by the unified RuntimeValue/CompileTimeValue core in ValuePolicies.hpp.
 *  The scalar policies below add only the named `size()` accessor, `empty()`, and the DEFAULT_VALUE member.
 */

#ifndef SLAM_POLICIES_SIZE_H_
#define SLAM_POLICIES_SIZE_H_

#include "axom/core/Macros.hpp"
#include "axom/slic.hpp"
#include "axom/slam/policies/ValuePolicies.hpp"

namespace axom::slam::policies
{
/**
 * \name OrderedSet_Size_Policies
 * \brief A few default policies for the size of an OrderedSet
 */

/// \{

/// \brief A policy class for the size of a set whose value can be set at runtime
template <typename IntType>
struct RuntimeSize : RuntimeValue<SizeTag<IntType>>
{
private:
  using BaseType = RuntimeValue<SizeTag<IntType>>;

public:
  static const IntType DEFAULT_VALUE;

  using BaseType::BaseType;

  AXOM_HOST_DEVICE constexpr IntType size() const { return this->value(); }
  AXOM_HOST_DEVICE constexpr IntType& size() { return this->value(); }

  AXOM_HOST_DEVICE constexpr bool empty() const { return this->m_value == IntType(); }
};

template <typename IntType>
const IntType RuntimeSize<IntType>::DEFAULT_VALUE = SizeTag<IntType>::defaultValue();

/// \brief A policy class for the size of a set that can be modified at runtime
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
      RuntimeSize<IntType>::m_value = s;
    }
  }

  inline void add(IntType s)
  {
    if(s >= 0)
    {
      RuntimeSize<IntType>::m_value += s;
    }
  }

  inline void subtract(IntType s)
  {
    if(s >= 0)
    {
      RuntimeSize<IntType>::m_value -= s;
    }
  }
};

/// \brief A policy class for a compile-time known set size
template <typename IntType, IntType INT_VAL>
struct CompileTimeSize : CompileTimeValue<INT_VAL, SizeTag<IntType>>
{
private:
  using BaseType = CompileTimeValue<INT_VAL, SizeTag<IntType>>;

public:
  static const IntType DEFAULT_VALUE = INT_VAL;

  using BaseType::BaseType;

  AXOM_HOST_DEVICE constexpr IntType size() const { return INT_VAL; }

  AXOM_HOST_DEVICE constexpr bool empty() const { return INT_VAL == IntType {}; }
};

/// \brief A policy class for an empty set (no size)
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

}  // end namespace axom::slam::policies

#endif  // SLAM_POLICIES_SIZE_H_
