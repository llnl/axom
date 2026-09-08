// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file ValuePolicies.hpp
 *
 * \brief Shared scalar storage for size, stride and offset policies.
 *
 * RuntimeValue stores one integer. CompileTimeValue supplies a constant without
 * storing a value in each object. Both provide value(), operator() and isValid().
 *
 * A tag defines the default value and validity predicate through defaultValue()
 * and isValidValue(value). The size, stride and offset policies add their named
 * accessors in SizePolicies.hpp, StridePolicies.hpp and OffsetPolicies.hpp.
 *
 * MultiDimStride, DynamicRuntimeSize and ZeroSize use separate implementations.
 */

#pragma once

#include "axom/core/Macros.hpp"

namespace axom::slam::policies
{
/// \name Value policy tags
/// \brief Default values and validity predicates for scalar policies.
/// \{

/*!
 * \brief Tag for set-size value policies.
 * \note Sizes must be nonnegative. The default is zero.
 */
template <typename IntType>
struct SizeTag
{
  AXOM_HOST_DEVICE static constexpr IntType defaultValue() { return IntType {}; }
  AXOM_HOST_DEVICE static constexpr bool isValidValue(IntType v) { return v >= IntType {}; }
};

/*!
 * \brief Tag for set-stride value policies.
 * \note Strides must be nonzero. The default is one. Maps also require positive strides.
 */
template <typename IntType>
struct StrideTag
{
  AXOM_HOST_DEVICE static constexpr IntType defaultValue() { return IntType(1); }
  AXOM_HOST_DEVICE static constexpr bool isValidValue(IntType v) { return v != IntType {}; }
};

/*!
 * \brief Tag for set-offset value policies.
 * \note Every offset is valid. The default is zero.
 */
template <typename IntType>
struct OffsetTag
{
  AXOM_HOST_DEVICE static constexpr IntType defaultValue() { return IntType {}; }
  AXOM_HOST_DEVICE static constexpr bool isValidValue(IntType) { return true; }
};

/// \}

/*!
 * \class RuntimeValue
 *
 * \brief Store a scalar policy value that can change at runtime.
 *
 * Stores a single \a IntType whose default is supplied by \a Tag.
 * Provides const and mutable value() access, operator(),
 * and an `isValid()` delegating to the tag's predicate.
 * 
 * Derived policies add size(), stride() or offset().
 *
 * \tparam Tag Supplies defaultValue() and isValidValue(). The default's type is IntType.
 */
template <typename Tag>
struct RuntimeValue
{
public:
  using TagType = Tag;
  using IntType = decltype(Tag::defaultValue());

  AXOM_HOST_DEVICE constexpr RuntimeValue(IntType val = Tag::defaultValue()) : m_value(val) { }

  AXOM_HOST_DEVICE constexpr auto value() const { return m_value; }
  AXOM_HOST_DEVICE constexpr auto& value() { return m_value; }

  constexpr auto operator()() const { return value(); }
  constexpr auto& operator()() { return value(); }

  constexpr bool isValid(bool) const { return Tag::isValidValue(m_value); }

protected:
  IntType m_value;
};

/*!
 * \class CompileTimeValue
 *
 * \brief Supply a scalar policy value fixed at compile time.
 *
 * The value \a V is fixed at compile time.
 * The constructor accepts an argument so that callers can construct runtime and
 * compile-time policies the same way. The argument must equal \a V.
 *
 * Provides value(), operator() and isValid(). Derived policies add their named accessor.
 *
 * \tparam V the compile-time value (its type is the policy's IntType).
 * \tparam Tag a value-policy tag carrying the default value and validity predicate.
 */
template <auto V, typename Tag>
struct CompileTimeValue
{
public:
  using TagType = Tag;
  using IntType = decltype(V);

  static constexpr IntType VALUE = V;

  AXOM_HOST_DEVICE constexpr CompileTimeValue(IntType val = V)
  {
    AXOM_UNUSED_VAR(val);
    AXOM_CONSTEXPR_ASSERT(val == V);
  }

  AXOM_HOST_DEVICE constexpr IntType value() const { return V; }

  constexpr IntType operator()() const { return value(); }

  constexpr bool isValid(bool) const { return Tag::isValidValue(V); }
};

}  // end namespace axom::slam::policies
