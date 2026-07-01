// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file ValuePolicies.hpp
 *
 * \brief Unified storage core for SLAM's scalar value policies.
 *
 * Slam's Size, Stride and Offset policies each historically defined
 * a near-identical `Runtime*` / `CompileTime*` pair, with the same: 
 *  (a) single-integer storage
 *  (b) defaulted/asserting constructors, and 
 *  (c) `operator()` accessors, 
 * and differing only in 
 *  (a) the name of the named accessor (`size()` / `stride()` / `offset()`), 
 *  (b) the validity predicate, and
 *  (c) the default value. 
 *
 * This file factors that shared substrate into one `RuntimeValue<Tag>` / `CompileTimeValue<auto, Tag>` family.
 *
 * A `Tag` type supplies the policy-specific knobs as static members:
 *   - `static constexpr IntType defaultValue();` -- the DEFAULT_VALUE
 *   - `static constexpr bool isValidValue(IntType);` -- validity predicate
 *
 * The named accessors (`size()`, `stride()`, `offset()`) are *not* provided here; 
 * they live in the thin leaf policies in SizePolicies.hpp / StridePolicies.hpp / OffsetPolicies.hpp
 * as one-line forwarders to `value()`, which keeps every existing call site spelling 
 * and signature unchanged while removing the storage/ctor/validity duplication.
 *
 * \note Multi-dimensional and dynamically-resizable policies (MultiDimStride,
 *  DynamicRuntimeSize) and the always-zero policies (ZeroSize) are not part of
 *  this scalar substrate and remain defined alongside their families.
 */

#ifndef SLAM_POLICIES_VALUE_H_
#define SLAM_POLICIES_VALUE_H_

#include "axom/core/Macros.hpp"

namespace axom::slam::policies
{
/// \name Value policy tags
/// \brief Tag types selecting the named-accessor family and the policy knobs
///  (default value, validity predicate) for the unified value-policy core.
/// \{

/*!
 * \brief Tag for set-size value policies.
 * \note Sizes may not be negative; the default size is zero.
 */
template <typename IntType>
struct SizeTag
{
  AXOM_HOST_DEVICE static constexpr IntType defaultValue() { return IntType {}; }
  AXOM_HOST_DEVICE static constexpr bool isValidValue(IntType v) { return v >= IntType {}; }
};

/*!
 * \brief Tag for set-stride value policies.
 * \note All non-zero strides are valid; the default stride is one.
 */
template <typename IntType>
struct StrideTag
{
  AXOM_HOST_DEVICE static constexpr IntType defaultValue() { return IntType(1); }
  AXOM_HOST_DEVICE static constexpr bool isValidValue(IntType v) { return v != IntType {}; }
};

/*!
 * \brief Tag for set-offset value policies.
 * \note Every offset is valid; the default offset is zero.
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
 * \brief Shared storage core for a runtime-settable scalar value policy.
 *
 * Stores a single \a IntType whose default is supplied by \a Tag.
 * Provides the generic `value()` accessors (const and mutable), `operator()`,
 * and an `isValid()` delegating to the tag's predicate.
 * 
 * Leaf policies derive from this and add their named accessor (`size()` / `stride()` / `offset()`).
 *
 * \tparam Tag a value-policy tag (SizeTag / StrideTag / OffsetTag)
 *  carrying the IntType, default value and validity predicate.
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
 * \brief Shared core for a compile-time-known scalar value policy.
 *
 * The value \a V is fixed at compile time.
 * The (defaulted) constructor argument exists only to satisfy 
 * the uniform policy-construction interface and is asserted to match \a V. 
 *
 * Provides the generic `value()` accessor, `operator()`, and an `isValid()` delegating to the tag's predicate.
 * Leaf policies derive from this and add their named accessor.
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

#endif  // SLAM_POLICIES_VALUE_H_
