// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file StridePolicies.hpp
 *
 * \brief Stride policies for SLAM
 *
 * Stride policies are meant to represent the fixed distance between consecutive
 * elements of an OrderedSet
 * A valid stride policy must support the following interface:
 *   [required]
 *    - DEFAULT_VALUE is a public static const IntType
 *    - IS_COMPILE_TIME is a public static const bool
 *    - stride() : IntType  -- returns the stride
 *    - isValid() : bool -- indicates whether the Stride policy of the set is valid
 *   [optional]
 *    - operator(): IntType -- alternate accessor for the stride value
 *
 * \note All non-zero stride values are valid.
 *
 * \note The single-stride Runtime/CompileTime storage, constructors and validity checking
 *  are provided by the unified RuntimeValue/CompileTimeValue core in ValuePolicies.hpp.
 *  The scalar policies below add only the stride-specific surface (named `stride()`/`shape()` accessors, 
 *  the dimensional typedefs, DEFAULT_VALUE/IS_COMPILE_TIME).
 *  MultiDimStride is a separate, inherently multi-dimensional policy and is unaffected.
 */

#pragma once

#include "axom/core/Macros.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/slam/policies/ValuePolicies.hpp"

namespace axom::slam::policies
{
/**
 * \name OrderedSet_Stride_Policies
 * \brief A few default policies for the stride of an OrderedSet
 */

/// \{

/**
 * \brief A policy class for the stride in a set.
 * When using this class, the stride can be set at runtime.
 */
template <typename IntType>
struct RuntimeStride : RuntimeValue<StrideTag<IntType>>
{
private:
  using BaseType = RuntimeValue<StrideTag<IntType>>;

public:
  static const IntType DEFAULT_VALUE;
  static const bool IS_COMPILE_TIME = false;
  constexpr static int NumDims = 1;

  using IndexType = IntType;
  using ShapeType = IntType;

  static constexpr IntType DefaultSize() { return StrideTag<IntType>::defaultValue(); }

  using BaseType::BaseType;

  /// \brief Returns the stride between consecutive elements.
  AXOM_HOST_DEVICE constexpr IntType stride() const { return this->value(); }
  AXOM_HOST_DEVICE constexpr IntType& stride() { return this->value(); }

  /*!
   * \brief Returns the shape of the inner data for a given stride.
   *  This only has meaning when used with Map-based types.
   */
  AXOM_HOST_DEVICE constexpr IntType shape() const { return this->value(); }

  void setStride(IntType str) { this->m_value = str; }
};

template <typename IntType>
const IntType RuntimeStride<IntType>::DEFAULT_VALUE = StrideTag<IntType>::defaultValue();

/// \brief A policy class for a compile-time known stride
template <typename IntType, IntType INT_VAL>
struct CompileTimeStride : CompileTimeValue<INT_VAL, StrideTag<IntType>>
{
private:
  using BaseType = CompileTimeValue<INT_VAL, StrideTag<IntType>>;

public:
  static const IntType DEFAULT_VALUE = INT_VAL;
  static const bool IS_COMPILE_TIME = true;
  constexpr static int NumDims = 1;

  using IndexType = IntType;
  using ShapeType = IntType;

  static constexpr IntType DefaultSize() { return DEFAULT_VALUE; }

  using BaseType::BaseType;

  AXOM_HOST_DEVICE constexpr IntType stride() const { return INT_VAL; }
  AXOM_HOST_DEVICE constexpr IntType shape() const { return INT_VAL; }

  AXOM_HOST_DEVICE void setStride(IntType AXOM_DEBUG_PARAM(val))
  {
    SLIC_ASSERT_MSG(val == INT_VAL,
                    "slam::CompileTimeStride -- tried to set a compile time stride"
                      << " with value (" << val << " ) that differs from the template"
                      << " parameter of " << INT_VAL << ".");
  }
};

/// \brief A policy class for a set with stride one (i.e. the default stride)
template <typename IntType>
using StrideOne = CompileTimeStride<IntType, 1>;

/// \brief A policy class for a set with multi-dimensional stride. Assumed layout is row-major.
template <typename IntType, int Dims>
struct MultiDimStride
{
  using IndexType = IntType;
  using ShapeType = StackArray<IntType, Dims>;
  constexpr static int NumDims = Dims;

  static ShapeType DefaultSize()
  {
    ShapeType array;
    for(int i = 0; i < Dims; i++)
    {
      array[i] = 1;
    }
    return array;
  }

  AXOM_HOST_DEVICE MultiDimStride(StackArray<IntType, Dims> shape) : m_shape(shape)
  {
    m_strides[Dims - 1] = 1;
    for(int i = Dims - 2; i >= 0; i--)
    {
      m_strides[i] = m_strides[i + 1] * m_shape[i + 1];
    }
  }

  /// \brief Returns the "flat" stride of all the subcomponents.
  AXOM_HOST_DEVICE inline IntType stride() const { return m_shape[0] * m_strides[0]; }

  inline IntType operator()() const { return stride(); }
  inline IntType& operator()() { return stride(); }

  /// \brief Returns the strides for each indexing dimension.
  AXOM_HOST_DEVICE inline ShapeType strides() const { return m_strides; }
  /*!
   * \brief Returns the multi-dimensional shape of the inner data.
   *  This only has meaning when used with Map-based types.
   */
  AXOM_HOST_DEVICE inline ShapeType shape() const { return m_shape; }

private:
  ShapeType m_shape;
  ShapeType m_strides;
};

/// \}

}  // end namespace axom::slam::policies
