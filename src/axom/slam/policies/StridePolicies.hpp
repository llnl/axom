// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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
 *   * [required]
 *   * DEFAULT_VALUE is a public static const IntType
 *   * IS_COMPILE_TIME is a public static const bool
 *   * stride() : IntType  -- returns the stride
 *   * isValid() : bool -- indicates whether the Stride policy of the set is valid
 *   * [optional]
 *   * operator(): IntType -- alternate accessor for the stride value
 *
 * \note All non-zero stride values are valid.
 */

#ifndef SLAM_POLICIES_STRIDE_H_
#define SLAM_POLICIES_STRIDE_H_

#include "axom/core/Macros.hpp"

namespace axom
{
namespace slam
{
namespace policies
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
struct RuntimeStride
{
public:
  static const IntType DEFAULT_VALUE = IntType(1);
  static const bool IS_COMPILE_TIME = false;
  constexpr static int NumDims = 1;

  using IndexType = IntType;
  using ShapeType = IntType;

  static constexpr IntType DefaultSize() { return DEFAULT_VALUE; }

  AXOM_HOST_DEVICE RuntimeStride(IntType stride = DEFAULT_VALUE) : m_stride(stride) { }

  /// \brief Returns the stride between consecutive elements.
  AXOM_HOST_DEVICE inline IntType stride() const { return m_stride; }
  /*!
   * \brief Returns the shape of the inner data for a given stride.
   *  This only has meaning when used with Map-based types.
   */
  AXOM_HOST_DEVICE inline IntType shape() const { return m_stride; }
  AXOM_HOST_DEVICE inline IntType& stride() { return m_stride; }

  void setStride(IntType str) { m_stride = str; }

  inline IntType operator()() const { return stride(); }
  inline IntType& operator()() { return stride(); }

  /** All non-zero strides are valid     */
  inline bool isValid(bool) const { return (m_stride != 0); }

  //inline bool hasStride() const       { return m_stride != IntType(); }

private:
  IntType m_stride;
};

/**
 * \brief A policy class for a compile-time known stride
 */
template <typename IntType, IntType INT_VAL>
struct CompileTimeStride
{
  static const IntType DEFAULT_VALUE = INT_VAL;
  static const bool IS_COMPILE_TIME = true;
  constexpr static int NumDims = 1;

  using IndexType = IntType;
  using ShapeType = IntType;

  static constexpr IntType DefaultSize() { return DEFAULT_VALUE; }

  AXOM_HOST_DEVICE CompileTimeStride(IntType val = DEFAULT_VALUE) { setStride(val); }

  AXOM_HOST_DEVICE inline IntType stride() const { return INT_VAL; }
  AXOM_HOST_DEVICE inline IntType shape() const { return INT_VAL; }
  inline IntType operator()() const { return stride(); }

  AXOM_HOST_DEVICE void setStride(IntType AXOM_DEBUG_PARAM(val))
  {
    SLIC_ASSERT_MSG(val == INT_VAL,
                    "slam::CompileTimeStride -- tried to set a compile time stride"
                      << " with value (" << val << " ) that differs from the template"
                      << " parameter of " << INT_VAL << ".");
  }

  /** All non-zero strides are valid     */
  inline bool isValid(bool) const { return (INT_VAL != 0); }
};

/**
 * \brief A policy class for a set with stride one (i.e. the default stride)
 */
template <typename IntType>
using StrideOne = CompileTimeStride<IntType, 1>;

/**
 * \brief A policy class for a set with multi-dimensional stride. Assumed
 *  layout is row-major.
 */
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

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_STRIDE_H_
