// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Optional.hpp
 *
 * \brief A minimal, host-device "maybe a value" type for slam.
 *
 * std::optional is a host-only facility in slam's portability model
 * and is not guaranteed to be available/usable in device code 
 * across all of slam's backends (SEQ/OMP/CUDA/HIP).
 *
 * slam::Optional<T> is a trivially-structured aggregate of {engaged flag, storage} 
 * that is AXOM_HOST_DEVICE throughout and has no throwing value() accessor
 * (querying an unengaged Optional in a kernel cannot throw, 
 *  so the contract is "check has_value() first"; in debug host builds a violation asserts).
 *
 * It has sufficient functionality for a kernel-side optional (i.e to check if something is found)
 * but is not a drop-in std::optional -- it does not have: exceptions, monadic and_then/transforms, 
 * or in-place construction machinery.
 */

#ifndef SLAM_OPTIONAL_H_
#define SLAM_OPTIONAL_H_

#include "axom/core/Macros.hpp"
#include "axom/slam/policies/ConstexprAssert.hpp"

#include <type_traits>

namespace axom::slam
{
/*!
 * \class Optional
 * \brief Host-device "maybe a value of type T".
 *
 * \tparam T the (literal/trivially-copyable) value type. Intended for the
 *  small value and index types that flow through slam's device find APIs.
 */
template <typename T>
struct Optional
{
  T m_value {};
  bool m_engaged {false};

  /// \brief Construct a disengaged Optional (no value).
  AXOM_HOST_DEVICE constexpr Optional() = default;

  /// \brief Construct an engaged Optional holding \a value.
  AXOM_HOST_DEVICE constexpr Optional(const T& value) : m_value(value), m_engaged(true) { }

  /// \brief Whether this Optional holds a value.
  AXOM_HOST_DEVICE constexpr bool has_value() const { return m_engaged; }

  /// \brief Whether this Optional holds a value (bool conversion).
  AXOM_HOST_DEVICE constexpr explicit operator bool() const { return m_engaged; }

  /*!
   * \brief Access the contained value.
   * \pre has_value() is true. There is no throwing accessor -- in host debug builds,
   *  a disengaged access asserts, and during constant evaluation it is a compile error.
   */
  AXOM_HOST_DEVICE constexpr const T& operator*() const
  {
    SLAM_CONSTEXPR_ASSERT(m_engaged);
    return m_value;
  }
  AXOM_HOST_DEVICE constexpr T& operator*()
  {
    SLAM_CONSTEXPR_ASSERT(m_engaged);
    return m_value;
  }

  /// \brief Return the contained value if engaged, otherwise \a fallback.
  AXOM_HOST_DEVICE constexpr T value_or(const T& fallback) const
  {
    return m_engaged ? m_value : fallback;
  }
};

}  // end namespace axom::slam

#endif  // SLAM_OPTIONAL_H_
