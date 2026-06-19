// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Optional.hpp
 *
 * \brief A minimal, host/device "maybe a value" type.
 *
 * `std::optional` is a host-only facility in Axom's portability model 
 * and is not guaranteed to be available/usable in device code across all backends (SEQ/OMP/CUDA/HIP).
 *
 * `axom::Optional<T>` is a trivially-structured aggregate of `{storage, engaged flag}`
 * that is `AXOM_HOST_DEVICE` throughout and has no throwing `value()` accessor.
 *
 * The contract is to check `has_value()` before accessing the value.
 * In host debug builds, a disengaged dereference asserts, and during constant evaluation it is a compile error.
 *
 * This type is intentionally small and is not a drop-in replacement for `std::optional`.
 * Specifically, it does not provide exceptions, monadic combinators, or in-place construction.
 */

#ifndef AXOM_OPTIONAL_HPP_
#define AXOM_OPTIONAL_HPP_

#include "axom/core/Macros.hpp"

#include <type_traits>

namespace axom
{
/*!
 * \class Optional
 * \brief Host/device "maybe a value of type T".
 *
 * \tparam T The value type. Intended for small, trivially-copyable types
 *  that are safe to capture and use in device kernels.
 */
template <typename T>
struct Optional
{
  static_assert(std::is_trivially_copyable_v<T>,
                "axom::Optional<T> is intended for trivially-copyable value types");

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
   * \pre has_value() is true. There is no throwing accessor. 
   *  In host debug builds, a disengaged access asserts,
   *  and during constant evaluation it is a compile error.
   */
  AXOM_HOST_DEVICE constexpr const T& operator*() const
  {
    AXOM_CONSTEXPR_ASSERT(m_engaged);
    return m_value;
  }
  AXOM_HOST_DEVICE constexpr T& operator*()
  {
    AXOM_CONSTEXPR_ASSERT(m_engaged);
    return m_value;
  }

  /// \brief Return the contained value if engaged, otherwise \a fallback.
  AXOM_HOST_DEVICE constexpr T value_or(const T& fallback) const
  {
    return m_engaged ? m_value : fallback;
  }
};

}  // namespace axom

#endif  // AXOM_OPTIONAL_HPP_
