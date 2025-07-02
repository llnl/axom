// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_DeviceHash_Hpp
#define Axom_Core_DeviceHash_Hpp

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

#include <type_traits>

namespace axom
{
namespace detail
{
template <typename T, typename Enable = void>
struct DeviceHashHelper;

/// \brief Specialization for integral types
template <typename T>
struct DeviceHashHelper<T, std::enable_if_t<std::is_integral<T>::value>>
{
  using argument_type = T;
  using result_type = axom::IndexType;
  AXOM_HOST_DEVICE axom::IndexType operator()(T value) const { return value; }
};

/// \brief Specialization for floating-point types
template <typename T>
struct DeviceHashHelper<T, std::enable_if_t<std::is_floating_point<T>::value>>
{
  using argument_type = T;
  using result_type = axom::IndexType;
  AXOM_HOST_DEVICE axom::IndexType operator()(T value) const
  {
    // Special case: -0.0 and 0.0 compare equal but have different byte representations.
    if(value == T {0.})
    {
      return 0;
    }
    return value;
  }
};

/// \brief SFINAE specialization for enum types
template <typename T>
struct DeviceHashHelper<T, std::enable_if_t<std::is_enum<T>::value>>
{
  using argument_type = T;
  using result_type = axom::IndexType;
  AXOM_HOST_DEVICE axom::IndexType operator()(T value) const
  {
    return static_cast<axom::IndexType>(value);
  }
};

/// \brief Specialization for pointer types
template <typename T>
struct DeviceHashHelper<T*, void>
{
  using argument_type = T*;
  using result_type = axom::IndexType;
  AXOM_HOST_DEVICE axom::IndexType operator()(T* ptr) const
  {
    return static_cast<axom::IndexType>(reinterpret_cast<std::uintptr_t>(ptr));
  }
};

/// \brief Default catch-all specialization. Passes through to std::hash.
template <typename T, typename Enable>
struct DeviceHashHelper
{
  using argument_type = T;
  using result_type = axom::IndexType;
  axom::IndexType operator()(const T& object) const
  {
    return static_cast<axom::IndexType>(std::hash<T> {}(object));
  }
};

}  // namespace detail

/*!
 * \class DeviceHash
 *
 * \brief Implements a host/device-callable hash function for supported types,
 *  and passes through to std::hash otherwise.
 */
template <typename T>
struct DeviceHash : public detail::DeviceHashHelper<T>
{
  using typename detail::DeviceHashHelper<T>::argument_type;
  using typename detail::DeviceHashHelper<T>::result_type;
};

}  // namespace axom

#endif
