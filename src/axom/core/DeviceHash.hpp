// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_DeviceHash_Hpp
#define Axom_Core_DeviceHash_Hpp

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

#include <cstdint>
#include <cstring>
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
  using result_type = std::uint64_t;
  AXOM_HOST_DEVICE std::uint64_t operator()(T value) const
  {
    return static_cast<std::uint64_t>(value);
  }
};

/// \brief Specialization for floating-point types
template <typename T>
struct DeviceHashHelper<T, std::enable_if_t<std::is_floating_point<T>::value>>
{
  using argument_type = T;
  using result_type = std::uint64_t;
  AXOM_HOST_DEVICE std::uint64_t operator()(T value) const
  {
    // -0.0 and 0.0 compare equal but have different bit patterns; normalize so both hash identically
    if(value == T {0.})
    {
      value = T {0.};
    }

    // Hash the bit pattern, not the converted value.
    // A float-to-integer value conversion collapses every key sharing an integer part,
    // e.g. all numbers between -1 and 1 converts to integer 0

    // NUM_WORDS is 1 for float or double, possibly 2 for long double
    constexpr std::size_t NUM_WORDS = (sizeof(T) + sizeof(std::uint64_t) - 1) / sizeof(std::uint64_t);
    // zero out words since we might only copy 4 bytes in for floats
    std::uint64_t words[NUM_WORDS] = {0};
    memcpy(words, &value, sizeof(T));

    std::uint64_t result = words[0];
    // Extra processing fortypes wider than 64 bits (long double).
    // Use an odd multiplier (2^64/golden-ratio-phi),
    // so the halves cannot cancel under a later XOR-style mixer
    for(std::size_t i = 1; i < NUM_WORDS; i++)
    {
      result = result * std::uint64_t {0x9e3779b97f4a7c15} + words[i];
    }

    return result;
  }
};

/// \brief SFINAE specialization for enum types
template <typename T>
struct DeviceHashHelper<T, std::enable_if_t<std::is_enum<T>::value>>
{
  using argument_type = T;
  using result_type = std::uint64_t;
  AXOM_HOST_DEVICE std::uint64_t operator()(T value) const
  {
    return static_cast<std::uint64_t>(value);
  }
};

/// \brief Specialization for pointer types
template <typename T>
struct DeviceHashHelper<T*, void>
{
  using argument_type = T*;
  using result_type = std::uint64_t;
  AXOM_HOST_DEVICE std::uint64_t operator()(T* ptr) const
  {
    return static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(ptr));
  }
};

/// \brief Default catch-all specialization. Passes through to std::hash.
template <typename T, typename Enable>
struct DeviceHashHelper
{
  using argument_type = T;
  using result_type = std::uint64_t;
  std::uint64_t operator()(const T& object) const
  {
    return static_cast<std::uint64_t>(std::hash<T> {}(object));
  }
};

}  // namespace detail

/*!
 * \class DeviceHash
 *
 * \brief Implements a host/device-callable hash function for supported types,
 *  and passes through to std::hash otherwise.
 *
 *  The result type is always std::uint64_t, independent of the configured axom::IndexType width.
 *  Hashes feed bit mixers and bucket selection, where truncating wide keys (e.g. 64-bit Morton codes)
 *  to a 32-bit IndexType before mixing would make keys equal mod 2^32 collide.
 */
template <typename T>
struct DeviceHash : public detail::DeviceHashHelper<T>
{
  using typename detail::DeviceHashHelper<T>::argument_type;
  using typename detail::DeviceHashHelper<T>::result_type;
};

}  // namespace axom

#endif
