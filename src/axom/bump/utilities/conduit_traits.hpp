// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_CONDUIT_TRAITS_HPP_
#define AXOM_BUMP_CONDUIT_TRAITS_HPP_

#include "axom/export/bump.h"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
namespace utilities
{
//------------------------------------------------------------------------------
/*!
 * \brief This class provides type traits that let us map C++ types
 *        to types / values useful in Conduit.
 */
template <typename T>
struct cpp2conduit
{ };

template <>
struct cpp2conduit<conduit::int8>
{
  using type = conduit::int8;
  static constexpr conduit::index_t id = conduit::DataType::INT8_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::int16>
{
  using type = conduit::int16;
  static constexpr conduit::index_t id = conduit::DataType::INT16_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::int32>
{
  using type = conduit::int32;
  static constexpr conduit::index_t id = conduit::DataType::INT32_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::int64>
{
  using type = conduit::int64;
  static constexpr conduit::index_t id = conduit::DataType::INT64_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::uint8>
{
  using type = conduit::uint8;
  static constexpr conduit::index_t id = conduit::DataType::UINT8_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::uint16>
{
  using type = conduit::uint16;
  static constexpr conduit::index_t id = conduit::DataType::UINT16_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::uint32>
{
  using type = conduit::uint32;
  static constexpr conduit::index_t id = conduit::DataType::UINT32_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::uint64>
{
  using type = conduit::uint64;
  static constexpr conduit::index_t id = conduit::DataType::UINT64_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::float32>
{
  using type = conduit::float32;
  static constexpr conduit::index_t id = conduit::DataType::FLOAT32_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::float64>
{
  using type = conduit::float64;
  static constexpr conduit::index_t id = conduit::DataType::FLOAT64_ID;
  AXOM_BUMP_EXPORT static const char *name;
};

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom

#endif
