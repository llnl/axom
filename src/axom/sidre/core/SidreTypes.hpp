// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 *  \file SidreTypes.hpp
 *
 *  \brief File containing types used in the Sidre component of axom.
 *
 */

#ifndef SIDRE_TYPES_HPP_
#define SIDRE_TYPES_HPP_

#include "SidreDataTypeIds.h"
#include "conduit.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

namespace axom
{
namespace sidre
{
// Type aliases to make Conduit usage easier and less visible in the Sidre API

/*!
 * \brief The Node class is the primary object in Conduit.
 */
using Node = conduit::Node;

/*!
 * \brief DataType is a general Conduit descriptor.
 */
using DataType = conduit::DataType;

/*!
 * \brief A Conduit Schema describes the data in a Node
 */
using Schema = conduit::Schema;

/*!
 * \brief IndexType is used for any labeling of a sidre object by an
 *        integer identifier.
 */
using IndexType = axom::IndexType;

/*!
 * \brief Common invalid index identifier used in sidre.
 */
constexpr IndexType InvalidIndex = axom::InvalidIndex;

/*!
 * \brief Common invalid name string usde in sidre.
 */
const std::string InvalidName = axom::utilities::string::InvalidName;

/*!
 * \brief Returns true if idx is valid, else false.
 *
 * Used for the loop test when iterating over the
 * Buffers or Attributes in a DataStore or the Views or Groups in a Group.
 */
inline bool indexIsValid(IndexType idx) { return idx != InvalidIndex; }

/*!
 * \brief Returns true if name is valid, else false.
 */
inline bool nameIsValid(const std::string& name)
{
  return name != axom::utilities::string::InvalidName;
}

/*!
 * \brief Enum that holds the numeric data type id options for sidre types.
 */

enum DataTypeId
{
  NO_TYPE_ID = SIDRE_NO_TYPE_ID,
  INT8_ID = SIDRE_INT8_ID,
  INT16_ID = SIDRE_INT16_ID,
  INT32_ID = SIDRE_INT32_ID,
  INT64_ID = SIDRE_INT64_ID,

  UINT8_ID = SIDRE_UINT8_ID,
  UINT16_ID = SIDRE_UINT16_ID,
  UINT32_ID = SIDRE_UINT32_ID,
  UINT64_ID = SIDRE_UINT64_ID,

  FLOAT32_ID = SIDRE_FLOAT32_ID,
  FLOAT64_ID = SIDRE_FLOAT64_ID,

  CHAR8_STR_ID = SIDRE_CHAR8_STR_ID,

  INT_ID = SIDRE_INT_ID,
  UINT_ID = SIDRE_UINT_ID,
  LONG_ID = SIDRE_LONG_ID,
  ULONG_ID = SIDRE_ULONG_ID,
  FLOAT_ID = SIDRE_FLOAT_ID,
  DOUBLE_ID = SIDRE_DOUBLE_ID
};

/// @cond INCLUDE_DETAIL

/*!
 * \brief The detail namespace contains code that is either used internally by
 *  the sidre implementation or is under evaluation.
 */
namespace detail
{
/*!
 * \brief Type traits to assist in converting compiler types to the appropriate
 *  data type ids.
 */
template <typename T>
struct SidreTT
{
  static const DataTypeId id = NO_TYPE_ID;
};

template <>
struct SidreTT<std::int8_t>
{
  static const DataTypeId id = INT8_ID;
};
template <>
struct SidreTT<std::int16_t>
{
  static const DataTypeId id = INT16_ID;
};
template <>
struct SidreTT<std::int32_t>
{
  static const DataTypeId id = INT32_ID;
};
template <>
struct SidreTT<std::int64_t>
{
  static const DataTypeId id = INT64_ID;
};

template <>
struct SidreTT<std::uint8_t>
{
  static const DataTypeId id = UINT8_ID;
};
template <>
struct SidreTT<std::uint16_t>
{
  static const DataTypeId id = UINT16_ID;
};
template <>
struct SidreTT<std::uint32_t>
{
  static const DataTypeId id = UINT32_ID;
};
template <>
struct SidreTT<std::uint64_t>
{
  static const DataTypeId id = UINT64_ID;
};

template <>
struct SidreTT<axom::float32>
{
  static const DataTypeId id = FLOAT32_ID;
};
template <>
struct SidreTT<axom::float64>
{
  static const DataTypeId id = FLOAT64_ID;
};
} /* end namespace detail */
/// @endcond

/*!
 * \brief TypeID is used to identify the type of a buffer (SIDRE_INT8_ID, etc).
 */
using TypeID = DataTypeId;

/*!
 * \brief Convenience function to convert int to TypeID type.
 *
 *  Used to convert C defines to C++ enumerations.
 */
inline TypeID getTypeID(const int typeID) { return static_cast<TypeID>(typeID); }

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_TYPES_HPP_ */
