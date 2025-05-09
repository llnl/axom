// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FIELDTYPES_HPP_
#define MINT_FIELDTYPES_HPP_

#include "axom/core/Types.hpp"  // for axom type definitions

namespace axom
{
namespace mint
{
/*!
 * \brief Enumerates the set of allowable field types.
 */
enum FieldType
{
  UNDEFINED_FIELD_TYPE = -1,  //!< UNDEFINED_FIELD_TYPE

  FLOAT_FIELD_TYPE,   //!< single-precision floating point field type
  DOUBLE_FIELD_TYPE,  //!< double-precision floating point field type
  INT32_FIELD_TYPE,   //!< fixed width 32-bit integer field type
  INT64_FIELD_TYPE,   //!< fixed width 64-bit integer field type

  NUMBER_OF_FIELD_TYPES  //!< NUMBER_OF_FIELD_TYPES
};

/*!
 * \brief Field traits struct to map a C++ primitive type to a FieldType
 */
template <typename FieldType>
struct field_traits
{
  static constexpr int type() { return UNDEFINED_FIELD_TYPE; };
};

/// \name Specialization of field_traits
/// @{

//------------------------------------------------------------------------------
template <>
struct field_traits<axom::float32>
{
  static constexpr int type() { return FLOAT_FIELD_TYPE; };
};

//------------------------------------------------------------------------------
template <>
struct field_traits<axom::float64>
{
  static constexpr int type() { return DOUBLE_FIELD_TYPE; };
};

//------------------------------------------------------------------------------
template <>
struct field_traits<std::int32_t>
{
  static constexpr int type() { return INT32_FIELD_TYPE; };
};

//------------------------------------------------------------------------------
template <>
struct field_traits<std::int64_t>
{
  static constexpr int type() { return INT64_FIELD_TYPE; };
};

/// @}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_FIELDTYPES_HPP_ */
