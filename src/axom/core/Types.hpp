// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *  \file Types.hpp
 *
 *  \brief Exposes some common types used by axom components.
 */

#ifndef AXOM_TYPES_HPP_
#define AXOM_TYPES_HPP_

// Axom includes
#include "axom/config.hpp"

// C/C++ includes
#include <cstdint>

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

namespace axom
{
using float32 = float;
using float64 = double;

/*
  Axom integer types are deprecated.
  Their use can trigger an warning or an error, or be quietly allowed,
  by using CMake -DAXOM_DEPRECATED_TYPES=<WARN|ERROR|ALLOW>
  Eventually, these types will be removed.
*/
#if AXOM_DEPRECATED_TYPES_N == 1 || AXOM_DEPRECATED_TYPES_N == 2
  #if AXOM_DEPRECATED_TYPES_N == 1
    #if defined(_MSC_VER)
      #pragma message( \
        "warning: Using deprecated Axom types.  Please see CMake variable AXOM_DEPRECATED_TYPES")
    #else
      #warning "Using deprecated Axom types.  Please see CMake variable AXOM_DEPRECATED_TYPES"
    #endif
  #endif
using int8 = std::int8_t;   /*!< 8-bit signed integer type      */
using uint8 = std::uint8_t; /*!< 8-bit unsigned integer type    */

using int16 = std::int16_t;   /*!< 16-bit signed integer type     */
using uint16 = std::uint16_t; /*!< 16-bit unsigned integer type   */

using int32 = std::int32_t;   /*!< 32-bit signed integer type     */
using uint32 = std::uint32_t; /*!< 32-bit unsigned integer type   */

  // Note: KW -- We assume that AXOM_NO_INT64_T will be defined
  // on systems/compilers that do not support 64 bit integer types
  #ifndef AXOM_NO_INT64_T
using int64 = std::int64_t;   /*!< 64-bit signed integer type     */
using uint64 = std::uint64_t; /*!< 64-bit unsigned integer type   */
  #endif
#endif

#if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
using IndexType = std::int64_t;
#else
using IndexType = std::int32_t;
#endif

static constexpr IndexType InvalidIndex = -1;

/*!
 * \brief Maps a type to itself in a non-deduced context.
 *
 * This is an Axom-local equivalent of C++20's std::type_identity for C++17.
 * Its primary use is to suppress template argument deduction for a function parameter.
 * A parameter declared as type_identity_t<T> does not participate in deducing T, 
 * so T is taken from the explicit template argument or from another parameter,
 * and the argument here is implicitly converted.
 */
template <typename T>
struct type_identity
{
  using type = T;
};

template <typename T>
using type_identity_t = typename type_identity<T>::type;

#ifdef AXOM_USE_MPI

// Note: MSVC complains about uninitialized static const integer class members,
// but our nvcc/mpi combination complains about static constexpr MPI_Datatypes.
// Since it's one line per traits class, implement both ways w/ an ifdef

/// Traits class to map Axom types to their corresponding MPI type.
template <class AxomType>
struct mpi_traits
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_DATATYPE_NULL;
  #else
  static const MPI_Datatype type;
  #endif
};

/// \name Specialization of mpi_traits
/// @{
template <>
struct mpi_traits<float64>
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_DOUBLE;
  #else
  static const MPI_Datatype type;
  #endif
};

template <>
struct mpi_traits<float32>
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_FLOAT;
  #else
  static const MPI_Datatype type;
  #endif
};

template <>
struct mpi_traits<std::int8_t>
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_INT8_T;
  #else
  static const MPI_Datatype type;
  #endif
};

template <>
struct mpi_traits<std::uint8_t>
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_UINT8_T;
  #else
  static const MPI_Datatype type;
  #endif
};

template <>
struct mpi_traits<std::int16_t>
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_INT16_T;
  #else
  static const MPI_Datatype type;
  #endif
};

template <>
struct mpi_traits<std::uint16_t>
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_UINT16_T;
  #else
  static const MPI_Datatype type;
  #endif
};

template <>
struct mpi_traits<std::int32_t>
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_INT32_T;
  #else
  static const MPI_Datatype type;
  #endif
};

template <>
struct mpi_traits<std::uint32_t>
{
  #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_UINT32_T;
  #else
  static const MPI_Datatype type;
  #endif
};

  #ifndef AXOM_NO_INT64_T
template <>
struct mpi_traits<std::int64_t>
{
    #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_INT64_T;
    #else
  static const MPI_Datatype type;
    #endif
};

template <>
struct mpi_traits<std::uint64_t>
{
    #ifdef _MSC_VER
  static constexpr MPI_Datatype type = MPI_UINT64_T;
    #else
  static const MPI_Datatype type;
    #endif
};
  #endif  // AXOM_NO_INT64_T

  /// @}

#endif  // AXOM_USE_MPI

}  // end namespace axom

#endif  // AXOM_TYPES_HPP_
