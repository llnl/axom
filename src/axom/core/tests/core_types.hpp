// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/NumericLimits.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <type_traits>  // for std::is_same, std::is_integral, etc.

#ifndef AXOM_USE_MPI

// MPI definitions so that the tests work without MPI
using MPI_Datatype = int;
constexpr int MPI_INT8_T = -1;
constexpr int MPI_UINT8_T = -1;
constexpr int MPI_INT16_T = -1;
constexpr int MPI_UINT16_T = -1;
constexpr int MPI_INT32_T = -1;
constexpr int MPI_UINT32_T = -1;
constexpr int MPI_INT64_T = -1;
constexpr int MPI_UINT64_T = -1;
constexpr int MPI_DOUBLE = -1;
constexpr int MPI_FLOAT = -1;

#endif

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------

template <typename AxomType>
void check_mpi_type(std::size_t expected_num_bytes, MPI_Datatype expected_mpi_type)
{
#ifdef AXOM_USE_MPI

  EXPECT_FALSE(axom::mpi_traits<AxomType>::type == MPI_DATATYPE_NULL);
  EXPECT_TRUE(axom::mpi_traits<AxomType>::type == expected_mpi_type);

  int actual_bytes = 0;
  MPI_Type_size(axom::mpi_traits<AxomType>::type, &actual_bytes);
  EXPECT_EQ(static_cast<std::size_t>(actual_bytes), expected_num_bytes);

#else

  AXOM_UNUSED_VAR(expected_num_bytes);
  AXOM_UNUSED_VAR(expected_mpi_type);

#endif /* AXOM_USE_MPI */
}

//------------------------------------------------------------------------------
template <typename RealType>
void check_real_type(std::size_t expected_num_bytes, MPI_Datatype expected_mpi_type)
{
  EXPECT_TRUE(std::is_floating_point<RealType>::value);
  EXPECT_TRUE(axom::numeric_limits<RealType>::is_signed);
  EXPECT_EQ(sizeof(RealType), expected_num_bytes);

  check_mpi_type<RealType>(expected_num_bytes, expected_mpi_type);
}

//------------------------------------------------------------------------------
template <typename IntegralType>
void check_integral_type(std::size_t expected_num_bytes,
                         bool is_signed,
                         int expected_num_digits,
                         MPI_Datatype expected_mpi_type)
{
  EXPECT_TRUE(axom::numeric_limits<IntegralType>::is_integer);
  EXPECT_EQ(axom::numeric_limits<IntegralType>::is_signed, is_signed);
  EXPECT_EQ(axom::numeric_limits<IntegralType>::digits, expected_num_digits);
  EXPECT_EQ(sizeof(IntegralType), expected_num_bytes);

  check_mpi_type<IntegralType>(expected_num_bytes, expected_mpi_type);
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

TEST(core_types, check_int8)
{
  constexpr std::size_t EXP_BYTES = 1;
  constexpr int NUM_DIGITS = 7;
  constexpr bool IS_SIGNED = true;

  check_integral_type<std::int8_t>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_INT8_T);
}

//------------------------------------------------------------------------------
TEST(core_types, check_uint8)
{
  constexpr std::size_t EXP_BYTES = 1;
  constexpr int NUM_DIGITS = 8;
  constexpr bool IS_SIGNED = false;

  check_integral_type<std::uint8_t>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_UINT8_T);
}

//------------------------------------------------------------------------------
TEST(core_types, check_int16)
{
  constexpr std::size_t EXP_BYTES = 2;
  constexpr int NUM_DIGITS = 15;
  constexpr bool IS_SIGNED = true;

  check_integral_type<std::int16_t>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_INT16_T);
}

//------------------------------------------------------------------------------
TEST(core_types, check_uint16)
{
  constexpr std::size_t EXP_BYTES = 2;
  constexpr int NUM_DIGITS = 16;
  constexpr bool IS_SIGNED = false;

  check_integral_type<std::uint16_t>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_UINT16_T);
}

//------------------------------------------------------------------------------
TEST(core_types, check_int32)
{
  constexpr std::size_t EXP_BYTES = 4;
  constexpr int NUM_DIGITS = 31;
  constexpr bool IS_SIGNED = true;

  check_integral_type<std::int32_t>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_INT32_T);
}

//------------------------------------------------------------------------------
TEST(core_types, check_uint32)
{
  constexpr std::size_t EXP_BYTES = 4;
  constexpr int NUM_DIGITS = 32;
  constexpr bool IS_SIGNED = false;

  check_integral_type<std::uint32_t>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_UINT32_T);
}

#ifndef AXOM_NO_INT64_T

//------------------------------------------------------------------------------
TEST(core_types, check_int64)
{
  constexpr std::size_t EXP_BYTES = 8;
  constexpr int NUM_DIGITS = 63;
  constexpr bool IS_SIGNED = true;

  check_integral_type<std::int64_t>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_INT64_T);
}

//------------------------------------------------------------------------------
TEST(core_types, check_uint64)
{
  constexpr std::size_t EXP_BYTES = 8;
  constexpr int NUM_DIGITS = 64;
  constexpr bool IS_SIGNED = false;

  check_integral_type<std::uint64_t>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_UINT64_T);
}

//------------------------------------------------------------------------------
TEST(core_types, check_float32)
{
  constexpr std::size_t EXP_BYTES = 4;
  check_real_type<axom::float32>(EXP_BYTES, MPI_FLOAT);
}

//------------------------------------------------------------------------------
TEST(core_types, check_float64)
{
  constexpr std::size_t EXP_BYTES = 8;
  check_real_type<axom::float64>(EXP_BYTES, MPI_DOUBLE);
}

#endif /* AXOM_NO_INT64_T */

//------------------------------------------------------------------------------
TEST(core_types, check_indextype)
{
  constexpr bool IS_SIGNED = true;

#ifdef AXOM_USE_64BIT_INDEXTYPE

  constexpr bool is_int64 = std::is_same<axom::IndexType, std::int64_t>::value;
  EXPECT_TRUE(is_int64);

  constexpr std::size_t EXP_BYTES = 8;
  constexpr int NUM_DIGITS = 63;
  check_integral_type<axom::IndexType>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_INT64_T);
#else

  constexpr bool is_int32 = std::is_same<axom::IndexType, std::int32_t>::value;
  EXPECT_TRUE(is_int32);

  constexpr std::size_t EXP_BYTES = 4;
  constexpr int NUM_DIGITS = 31;
  check_integral_type<axom::IndexType>(EXP_BYTES, IS_SIGNED, NUM_DIGITS, MPI_INT32_T);

#endif
}
