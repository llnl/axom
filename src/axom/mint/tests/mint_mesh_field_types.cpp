// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mint/mesh/FieldTypes.hpp"

#include "axom/core/Types.hpp"
#include "axom/slic.hpp"

#include "gtest/gtest.h"

// namespace aliases
namespace mint = axom::mint;

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_mesh_fieldtypes, field_traits)
{
  using int32 = std::int32_t;
  using int64 = std::int64_t;

  EXPECT_EQ(mint::field_traits<char>::type(), mint::UNDEFINED_FIELD_TYPE);
  EXPECT_EQ(mint::field_traits<float>::type(), mint::FLOAT_FIELD_TYPE);
  EXPECT_EQ(mint::field_traits<double>::type(), mint::DOUBLE_FIELD_TYPE);
  EXPECT_EQ(mint::field_traits<int>::type(), mint::INT32_FIELD_TYPE);
  EXPECT_EQ(mint::field_traits<int32>::type(), mint::INT32_FIELD_TYPE);
  EXPECT_EQ(mint::field_traits<int64>::type(), mint::INT64_FIELD_TYPE);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
