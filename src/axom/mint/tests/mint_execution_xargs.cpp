// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mint/execution/xargs.hpp"  // execution xargs types/traits
#include "axom/slic.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <cstring>

namespace mint = axom::mint;
namespace xargs = mint::xargs;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
template <typename ArgType>
void check_valid()
{
  SLIC_INFO("checking [" << mint::xargs_traits<ArgType>::name() << "]");
  EXPECT_TRUE(mint::xargs_traits<ArgType>::valid());
  EXPECT_TRUE(strlen(mint::xargs_traits<ArgType>::name()) > 0);
}

//------------------------------------------------------------------------------
template <typename ArgType>
void check_invalid()
{
  SLIC_INFO("checking invalid policy!");
  EXPECT_FALSE(mint::xargs_traits<ArgType>::valid());
  EXPECT_TRUE(strlen(mint::xargs_traits<ArgType>::name()) == 0);
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_execution_xargs, check_types)
{
  check_valid<xargs::index>();
  check_valid<xargs::ij>();
  check_valid<xargs::ijk>();
  check_valid<xargs::x>();
  check_valid<xargs::xy>();
  check_valid<xargs::xyz>();
  check_valid<xargs::nodeids>();
  check_valid<xargs::coords>();
  check_valid<xargs::faceids>();
  check_valid<xargs::cellids>();

  struct invalid_type
  { };
  check_invalid<invalid_type>();
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
