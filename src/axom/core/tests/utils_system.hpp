// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/System.hpp"

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(utils_system, getUserName)
{
  // Returns empty string on failure
  std::string user_name = axom::utilities::getUserName();
  std::cout << "user name = " << user_name << std::endl;
  EXPECT_TRUE(user_name != "");
}

TEST(utils_system, getHostName)
{
  // Returns empty string on failure
  std::string host_name = axom::utilities::getHostName();
  std::cout << "host name = " << host_name << std::endl;
  EXPECT_TRUE(host_name != "");
}
