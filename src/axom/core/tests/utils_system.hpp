// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
