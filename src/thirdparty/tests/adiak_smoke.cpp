// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "adiak.hpp"

#include <iostream>

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(adiak_smoke, check_version)
{
  std::cout << "Using adiak version: "     //
            << ADIAK_VERSION << "."        //
            << ADIAK_MINOR_VERSION << "."  //
            << ADIAK_POINT_VERSION << std::endl;

  EXPECT_TRUE(ADIAK_VERSION >= 0);
  EXPECT_TRUE(ADIAK_MINOR_VERSION >= 0);
  EXPECT_TRUE(ADIAK_POINT_VERSION >= 0);
}

TEST(adiak_smoke, basic)
{
  adiak::init(nullptr);

  // use a built-in value
  {
    bool ret = adiak::user();
    EXPECT_TRUE(ret);
  }

  // set a custom value
  {
    bool ret = adiak::value("test_type", "smoke");
    EXPECT_TRUE(ret);
  }

  adiak::fini();
}
