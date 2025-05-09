// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "caliper/caliper-config.h"
#include "caliper/cali.h"
#include "caliper/Caliper.h"
#include "caliper/ConfigManager.h"

#include <iostream>

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(caliper_smoke, check_version)
{
  std::cout << "Using caliper version: "     //
            << CALIPER_MAJOR_VERSION << "."  //
            << CALIPER_MINOR_VERSION << "."  //
            << CALIPER_PATCH_VERSION << std::endl;

  EXPECT_TRUE(CALIPER_MAJOR_VERSION >= 0);
  EXPECT_TRUE(CALIPER_MINOR_VERSION >= 0);
  EXPECT_TRUE(CALIPER_PATCH_VERSION >= 0);
}

TEST(caliper_smoke, config_manager)
{
  cali::ConfigManager mgr;
  mgr.add("runtime-report(output=stdout,calc.inclusive=true)");

  // Check for configuration string parse errors
  if(mgr.error())
  {
    FAIL() << "ConfigManager: " << mgr.error_msg() << std::endl;
  }

  mgr.start();

  // let's do something with annotations
  {
    CALI_MARK_BEGIN("my region");

    for(int i = 0; i < 10; ++i)
    {
      CALI_CXX_MARK_SCOPE("inner1");
    }

    for(int i = 0; i < 100; ++i)
    {
      CALI_CXX_MARK_SCOPE("inner2");
    }

    CALI_MARK_END("my region");
  }

  mgr.flush();
}
