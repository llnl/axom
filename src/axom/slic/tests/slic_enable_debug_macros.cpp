// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Exercise the Slic debug macro guard independently of AXOM_DEBUG.
#ifdef AXOM_DEBUG
  #undef AXOM_DEBUG
#endif

#define AXOM_ENABLE_SLIC_DEBUG_MACROS 1

#include "axom/slic.hpp"

#include "gtest/gtest.h"

TEST(slic_enable_debug_macros, enabled_without_axom_debug)
{
  axom::slic::SimpleLogger logger;
  int evaluation_count = 0;

  SLIC_ASSERT(++evaluation_count == 1);
  SLIC_CHECK(++evaluation_count == 2);
  SLIC_DEBUG(++evaluation_count);

  EXPECT_EQ(evaluation_count, 3);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
