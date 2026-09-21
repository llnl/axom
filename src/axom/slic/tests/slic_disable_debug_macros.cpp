// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Exercise the explicit Slic override when AXOM_DEBUG is enabled.
#ifndef AXOM_DEBUG
  #define AXOM_DEBUG 1
#endif

#ifdef AXOM_ENABLE_SLIC_DEBUG_MACROS
  #undef AXOM_ENABLE_SLIC_DEBUG_MACROS
#endif
#define AXOM_ENABLE_SLIC_DEBUG_MACROS 0

#include "axom/slic.hpp"

#include "gtest/gtest.h"

TEST(slic_disable_debug_macros, disabled_with_axom_debug)
{
  axom::slic::SimpleLogger logger;
  int evaluation_count = 0;

  SLIC_ASSERT(++evaluation_count == 1);
  SLIC_CHECK(++evaluation_count == 1);
  SLIC_DEBUG(++evaluation_count);

  EXPECT_EQ(evaluation_count, 0);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
