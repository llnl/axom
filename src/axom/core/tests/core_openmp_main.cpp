// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"  // for compile-time definitions

#include "core_openmp_map.hpp"

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  return RUN_ALL_TESTS();
}
