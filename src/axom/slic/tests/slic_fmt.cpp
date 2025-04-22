// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file slic_fmt.cpp
 *
 * \brief A simple test to see if we can use the fmt library within slic macros.
 */

#include "axom/config.hpp"
#include "axom/fmt.hpp"
#include "axom/slic.hpp"

#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(slic_fmt, basic_use)
{
  SLIC_INFO("Formatting with C++ streams and with "
            << axom::fmt::format("the '{}' string formatting library", "fmt")
            << " can be used together in slic macros.");

  // Style similar to Python format
  SLIC_INFO(axom::fmt::format("{1} {0}", "world", "Hello"));
  EXPECT_EQ("Hello world", axom::fmt::format("{1} {0}", "world", "Hello"));

#if !defined(__INTEL_COMPILER) && __INTEL_COMPILER != 1800
  // Note: This test seg faults on intel 18

  // Python style w/ 'dictionary'
  SLIC_INFO(axom::fmt::format("Hello, {name}! Goodbye, {name}.", axom::fmt::arg("name", "Axom")));
  EXPECT_EQ("Hello, Axom! Goodbye, Axom.",
            axom::fmt::format("Hello, {name}! Goodbye, {name}.", axom::fmt::arg("name", "Axom")));
#endif

  // Python style with additional formatting
  SLIC_INFO(axom::fmt::format("int:{0:d};  hex:{0:x};  oct:{0:o}; bin:{0:b}", 42));
  EXPECT_EQ("int:42;  hex:2a;  oct:52; bin:101010",
            axom::fmt::format("int:{0:d};  hex:{0:x};  oct:{0:o}; bin:{0:b}", 42));

  // sprintf style
  SLIC_INFO(axom::fmt::sprintf("Two significant digits: %.2f", 1.234567));
  EXPECT_EQ("1.23", axom::fmt::sprintf("%.2f", 1.234567));
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
