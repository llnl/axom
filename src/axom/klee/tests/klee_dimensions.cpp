// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Dimensions.hpp"

#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "gtest/gtest.h"
#include <iostream>

namespace klee = axom::klee;

TEST(Dimensions, two)
{
  klee::Dimensions two {klee::Dimensions::Two};
  EXPECT_EQ("Two", axom::fmt::format("{}", two));
  std::cout << "The dimension is " << axom::fmt::format("{}", two) << "\n";
}

TEST(Dimensions, three)
{
  klee::Dimensions three {klee::Dimensions::Three};
  EXPECT_EQ("Three", axom::fmt::format("{}", three));
  std::cout << "The dimension is " << axom::fmt::format("{}", three) << "\n";
}

TEST(Dimensions, unspecified)
{
  klee::Dimensions def {};
  EXPECT_EQ("Unspecified", axom::fmt::format("{}", def));
  std::cout << "The dimension is " << axom::fmt::format("{}", def) << "\n";

  klee::Dimensions unspec {klee::Dimensions::Unspecified};
  EXPECT_EQ("Unspecified", axom::fmt::format("{}", unspec));
  std::cout << "The dimension is " << axom::fmt::format("{}", unspec) << "\n";
}

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  int result = RUN_ALL_TESTS();
  return result;
}
