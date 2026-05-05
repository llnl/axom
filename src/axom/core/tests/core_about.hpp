// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/utilities/About.hpp"

#include <sstream>
#include <string>

TEST(core_about, print_about) { axom::about(); }

//-----------------------------------------------------------------------------
TEST(core_about, get_version)
{
  std::ostringstream expected_version;
  expected_version << AXOM_VERSION_FULL;
  std::string sha = axom::gitSHA();
  if(!sha.empty())
  {
    expected_version << "-" << sha;
  }

  std::string axom_version = axom::getVersion();
  EXPECT_EQ(expected_version.str(), axom_version);

  std::cout << "Version: " << axom::getVersion() << std::endl;
}
