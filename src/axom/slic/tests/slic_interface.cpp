// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic/interface/slic.hpp"

TEST(slic_interface, initialize_finalize)
{
  EXPECT_FALSE(axom::slic::isInitialized());

  // Initialize and finalize slic
  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());

  // Check that we can Initialize and finalize slic a second time
  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());
}

TEST(slic_interface, logging_level)
{
  axom::slic::initialize();

  // test get and set for each level
  for(int i = 0; i < axom::slic::message::Num_Levels; ++i)
  {
    axom::slic::message::Level lev = static_cast<axom::slic::message::Level>(i);

    axom::slic::setLoggingMsgLevel(lev);
    EXPECT_EQ(lev, axom::slic::getLoggingMsgLevel());
  }

  axom::slic::finalize();
}
