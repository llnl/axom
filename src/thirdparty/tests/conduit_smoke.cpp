// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: conduit_smoke.cpp
///
//-----------------------------------------------------------------------------

#include "conduit.hpp"

#include <iostream>
#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(conduit_smoke, basic_use)
{
  EXPECT_EQ(sizeof(conduit::uint32), (size_t)4);
  EXPECT_EQ(sizeof(conduit::uint64), (size_t)8);
  EXPECT_EQ(sizeof(conduit::float64), (size_t)8);

  std::cout << conduit::about() << std::endl;
}
