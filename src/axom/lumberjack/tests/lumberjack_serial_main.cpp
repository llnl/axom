// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "lumberjack_Lumberjack.hpp"
#include "lumberjack_Message.hpp"
#include "lumberjack_TextEqualityCombiner.hpp"
#include "lumberjack_TextTagCombiner.hpp"
#include "lumberjack_LineFileTagCombiner.hpp"

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
