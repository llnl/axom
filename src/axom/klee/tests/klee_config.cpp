// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "gtest/gtest.h"

namespace axom
{
namespace klee
{
TEST(KleeConfig, defines)
{
#ifndef AXOM_USE_KLEE
  FAIL() << "AXOM_ENABLE_KLEE not defined";
#endif
}

}  // namespace klee
}  // namespace axom
