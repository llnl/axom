// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_TESTS_CORE_CONSTEXPR_ASSERT_HPP_
#define AXOM_CORE_TESTS_CORE_CONSTEXPR_ASSERT_HPP_

#include "gtest/gtest.h"

#include "axom/core/Macros.hpp"

namespace
{
constexpr int checked_add(int a, int b)
{
  AXOM_CONSTEXPR_ASSERT(a >= 0);
  AXOM_CONSTEXPR_ASSERT(b >= 0);
  return a + b;
}

static_assert(checked_add(1, 2) == 3, "AXOM_CONSTEXPR_ASSERT works in constant evaluation");
}  // namespace

TEST(core_constexpr_assert, usable_in_constexpr)
{
  constexpr int v = checked_add(4, 5);
  EXPECT_EQ(v, 9);
}

TEST(core_constexpr_assert, runtime_true_noop)
{
  AXOM_CONSTEXPR_ASSERT(true);
  SUCCEED();
}

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
TEST(core_constexpr_assert, runtime_false_death)
{
  EXPECT_DEATH_IF_SUPPORTED([]() { AXOM_CONSTEXPR_ASSERT(false); }(), ".*");
}
#endif

#endif  // AXOM_CORE_TESTS_CORE_CONSTEXPR_ASSERT_HPP_
