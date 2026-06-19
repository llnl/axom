// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_TESTS_CORE_OPTIONAL_HPP_
#define AXOM_CORE_TESTS_CORE_OPTIONAL_HPP_

#include "gtest/gtest.h"

#include "axom/core/Optional.hpp"

namespace
{
static_assert(!axom::Optional<int>().has_value(), "default Optional is disengaged");
static_assert(axom::Optional<int>(42).has_value(), "value-constructed Optional is engaged");
static_assert(*axom::Optional<int>(42) == 42, "engaged Optional yields its value");
static_assert(axom::Optional<int>(42).value_or(-1) == 42, "value_or returns the value when engaged");
static_assert(axom::Optional<int>().value_or(-1) == -1, "value_or returns fallback when empty");
static_assert(static_cast<bool>(axom::Optional<int>(0)), "engaged-with-zero is still engaged");
static_assert(!static_cast<bool>(axom::Optional<int>()), "disengaged converts to false");

static_assert(std::is_trivially_copyable_v<axom::Optional<int>>,
              "axom::Optional<int> is trivially copyable (device-capturable)");
}  // namespace

TEST(core_optional, engaged_and_disengaged)
{
  axom::Optional<double> empty;
  EXPECT_FALSE(empty.has_value());
  EXPECT_FALSE(static_cast<bool>(empty));
  EXPECT_DOUBLE_EQ(empty.value_or(2.5), 2.5);

  axom::Optional<double> full(3.25);
  EXPECT_TRUE(full.has_value());
  EXPECT_TRUE(static_cast<bool>(full));
  EXPECT_DOUBLE_EQ(*full, 3.25);
  EXPECT_DOUBLE_EQ(full.value_or(2.5), 3.25);

  *full = 9.0;
  EXPECT_DOUBLE_EQ(*full, 9.0);
}

#endif  // AXOM_CORE_TESTS_CORE_OPTIONAL_HPP_
