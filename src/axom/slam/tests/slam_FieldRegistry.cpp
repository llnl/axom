// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_FieldRegistry.cpp
 *
 * \brief Unit tests for slam::FieldRegistry, covering the std::optional find APIs
 *  and transparent (string_view) heterogeneous lookup
 */

#include "gtest/gtest.h"

#include "axom/slam/RangeSet.hpp"
#include "axom/slam/FieldRegistry.hpp"

#include <optional>
#include <string>
#include <string_view>

namespace
{
namespace slam = axom::slam;

using SetType = slam::PositionSet<>;
using ScalarRegistry = slam::FieldRegistry<SetType, double>;
using IndexRegistry = slam::FieldRegistry<SetType, int>;

}  // anonymous namespace

TEST(slam_FieldRegistry, scalar_add_get_has)
{
  ScalarRegistry reg;
  EXPECT_FALSE(reg.hasScalar("gravity"));

  reg.addScalar("gravity", 9.81);
  EXPECT_TRUE(reg.hasScalar("gravity"));
  EXPECT_DOUBLE_EQ(reg.getScalar("gravity"), 9.81);

  // getScalar returns a mutable reference
  reg.getScalar("gravity") = 9.80665;
  EXPECT_DOUBLE_EQ(reg.getScalar("gravity"), 9.80665);
}

TEST(slam_FieldRegistry, find_scalar_optional)
{
  ScalarRegistry reg;
  reg.addScalar("dt", 0.5);

  std::optional<double> hit = reg.findScalar("dt");
  ASSERT_TRUE(hit.has_value());
  EXPECT_DOUBLE_EQ(*hit, 0.5);

  std::optional<double> miss = reg.findScalar("missing");
  EXPECT_FALSE(miss.has_value());

  // find does not insert: a missing query leaves the registry unchanged.
  EXPECT_FALSE(reg.hasScalar("missing"));
}

TEST(slam_FieldRegistry, find_field_optional)
{
  SetType s(10);
  ScalarRegistry reg;
  reg.addField("temperature", &s);

  auto hit = reg.findField("temperature");
  ASSERT_TRUE(hit.has_value());
  EXPECT_EQ(hit->get().size(), 10);

  auto miss = reg.findField("pressure");
  EXPECT_FALSE(miss.has_value());
  EXPECT_FALSE(reg.hasField("pressure"));
}

TEST(slam_FieldRegistry, find_buffer_optional)
{
  IndexRegistry reg;
  reg.addBuffer("indices", 7);

  auto hit = reg.findBuffer("indices");
  ASSERT_TRUE(hit.has_value());

  auto miss = reg.findBuffer("absent");
  EXPECT_FALSE(miss.has_value());
}

TEST(slam_FieldRegistry, heterogeneous_lookup_no_allocation)
{
  ScalarRegistry reg;
  reg.addScalar("energy", 1.0);

  // Query with a string_view and a const char* -- transparent comparison means
  // neither constructs a temporary std::string key.
  std::string_view sv = "energy";
  EXPECT_TRUE(reg.hasScalar(sv));
  EXPECT_TRUE(reg.hasScalar("energy"));
  EXPECT_DOUBLE_EQ(reg.getScalar(sv), 1.0);

  std::optional<double> hit = reg.findScalar(sv);
  ASSERT_TRUE(hit.has_value());
  EXPECT_DOUBLE_EQ(*hit, 1.0);
}

TEST(slam_FieldRegistry, nameless_keys_are_unique)
{
  SetType s(3);
  ScalarRegistry reg;

  auto& f0 = reg.addNamelessField(&s);
  auto& f1 = reg.addNamelessField(&s);
  // Distinct fields were created (distinct storage).
  EXPECT_NE(&f0, &f1);

  IndexRegistry ireg;
  auto& b0 = ireg.addNamelessBuffer(2);
  auto& b1 = ireg.addNamelessBuffer(2);
  EXPECT_NE(&b0, &b1);
}
