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
#include <type_traits>
#include <utility>

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

TEST(slam_FieldRegistry, view_buffer_add_get_find)
{
  IndexRegistry reg;

  int data[4] = {10, 20, 30, 40};
  reg.addBufferView("view", data, SetType::PositionType {4});

  EXPECT_TRUE(reg.hasBufferView("view"));
  EXPECT_EQ(reg.getBufferView("view")[2], 30);

  auto hit = reg.findBufferView("view");
  ASSERT_TRUE(hit.has_value());
  EXPECT_EQ(hit->get()[3], 40);

  const IndexRegistry& creg = reg;
  auto chit = creg.findBufferView("view");
  ASSERT_TRUE(chit.has_value());
  EXPECT_EQ(chit->get()[0], 10);
}

TEST(slam_FieldRegistry, view_field_is_non_owning)
{
  SetType s(5);
  ScalarRegistry reg;

  double data[5] = {1., 2., 3., 4., 5.};
  auto& field = reg.addFieldView("temp", &s, data);

  static_assert(
    std::is_same_v<typename std::remove_reference_t<decltype(field)>::IndirectionPolicy::IndirectionBufferType,
                   axom::ArrayView<double>>,
    "view fields use ArrayView-backed indirection");
  using ViewFieldType = std::remove_reference_t<decltype(field)>;
  static_assert(std::is_same_v<typename ViewFieldType::ConstValueType, double&>,
                "const ArrayView-backed fields preserve mutable view semantics");

  using ConstViewIndirection =
    slam::policies::ArrayViewIndirection<SetType::PositionType, const double>;
  static_assert(std::is_same_v<typename ConstViewIndirection::IndirectionResult, const double&>,
                "ArrayView<const T> exposes immutable values");
  static_assert(std::is_same_v<typename ConstViewIndirection::ConstIndirectionResult, const double&>,
                "const ArrayView<const T> fields expose immutable values");

  EXPECT_TRUE(reg.hasFieldView("temp"));
  EXPECT_DOUBLE_EQ(reg.getFieldView("temp")[1], 2.0);

  // Mutating through the map updates the backing buffer.
  reg.getFieldView("temp")[2] = 42.0;
  EXPECT_DOUBLE_EQ(data[2], 42.0);

  // Mutating the backing buffer is visible through the map.
  data[0] = -3.0;
  EXPECT_DOUBLE_EQ(reg.getFieldView("temp")[0], -3.0);

  const ScalarRegistry& creg = reg;
  const auto& constField = creg.getFieldView("temp");
  EXPECT_DOUBLE_EQ(constField[0], -3.0);

  auto hit = creg.findFieldView("temp");
  ASSERT_TRUE(hit.has_value());
  EXPECT_EQ(hit->get().size(), 5);
  EXPECT_DOUBLE_EQ(hit->get()[2], 42.0);
}
