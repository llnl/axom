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

#include "axom/slic.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/FieldRegistry.hpp"

#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

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

TEST(slam_FieldRegistry, heterogeneous_lookup_with_string_view_key)
{
  ScalarRegistry reg;
  reg.addScalar("energy", 1.0);

  // The lookup tables use transparent comparison (std::less<>), so a std::string_view or a const char*
  // resolves directly against the std::string keys. (A non-transparent std::map would require
  // an implicit std::string_view -> std::string conversion, which is explicit and would not compile.
  // So the fact that these calls compile and resolve is the guarantees that no temporary std::string key
  // is being constructed for lookup.)
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
  auto& field = reg.addField("temp", &s, data);

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

  EXPECT_TRUE(reg.hasField("temp"));
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

TEST(slam_FieldRegistry, field_storage_modes_share_keyspace)
{
  SetType s(3);
  ScalarRegistry reg;

  double data[3] = {10., 20., 30.};
  auto& viewField = reg.addField("density", &s, data);

  EXPECT_TRUE(reg.hasField("density"));
  EXPECT_TRUE(reg.hasFieldView("density"));
  EXPECT_FALSE(reg.findField("density").has_value());
  ASSERT_TRUE(reg.findFieldView("density").has_value());

  viewField[1] = 21.;
  EXPECT_DOUBLE_EQ(data[1], 21.);

  auto& owningField = reg.addField("density", &s);
  EXPECT_TRUE(reg.hasField("density"));
  EXPECT_FALSE(reg.hasFieldView("density"));
  EXPECT_FALSE(reg.findFieldView("density").has_value());
  ASSERT_TRUE(reg.findField("density").has_value());

  owningField[1] = 99.;
  EXPECT_DOUBLE_EQ(reg.getField("density")[1], 99.);
  EXPECT_DOUBLE_EQ(data[1], 21.);

  auto& viewFieldAgain = reg.addFieldView("density", &s, data);
  EXPECT_TRUE(reg.hasField("density"));
  EXPECT_TRUE(reg.hasFieldView("density"));
  EXPECT_FALSE(reg.findField("density").has_value());
  EXPECT_DOUBLE_EQ(viewFieldAgain[1], 21.);
}

TEST(slam_FieldRegistry, get_wrong_storage_mode_is_a_precondition_violation)
{
  // getField()/getFieldView() fetch a specific std::variant alternative. Asking for the
  // wrong alternative for an existing key is a precondition violation: it asserts in debug
  // builds (via the verify*Key helpers) and throws std::bad_variant_access in release builds.
  SetType s(3);
  ScalarRegistry reg;

  double data[3] = {10., 20., 30.};
  reg.addFieldView("density", &s, data);  // view-backed
  reg.addField("temperature", &s);        // owning

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so these checks are skipped there.
  EXPECT_DEATH_IF_SUPPORTED(reg.getField("density"), "");
  EXPECT_DEATH_IF_SUPPORTED(reg.getFieldView("temperature"), "");

  const ScalarRegistry& creg = reg;
  EXPECT_DEATH_IF_SUPPORTED(creg.getField("density"), "");
  EXPECT_DEATH_IF_SUPPORTED(creg.getFieldView("temperature"), "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
  // In release builds the underlying std::get throws on the wrong alternative.
  EXPECT_THROW(reg.getField("density"), std::bad_variant_access);
  EXPECT_THROW(reg.getFieldView("temperature"), std::bad_variant_access);
#endif

  // The safe, non-throwing discriminators agree on the storage mode in all build types.
  EXPECT_FALSE(reg.findField("density").has_value());
  EXPECT_TRUE(reg.findFieldView("density").has_value());
  EXPECT_TRUE(reg.findField("temperature").has_value());
  EXPECT_FALSE(reg.findFieldView("temperature").has_value());
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
