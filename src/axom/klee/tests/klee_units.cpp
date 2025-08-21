// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Units.hpp"
#include "axom/klee/KleeError.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace axom
{
namespace klee
{
using ::testing::DoubleEq;
using ::testing::ElementsAre;
using ::testing::HasSubstr;
using ::testing::Matches;

struct Length
{
  double value;
  LengthUnit units;
};

bool convertsTo(Length length1, Length length2)
{
  double length2InLength1Units = convert(length2.value, length2.units, length1.units);
  return Matches(DoubleEq(length1.value))(length2InLength1Units);
}

bool areEquivalent(Length length1, Length length2)
{
  return convertsTo(length1, length2) && convertsTo(length2, length1);
}

TEST(Units, parseLengthUnits)
{
  EXPECT_EQ(LengthUnit::km, parseLengthUnits("km", Path {}));
  EXPECT_EQ(LengthUnit::m, parseLengthUnits("m", Path {}));
  EXPECT_EQ(LengthUnit::dm, parseLengthUnits("dm", Path {}));
  EXPECT_EQ(LengthUnit::cm, parseLengthUnits("cm", Path {}));
  EXPECT_EQ(LengthUnit::mm, parseLengthUnits("mm", Path {}));
  EXPECT_EQ(LengthUnit::um, parseLengthUnits("um", Path {}));
  EXPECT_EQ(LengthUnit::nm, parseLengthUnits("nm", Path {}));
  EXPECT_EQ(LengthUnit::angstrom, parseLengthUnits("A", Path {}));
  EXPECT_EQ(LengthUnit::miles, parseLengthUnits("miles", Path {}));
  EXPECT_EQ(LengthUnit::feet, parseLengthUnits("ft", Path {}));
  EXPECT_EQ(LengthUnit::feet, parseLengthUnits("feet", Path {}));
  EXPECT_EQ(LengthUnit::inches, parseLengthUnits("in", Path {}));
  EXPECT_EQ(LengthUnit::inches, parseLengthUnits("inches", Path {}));
  EXPECT_EQ(LengthUnit::mils, parseLengthUnits("mils", Path {}));

  Path errorPath {"some/path"};
  try
  {
    parseLengthUnits("bad_units", errorPath);
    FAIL() << "Should have thrown";
  }
  catch(const KleeError &error)
  {
    ASSERT_EQ(1u, error.getErrors().size());
    EXPECT_THAT(error.getErrors()[0].message, HasSubstr("bad_units"));
    EXPECT_EQ(errorPath, error.getErrors()[0].path);
  }
}

TEST(Units, getConversionFactor)
{
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::km, LengthUnit::m));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::m, LengthUnit::dm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::dm, LengthUnit::cm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::cm, LengthUnit::mm));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::mm, LengthUnit::um));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::um, LengthUnit::nm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::nm, LengthUnit::angstrom));
  EXPECT_DOUBLE_EQ(2.54, getConversionFactor(LengthUnit::inches, LengthUnit::cm));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::inches, LengthUnit::mils));
  EXPECT_DOUBLE_EQ(12, getConversionFactor(LengthUnit::feet, LengthUnit::inches));
  EXPECT_DOUBLE_EQ(5280, getConversionFactor(LengthUnit::miles, LengthUnit::feet));

  EXPECT_DOUBLE_EQ(1, getConversionFactor(LengthUnit::miles, LengthUnit::miles));

  EXPECT_THROW(getConversionFactor(LengthUnit::cm, LengthUnit::unspecified), std::invalid_argument);
  EXPECT_THROW(getConversionFactor(LengthUnit::unspecified, LengthUnit::cm), std::invalid_argument);
}

TEST(Units, convert_adjacent)
{
  EXPECT_TRUE(areEquivalent({1, LengthUnit::km}, {1000, LengthUnit::m}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::m}, {10, LengthUnit::dm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::dm}, {10, LengthUnit::cm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::cm}, {10, LengthUnit::mm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::mm}, {1000, LengthUnit::um}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::um}, {1000, LengthUnit::nm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::nm}, {10, LengthUnit::angstrom}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::inches}, {2.54, LengthUnit::cm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::inches}, {1000, LengthUnit::mils}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::feet}, {12, LengthUnit::inches}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::miles}, {5280, LengthUnit::feet}));
}

TEST(Units, convert_all)
{
  Length equivalentLengths[] = {
    {0.254254e-3, LengthUnit::km},
    {0.254254, LengthUnit::m},
    {0.254254e1, LengthUnit::dm},
    {0.254254e2, LengthUnit::cm},
    {0.254254e3, LengthUnit::mm},
    {0.254254e6, LengthUnit::um},
    {0.254254e9, LengthUnit::nm},
    {0.254254e10, LengthUnit::angstrom},
    {10.01, LengthUnit::inches},
    {10.01e3, LengthUnit::mils},
    {10.01 / 12, LengthUnit::feet},
    {10.01 / 12 / 5280, LengthUnit::miles},
  };

  for(auto source : equivalentLengths)
  {
    for(auto target : equivalentLengths)
    {
      EXPECT_TRUE(areEquivalent(source, target));
    }
  }
}

TEST(Units, convertAll)
{
  std::vector<double> values {1.0, 2.0, 3.0};
  convertAll(values, LengthUnit::m, LengthUnit::cm);
  EXPECT_THAT(values, ElementsAre(100, 200, 300));
}

}  // namespace klee
}  // namespace axom