// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/Units.hpp"

#include "gtest/gtest.h"

#include <stdexcept>
#include <string>
#include <vector>

namespace axom
{
namespace utilities
{
struct Length
{
  double value;
  LengthUnit units;
};

void expectConvertsTo(Length length1, Length length2)
{
  double length2InLength1Units = convert(length2.value, length2.units, length1.units);
  EXPECT_DOUBLE_EQ(length1.value, length2InLength1Units);
}

void expectEquivalent(Length length1, Length length2)
{
  expectConvertsTo(length1, length2);
  expectConvertsTo(length2, length1);
}

TEST(Units, getLengthUnit)
{
  EXPECT_EQ(LengthUnit::km, getLengthUnit("km"));
  EXPECT_EQ(LengthUnit::hm, getLengthUnit("hectometer"));
  EXPECT_EQ(LengthUnit::dam, getLengthUnit("decametre"));
  EXPECT_EQ(LengthUnit::m, getLengthUnit("m"));
  EXPECT_EQ(LengthUnit::m, getLengthUnit("METRES"));
  EXPECT_EQ(LengthUnit::dm, getLengthUnit("dm"));
  EXPECT_EQ(LengthUnit::cm, getLengthUnit("cm"));
  EXPECT_EQ(LengthUnit::cm, getLengthUnit("centimeter"));
  EXPECT_EQ(LengthUnit::mm, getLengthUnit("mm"));
  EXPECT_EQ(LengthUnit::um, getLengthUnit("um"));
  EXPECT_EQ(LengthUnit::um, getLengthUnit("micron"));
  EXPECT_EQ(LengthUnit::nm, getLengthUnit("nm"));
  EXPECT_EQ(LengthUnit::pm, getLengthUnit("picometers"));
  EXPECT_EQ(LengthUnit::fm, getLengthUnit("femtometre"));
  EXPECT_EQ(LengthUnit::fm, getLengthUnit("femtometres"));
  EXPECT_EQ(LengthUnit::am, getLengthUnit("attometer"));
  EXPECT_EQ(LengthUnit::angstrom, getLengthUnit("A"));
  EXPECT_EQ(LengthUnit::miles, getLengthUnit("miles"));
  EXPECT_EQ(LengthUnit::miles, getLengthUnit("mi"));
  EXPECT_EQ(LengthUnit::feet, getLengthUnit("ft"));
  EXPECT_EQ(LengthUnit::feet, getLengthUnit("feet"));
  EXPECT_EQ(LengthUnit::inches, getLengthUnit("in"));
  EXPECT_EQ(LengthUnit::inches, getLengthUnit("inches"));
  EXPECT_EQ(LengthUnit::mils, getLengthUnit("mils"));

  try
  {
    getLengthUnit("bad_units");
    FAIL() << "Should have thrown";
  }
  catch(const std::invalid_argument& error)
  {
    EXPECT_NE(std::string::npos, std::string(error.what()).find("bad_units"));
  }
}

TEST(Units, getConversionFactor)
{
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::km, LengthUnit::m));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::km, LengthUnit::hm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::hm, LengthUnit::dam));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::dam, LengthUnit::m));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::m, LengthUnit::dm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::dm, LengthUnit::cm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::cm, LengthUnit::mm));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::mm, LengthUnit::um));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::um, LengthUnit::nm));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::nm, LengthUnit::pm));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::pm, LengthUnit::fm));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::fm, LengthUnit::am));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::nm, LengthUnit::angstrom));
  EXPECT_DOUBLE_EQ(2.54, getConversionFactor(LengthUnit::inches, LengthUnit::cm));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::inches, LengthUnit::mils));
  EXPECT_DOUBLE_EQ(12, getConversionFactor(LengthUnit::feet, LengthUnit::inches));
  EXPECT_DOUBLE_EQ(5280, getConversionFactor(LengthUnit::miles, LengthUnit::feet));

  EXPECT_DOUBLE_EQ(1, getConversionFactor(LengthUnit::miles, LengthUnit::miles));

  EXPECT_THROW(getConversionFactor(LengthUnit::cm, LengthUnit::unspecified), std::invalid_argument);
  EXPECT_THROW(getConversionFactor(LengthUnit::unspecified, LengthUnit::cm), std::invalid_argument);
}

TEST(Units, getLengthUnitName)
{
  EXPECT_EQ("km", getLengthUnitName(LengthUnit::km));
  EXPECT_EQ("hm", getLengthUnitName(LengthUnit::hm));
  EXPECT_EQ("dam", getLengthUnitName(LengthUnit::dam));
  EXPECT_EQ("m", getLengthUnitName(LengthUnit::m));
  EXPECT_EQ("dm", getLengthUnitName(LengthUnit::dm));
  EXPECT_EQ("cm", getLengthUnitName(LengthUnit::cm));
  EXPECT_EQ("mm", getLengthUnitName(LengthUnit::mm));
  EXPECT_EQ("um", getLengthUnitName(LengthUnit::um));
  EXPECT_EQ("nm", getLengthUnitName(LengthUnit::nm));
  EXPECT_EQ("pm", getLengthUnitName(LengthUnit::pm));
  EXPECT_EQ("fm", getLengthUnitName(LengthUnit::fm));
  EXPECT_EQ("am", getLengthUnitName(LengthUnit::am));
  EXPECT_EQ("A", getLengthUnitName(LengthUnit::angstrom));
  EXPECT_EQ("ft", getLengthUnitName(LengthUnit::feet));
  EXPECT_EQ("in", getLengthUnitName(LengthUnit::inches));
  EXPECT_EQ("mils", getLengthUnitName(LengthUnit::mils));
  EXPECT_EQ("miles", getLengthUnitName(LengthUnit::miles));
  EXPECT_EQ("unspecified", getLengthUnitName(LengthUnit::unspecified));
}

TEST(Units, convert_adjacent)
{
  expectEquivalent({1, LengthUnit::km}, {1000, LengthUnit::m});
  expectEquivalent({1, LengthUnit::hm}, {100, LengthUnit::m});
  expectEquivalent({1, LengthUnit::dam}, {10, LengthUnit::m});
  expectEquivalent({1, LengthUnit::m}, {10, LengthUnit::dm});
  expectEquivalent({1, LengthUnit::dm}, {10, LengthUnit::cm});
  expectEquivalent({1, LengthUnit::cm}, {10, LengthUnit::mm});
  expectEquivalent({1, LengthUnit::mm}, {1000, LengthUnit::um});
  expectEquivalent({1, LengthUnit::um}, {1000, LengthUnit::nm});
  expectEquivalent({1, LengthUnit::nm}, {1000, LengthUnit::pm});
  expectEquivalent({1, LengthUnit::pm}, {1000, LengthUnit::fm});
  expectEquivalent({1, LengthUnit::fm}, {1000, LengthUnit::am});
  expectEquivalent({1, LengthUnit::nm}, {10, LengthUnit::angstrom});
  expectEquivalent({1, LengthUnit::inches}, {2.54, LengthUnit::cm});
  expectEquivalent({1, LengthUnit::inches}, {1000, LengthUnit::mils});
  expectEquivalent({1, LengthUnit::feet}, {12, LengthUnit::inches});
  expectEquivalent({1, LengthUnit::miles}, {5280, LengthUnit::feet});
}

TEST(Units, convert_all)
{
  Length equivalentLengths[] = {
    {0.254254e-3, LengthUnit::km},
    {0.254254e-2, LengthUnit::hm},
    {0.254254e-1, LengthUnit::dam},
    {0.254254, LengthUnit::m},
    {0.254254e1, LengthUnit::dm},
    {0.254254e2, LengthUnit::cm},
    {0.254254e3, LengthUnit::mm},
    {0.254254e6, LengthUnit::um},
    {0.254254e9, LengthUnit::nm},
    {0.254254e12, LengthUnit::pm},
    {0.254254e15, LengthUnit::fm},
    {0.254254e18, LengthUnit::am},
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
      expectEquivalent(source, target);
    }
  }
}

TEST(Units, convert_throws_on_unspecified)
{
  EXPECT_THROW(convert(1.0, LengthUnit::cm, LengthUnit::unspecified), std::invalid_argument);
  EXPECT_THROW(convert(1.0, LengthUnit::unspecified, LengthUnit::cm), std::invalid_argument);
}

TEST(Units, convertAll)
{
  std::vector<double> values {1.0, 2.0, 3.0};
  convertAll(values, LengthUnit::m, LengthUnit::cm);
  ASSERT_EQ(3u, values.size());
  EXPECT_DOUBLE_EQ(100, values[0]);
  EXPECT_DOUBLE_EQ(200, values[1]);
  EXPECT_DOUBLE_EQ(300, values[2]);
}

TEST(Units, convertAll_same_units_is_noop)
{
  std::vector<double> values {1.0, 2.0, 3.0};
  convertAll(values, LengthUnit::feet, LengthUnit::feet);
  ASSERT_EQ(3u, values.size());
  EXPECT_DOUBLE_EQ(1.0, values[0]);
  EXPECT_DOUBLE_EQ(2.0, values[1]);
  EXPECT_DOUBLE_EQ(3.0, values[2]);
}

TEST(Units, convertAll_imperial_to_metric)
{
  std::vector<double> values {1.0, 2.5, 10.0};
  convertAll(values, LengthUnit::feet, LengthUnit::inches);
  ASSERT_EQ(3u, values.size());
  EXPECT_DOUBLE_EQ(12.0, values[0]);
  EXPECT_DOUBLE_EQ(30.0, values[1]);
  EXPECT_DOUBLE_EQ(120.0, values[2]);
}

TEST(Units, convertAll_empty)
{
  std::vector<double> values;
  convertAll(values, LengthUnit::m, LengthUnit::cm);
  EXPECT_TRUE(values.empty());
}

TEST(Units, convertAll_throws_on_unspecified)
{
  std::vector<double> values {1.0, 2.0, 3.0};
  EXPECT_THROW(convertAll(values, LengthUnit::cm, LengthUnit::unspecified), std::invalid_argument);
  EXPECT_THROW(convertAll(values, LengthUnit::unspecified, LengthUnit::cm), std::invalid_argument);
}

}  // namespace utilities
}  // namespace axom
