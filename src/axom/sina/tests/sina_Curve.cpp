// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/Curve.hpp"
#include "axom/sina/core/ConduitUtil.hpp"
#include "axom/sina/tests/SinaMatchers.hpp"

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

using ::testing::ContainerEq;
using ::testing::ElementsAre;

TEST(Curve, createFromVector)
{
  std::vector<double> values {1, 2, 3, 4, 5, 6};
  Curve const curve {"theName", values};
  EXPECT_EQ("theName", curve.getName());
  EXPECT_THAT(curve.getValues(), ElementsAre(1, 2, 3, 4, 5, 6));
  EXPECT_EQ("", curve.getUnits());
  EXPECT_THAT(curve.getTags(), ElementsAre());
}

TEST(Curve, createFromPointer)
{
  double values[] {1, 2, 3, 4, 5, 6};
  Curve const curve {"theName", values, sizeof(values) / sizeof(double)};
  EXPECT_EQ("theName", curve.getName());
  EXPECT_THAT(curve.getValues(), ElementsAre(1, 2, 3, 4, 5, 6));
  EXPECT_EQ("", curve.getUnits());
  EXPECT_THAT(curve.getTags(), ElementsAre());
}

TEST(Curve, createFromInitializerList)
{
  Curve const curve {"theName", {1, 2, 3, 4, 5, 6}};
  EXPECT_EQ("theName", curve.getName());
  EXPECT_THAT(curve.getValues(), ElementsAre(1, 2, 3, 4, 5, 6));
  EXPECT_EQ("", curve.getUnits());
  EXPECT_THAT(curve.getTags(), ElementsAre());
}

TEST(Curve, setUnits)
{
  Curve curve {"theName", {1, 2, 3}};
  EXPECT_EQ("", curve.getUnits());
  curve.setUnits("cm");
  EXPECT_EQ("cm", curve.getUnits());
}

TEST(Curve, setTags)
{
  Curve curve {"theName", {1, 2, 3}};
  EXPECT_THAT(curve.getTags(), ElementsAre());
  curve.setTags({"t1", "t2", "t3"});
  EXPECT_THAT(curve.getTags(), ElementsAre("t1", "t2", "t3"));
}

TEST(Curve, createFromNode_requiredOnly)
{
  conduit::Node curveAsNode = parseJsonValue(R"(
    {
        "value": [1.0, 2.0, 3.0]
    }
    )");
  Curve curve {"theName", curveAsNode};
  EXPECT_EQ("theName", curve.getName());
  EXPECT_THAT(curve.getValues(), ElementsAre(1, 2, 3));
  EXPECT_EQ("", curve.getUnits());
  EXPECT_THAT(curve.getTags(), ElementsAre());
}

TEST(Curve, createFromNode_optionalFields)
{
  conduit::Node curveAsNode = parseJsonValue(R"(
    {
        "value": [1.0, 2.0, 3.0],
        "units": "cm",
        "tags": ["t1", "t2", "t3"]
    }
    )");
  Curve curve {"theName", curveAsNode};
  EXPECT_EQ("theName", curve.getName());
  EXPECT_THAT(curve.getValues(), ElementsAre(1, 2, 3));
  EXPECT_EQ("cm", curve.getUnits());
  EXPECT_THAT(curve.getTags(), ElementsAre("t1", "t2", "t3"));
}

TEST(Curve, toNode_requiredOnly)
{
  Curve const curve {"theName", {1, 2, 3, 4}};
  std::string expected = (R"({
        "value": [1.0, 2.0, 3.0, 4.0]
    })");
  EXPECT_THAT(curve.toNode(), MatchesJsonMatcher(expected));
}

TEST(Curve, toNode_optionalFields)
{
  Curve curve {"theName", {1, 2, 3, 4}};
  curve.setUnits("cm");
  curve.setTags({"t1", "t2", "t3"});
  std::string expected = R"({
        "value": [1.0, 2.0, 3.0, 4.0],
        "units": "cm",
        "tags": ["t1", "t2", "t3"]
    })";
  EXPECT_THAT(curve.toNode(), MatchesJsonMatcher(expected));
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
