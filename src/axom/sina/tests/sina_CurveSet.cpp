// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/CurveSet.hpp"
#include "axom/sina/tests/SinaMatchers.hpp"

#include <utility>
#include <unordered_map>

namespace axom
{
namespace sina
{

// NOTE: We need an operator== for the tests. For it to be able to be found
// by different matchers, it needs to be in the same namespace as Curve.
// If we need up needing it in another test, we'll have to move it to another
// file.
//
// NOTE: Since this isn't in an unnamed namespace, we need a forward
// declaration to satisfy strict compiler warnings.
bool operator==(Curve const &lhs, Curve const &rhs);

/**
 * Compare two curves for equality. All fields must be equal, including the
 * doubles in the lists of values. This is not suitable for checking any
 * calculated values.
 *
 * @param lhs the left-hand-side operand
 * @param rhs the right-hand-side operand
 * @return whether the curves are equal
 */
bool operator==(Curve const &lhs, Curve const &rhs)
{
  bool r = lhs.getName() == rhs.getName() && lhs.getUnits() == rhs.getUnits() &&
    lhs.getTags() == rhs.getTags() && lhs.getValues() == rhs.getValues();
  return r;
}

namespace testing
{
namespace
{

using ::testing::ContainerEq;
using ::testing::ElementsAre;

TEST(CurveSet, initialState)
{
  CurveSet const cs {"theName"};
  ASSERT_EQ("theName", cs.getName());
  ASSERT_TRUE(cs.getIndependentCurves().empty());
  ASSERT_TRUE(cs.getDependentCurves().empty());
  ASSERT_NE(&cs.getIndependentCurves(), &cs.getDependentCurves());
}

TEST(CurveSet, addIndependentCurves)
{
  CurveSet cs {"testSet"};
  std::unordered_map<std::string, Curve> expectedCurves;

  Curve i1 {"i1", {1, 2, 3}};
  cs.addIndependentCurve(i1);
  expectedCurves.insert(std::make_pair(i1.getName(), i1));
  EXPECT_THAT(cs.getIndependentCurves(), ContainerEq(expectedCurves));

  Curve i2 {"i2", {4, 5, 6}};
  cs.addIndependentCurve(i2);
  expectedCurves.insert(std::make_pair(i2.getName(), i2));
  EXPECT_THAT(cs.getIndependentCurves(), ContainerEq(expectedCurves));
}

TEST(CurveSet, addIndependentCurves_replaceExisting)
{
  CurveSet cs {"testSet"};

  Curve i1 {"theName", {1, 2, 3}};
  cs.addIndependentCurve(i1);
  EXPECT_THAT(cs.getIndependentCurves().at("theName").getValues(), ElementsAre(1, 2, 3));

  Curve i2 {"theName", {4, 5, 6}};
  cs.addIndependentCurve(i2);
  EXPECT_THAT(cs.getIndependentCurves().at("theName").getValues(), ElementsAre(4, 5, 6));
}

TEST(CurveSet, addDependentCurves)
{
  CurveSet cs {"testSet"};
  std::unordered_map<std::string, Curve> expectedCurves;

  Curve i1 {"i1", {1, 2, 3}};
  cs.addDependentCurve(i1);
  expectedCurves.insert(std::make_pair(i1.getName(), i1));
  EXPECT_THAT(cs.getDependentCurves(), ContainerEq(expectedCurves));

  Curve i2 {"i2", {4, 5, 6}};
  cs.addDependentCurve(i2);
  expectedCurves.insert(std::make_pair(i2.getName(), i2));
  EXPECT_THAT(cs.getDependentCurves(), ContainerEq(expectedCurves));
}

TEST(CurveSet, addDependentCurves_replaceExisting)
{
  CurveSet cs {"testSet"};

  Curve d1 {"theName", {1, 2, 3}};
  cs.addDependentCurve(d1);
  EXPECT_THAT(cs.getDependentCurves().at("theName").getValues(), ElementsAre(1, 2, 3));

  Curve d2 {"theName", {4, 5, 6}};
  cs.addDependentCurve(d2);
  EXPECT_THAT(cs.getDependentCurves().at("theName").getValues(), ElementsAre(4, 5, 6));
}

TEST(CurveSet, createFromNode_empty)
{
  conduit::Node curveSetAsNode = parseJsonValue(R"({})");
  CurveSet curveSet {"theName", curveSetAsNode};
  EXPECT_EQ("theName", curveSet.getName());
  std::unordered_map<std::string, Curve> emptyMap;
  EXPECT_THAT(curveSet.getDependentCurves(), ContainerEq(emptyMap));
  EXPECT_THAT(curveSet.getIndependentCurves(), ContainerEq(emptyMap));
}

TEST(CurveSet, createFromNode_emptySets)
{
  conduit::Node curveSetAsNode = parseJsonValue(R"({
      "dependent": {},
      "independent": {}
    })");
  CurveSet curveSet {"theName", curveSetAsNode};
  EXPECT_EQ("theName", curveSet.getName());
  std::unordered_map<std::string, Curve> emptyMap;
  EXPECT_THAT(curveSet.getDependentCurves(), ContainerEq(emptyMap));
  EXPECT_THAT(curveSet.getIndependentCurves(), ContainerEq(emptyMap));
}

TEST(CurveSet, createFromNode_curveSetsDefined)
{
  conduit::Node curveSetAsNode = parseJsonValue(R"({
      "independent": {
        "indep1": { "value": [10, 20, 30]},
        "indep2/with/slash": { "value": [40, 50, 60]}
      },
      "dependent": {
        "dep1": { "value": [1, 2, 3]},
        "dep2/with/slash": { "value": [4, 5, 6]}
      }
    })");
  CurveSet curveSet {"theName", curveSetAsNode};
  EXPECT_EQ("theName", curveSet.getName());

  std::unordered_map<std::string, Curve> expectedDependents {
    {"dep1", Curve {"dep1", {1, 2, 3}}},
    {"dep2/with/slash", Curve {"dep2/with/slash", {4, 5, 6}}},
  };
  EXPECT_THAT(curveSet.getDependentCurves(), ContainerEq(expectedDependents));

  std::unordered_map<std::string, Curve> expectedIndependents {
    {"indep1", Curve {"indep1", {10, 20, 30}}},
    {"indep2/with/slash", Curve {"indep2/with/slash", {40, 50, 60}}},
  };
  EXPECT_THAT(curveSet.getIndependentCurves(), ContainerEq(expectedIndependents));
}

TEST(CurveSet, createFromNode_orderPreserved)
{
  conduit::Node curveSetAsNode = parseJsonValue(R"({
      "independent": {
        "blue": { "value": [1, 2, 3]},
        "purple": { "value": [1, 2, 3]},
        "green": { "value": [1, 2, 3]},
        "gold": { "value": [1, 2, 3]},
        "red": { "value": [1, 2, 3]}
      },
      "dependent": {
        "purple": { "value": [1, 2, 3]},
        "red": { "value": [1, 2, 3]},
        "blue": { "value": [1, 2, 3]},
        "gold": { "value": [1, 2, 3]},
        "green": { "value": [1, 2, 3]}
      }
    })");
  CurveSet curveSet {"theName", curveSetAsNode};

  std::vector<std::string> expectedDependentsOrder {"purple", "red", "blue", "gold", "green"};

  EXPECT_THAT(curveSet.getOrderedDependentCurveNames(), ContainerEq(expectedDependentsOrder));

  std::vector<std::string> expectedIndependentsOrder {"blue", "purple", "green", "gold", "red"};

  EXPECT_THAT(curveSet.getOrderedIndependentCurveNames(), ContainerEq(expectedIndependentsOrder));
}

TEST(CurveSet, toNode_empty)
{
  CurveSet curveSet {"theName"};
  std::string expected = R"({
        "independent": {},
        "dependent": {}
    })";
  EXPECT_THAT(curveSet.toNode(CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST),
              MatchesJsonMatcher(expected));
}

TEST(CurveSet, toNode_withCurves)
{
  CurveSet curveSet {"theName"};
  curveSet.addIndependentCurve(Curve {"i1", {1, 2, 3}});
  curveSet.addIndependentCurve(Curve {"i2/with/slash", {4, 5, 6}});
  curveSet.addDependentCurve(Curve {"d1", {10, 20, 30}});
  curveSet.addDependentCurve(Curve {"d2/with/slash", {40, 50, 60}});
  std::string expected = R"({
        "independent": {
            "i1": {
                "value": [1.0, 2.0, 3.0]
            },
            "i2/with/slash": {
                "value": [4.0, 5.0, 6.0]
            }
        },
        "dependent": {
            "d1": {
                "value": [10.0, 20.0, 30.0]
            },
            "d2/with/slash": {
                "value": [40.0, 50.0, 60.0]
            }
        }
    })";
  EXPECT_THAT(curveSet.toNode(CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST),
              MatchesJsonMatcher(expected));
}

TEST(CurveSet, toNode_withCurves_orderIsPreserved)
{
  CurveSet curveSet {"theName"};
  curveSet.addIndependentCurve(Curve {"blue", {1, 2, 3}});
  curveSet.addIndependentCurve(Curve {"purple", {1, 2, 3}});
  curveSet.addIndependentCurve(Curve {"green", {1, 2, 3}});
  curveSet.addIndependentCurve(Curve {"gold", {1, 2, 3}});
  curveSet.addIndependentCurve(Curve {"red", {1, 2, 3}});
  curveSet.addDependentCurve(Curve {"purple", {10, 20, 30}});
  curveSet.addDependentCurve(Curve {"red", {10, 20, 30}});
  curveSet.addDependentCurve(Curve {"blue", {10, 20, 30}});
  curveSet.addDependentCurve(Curve {"gold", {10, 20, 30}});
  curveSet.addDependentCurve(Curve {"green", {10, 20, 30}});
  auto expected = R"({
        "independent": {
            "blue": { "value": [1.0, 2.0, 3.0] },
            "purple": { "value": [1.0, 2.0, 3.0] },
            "green": { "value": [1.0, 2.0, 3.0] },
            "gold": { "value": [1.0, 2.0, 3.0] },
            "red": { "value": [1.0, 2.0, 3.0] }
        },
        "dependent": {
            "purple": { "value": [10.0, 20.0, 30.0] },
            "red": { "value": [10.0, 20.0, 30.0] },
            "blue": { "value": [10.0, 20.0, 30.0] },
            "gold": { "value": [10.0, 20.0, 30.0] },
            "green": { "value": [10.0, 20.0, 30.0] }
        }
    })";
  EXPECT_THAT(curveSet.toNode(CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST),
              MatchesJsonMatcher(expected));

  // The order needs to be preserved as we hand it back and forth
  std::vector<std::string> expectedIndependentsOrder {"blue", "purple", "green", "gold", "red"};

  EXPECT_THAT(curveSet.getOrderedIndependentCurveNames(), ContainerEq(expectedIndependentsOrder));

  CurveSet curveSetBack {"theName", curveSet.toNode(CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST)};

  std::vector<std::string> expectedDependentsOrder {"purple", "red", "blue", "gold", "green"};

  EXPECT_THAT(curveSet.getOrderedDependentCurveNames(), ContainerEq(expectedDependentsOrder));
}

TEST(CurveSet, toNode_withCurves_sortOrders)
{
  CurveSet curveSet {"theName"};
  curveSet.addIndependentCurve(Curve {"black", {1, 2, 3}});
  curveSet.addIndependentCurve(Curve {"white", {1, 2, 3}});
  curveSet.addIndependentCurve(Curve {"lime", {1, 2, 3}});
  auto expectedAlphabetic = R"({
        "independent": {
            "black": { "value": [1.0, 2.0, 3.0] },
            "lime": { "value": [1.0, 2.0, 3.0] },
            "white": { "value": [1.0, 2.0, 3.0] }
    }, "dependent": {}})";
  EXPECT_THAT(curveSet.toNode(CurveSet::CurveOrder::ALPHABETIC),
              MatchesJsonMatcher(expectedAlphabetic));
  auto expectedReverseAlphabetic = R"({
        "independent": {
            "white": { "value": [1.0, 2.0, 3.0] },
            "lime": { "value": [1.0, 2.0, 3.0] },
            "black": { "value": [1.0, 2.0, 3.0] }
    }, "dependent": {}})";
  EXPECT_THAT(curveSet.toNode(CurveSet::CurveOrder::REVERSE_ALPHABETIC),
              MatchesJsonMatcher(expectedReverseAlphabetic));
  auto expectedNewestFirst = R"({
        "independent": {
            "lime": { "value": [1.0, 2.0, 3.0] },
            "white": { "value": [1.0, 2.0, 3.0] },
            "black": { "value": [1.0, 2.0, 3.0] }
    }, "dependent": {}})";
  EXPECT_THAT(curveSet.toNode(CurveSet::CurveOrder::REGISTRATION_NEWEST_FIRST),
              MatchesJsonMatcher(expectedNewestFirst));
  // Note this is the "default" order from the user perspective, but Record's in charge of deciding default behavior
  auto expectedOldestFirst = R"({
        "independent": {
            "black": { "value": [1.0, 2.0, 3.0] },
            "white": { "value": [1.0, 2.0, 3.0] },
            "lime": { "value": [1.0, 2.0, 3.0] }
    }, "dependent": {}})";
  EXPECT_THAT(curveSet.toNode(CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST),
              MatchesJsonMatcher(expectedOldestFirst));
}

TEST(CurveSet, customSortOrderReturnsSuccess)
{
  CurveSet curveSet {"theName"};
  curveSet.addDependentCurve(Curve {"puce", {1, 2, 3}});
  curveSet.addDependentCurve(Curve {"blush", {4, 5, 6}});
  curveSet.addDependentCurve(Curve {"seaglass", {10, 11, 12}});
  bool success;
  std::vector<std::string> initialOrder = {"puce", "blush", "seaglass"};
  std::vector<std::string> tooShort = {"blush", "seaglass"};
  std::vector<std::string> repeats = {"blush", "blush", "seaglass"};
  std::vector<std::string> wrongNames = {"tealy", "puce", "blush", "seaglass"};
  std::vector<std::string> justRight = {"seaglass", "puce", "blush"};
  success = curveSet.applyCustomDependentCurveOrder(tooShort);
  EXPECT_THAT(curveSet.getOrderedDependentCurveNames(), ContainerEq(initialOrder));
  EXPECT_THAT(success, false);
  success = curveSet.applyCustomDependentCurveOrder(repeats);
  EXPECT_THAT(curveSet.getOrderedDependentCurveNames(), ContainerEq(initialOrder));
  EXPECT_THAT(success, false);
  success = curveSet.applyCustomDependentCurveOrder(wrongNames);
  EXPECT_THAT(curveSet.getOrderedDependentCurveNames(), ContainerEq(initialOrder));
  EXPECT_THAT(success, false);
  success = curveSet.applyCustomDependentCurveOrder(justRight);
  EXPECT_THAT(curveSet.getOrderedDependentCurveNames(), ContainerEq(justRight));
  EXPECT_THAT(success, true);
}

TEST(CurveSet, customDependentSortOrder)
{
  CurveSet curveSet {"theName"};
  curveSet.addDependentCurve(Curve {"teal", {1, 2, 3}});
  curveSet.addDependentCurve(Curve {"cyan", {4, 5, 6}});
  curveSet.addDependentCurve(Curve {"magenta", {10, 11, 12}});
  curveSet.addDependentCurve(Curve {"salmon", {7, 8, 9}});
  std::vector<std::string> initialOrder = {"teal", "cyan", "magenta", "salmon"};
  std::vector<std::string> newOrder = {"teal", "cyan", "salmon", "magenta"};
  EXPECT_THAT(curveSet.getOrderedDependentCurveNames(), ContainerEq(initialOrder));
  curveSet.applyCustomDependentCurveOrder(newOrder);
  EXPECT_THAT(curveSet.getOrderedDependentCurveNames(), ContainerEq(newOrder));
}

TEST(CurveSet, customIndependentSortOrder)
{
  CurveSet curveSet {"theName"};
  curveSet.addIndependentCurve(Curve {"brown", {-2, -1, 0}});
  curveSet.addIndependentCurve(Curve {"lightgrey", {4, 5, 6}});
  curveSet.addIndependentCurve(Curve {"bronze", {1, 2, 3}});
  std::vector<std::string> initialOrder = {"brown", "lightgrey", "bronze"};
  std::vector<std::string> newOrder = {"bronze", "brown", "lightgrey"};
  EXPECT_THAT(curveSet.getOrderedIndependentCurveNames(), ContainerEq(initialOrder));
  curveSet.applyCustomIndependentCurveOrder(newOrder);
  EXPECT_THAT(curveSet.getOrderedIndependentCurveNames(), ContainerEq(newOrder));
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
