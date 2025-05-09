// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/lumberjack/TextEqualityCombiner.hpp"

TEST(lumberjack_TextEqualityCombiner, case01)
{
  //Test most basic case: two equal texts
  std::string text = "I never wanted to do this job in the first place!";
  axom::lumberjack::Message m1;
  m1.text(text);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);

  axom::lumberjack::Message m2;
  m2.text(text);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);

  axom::lumberjack::TextEqualityCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  c.combine(m1, m2, 5);

  EXPECT_EQ(shouldMessagesBeCombined, true);
  EXPECT_EQ(m1.text().compare(text), 0);
  EXPECT_EQ(m1.count(), 2);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.ranks()[1], 14);
}

TEST(lumberjack_TextEqualityCombiner, case02)
{
  //Test most basic negative case: two not-equal texts
  std::string text1 = "I never wanted to do this job in the first place!";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);

  std::string text2 = "This text is not equal to the first.";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);

  axom::lumberjack::TextEqualityCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  EXPECT_EQ(shouldMessagesBeCombined, false);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);

  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
}
