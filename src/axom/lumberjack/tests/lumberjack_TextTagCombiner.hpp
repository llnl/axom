// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/lumberjack/TextTagCombiner.hpp"

TEST(lumberjack_TextTagCombiner, case01)
{
  //Test positive case: two equal texts and tags
  std::string text = "I never wanted to do this job in the first place!";
  axom::lumberjack::Message m1;
  m1.text(text);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);
  m1.tag("myTag");

  axom::lumberjack::Message m2;
  m2.text(text);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);
  m2.tag("myTag");

  axom::lumberjack::TextTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  c.combine(m1, m2, 5);

  EXPECT_EQ(shouldMessagesBeCombined, true);
  EXPECT_EQ(m1.text().compare(text), 0);
  EXPECT_EQ(m1.count(), 2);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.ranks()[1], 14);
  EXPECT_EQ(m2.tag().compare("myTag"), 0);
}

TEST(lumberjack_TextTagCombiner, case02)
{
  //Test negative case: two not-equal texts, but equal tags
  std::string text1 = "I never wanted to do this job in the first place!";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);
  m1.tag("myTag");

  std::string text2 = "This text is not equal to the first.";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);
  m2.tag("myTag");

  axom::lumberjack::TextTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  EXPECT_EQ(shouldMessagesBeCombined, false);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myTag"), 0);
}

TEST(lumberjack_TextTagCombiner, case03)
{
  //Test negative case: two equal texts, but not-equal tags
  std::string text = "I never wanted to do this job in the first place!";

  axom::lumberjack::Message m1;
  m1.text(text);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);
  m1.tag("tag1");

  axom::lumberjack::Message m2;
  m2.text(text);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);
  m2.tag("tag2");

  axom::lumberjack::TextTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  EXPECT_EQ(shouldMessagesBeCombined, false);
  EXPECT_EQ(m1.text().compare(text), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("tag1"), 0);

  EXPECT_EQ(m2.text().compare(text), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("tag2"), 0);
}
