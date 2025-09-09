// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <cstring>
#include <string>
#include <vector>

#include "axom/lumberjack/Message.hpp"

struct TestData
{
  std::string text;
  int rank;
  int lineNumber;
  int level;
  double creationTime;
  std::string tag;
  std::string fileName;
  std::string packed;

  TestData(std::string text_,
           int rank_,
           int lineNumber_,
           int level_,
           double creationTime_,
           std::string tag_,
           std::string fileName_,
           std::string packed_)
    : text(text_)
    , rank(rank_)
    , lineNumber(lineNumber_)
    , level(level_)
    , creationTime(creationTime_)
    , tag(tag_)
    , fileName(fileName_)
    , packed(packed_)
  { }
};

inline std::vector<TestData> getTestData()
{
  std::vector<TestData> testStrings;
  //<ranks delimited by ,>*<rank count>*<file name>*
  //  <line number>*<level>*<creation time>*<tag>*<text>
  testStrings.emplace_back("test empty filename",
                           123,
                           5,
                           0,
                           0,
                           "tag",
                           "",
                           "123*1**5*0*0.000000*tag*test empty filename");
  testStrings.emplace_back("",
                           123,
                           5,
                           1,
                           0,
                           "tag",
                           "test empty message",
                           "123*1*test empty message*5*1*0.000000*tag*");
  testStrings.emplace_back("test",
                           123,
                           5,
                           1,
                           0,
                           "",
                           "test tag message",
                           "123*1*test tag message*5*1*0.000000**test");
  testStrings.emplace_back("test",
                           123,
                           5,
                           1,
                           0,
                           "tag123",
                           "foo.cpp",
                           "123*1*foo.cpp*5*1*0.000000*tag123*test");
  testStrings.emplace_back("abcdef",
                           1,
                           164,
                           1,
                           0,
                           "123tag",
                           "bar/baz.cpp",
                           "1*1*bar/baz.cpp*164*1*0.000000*123tag*abcdef");
  testStrings.emplace_back("this is a test string",
                           0,
                           999999,
                           3,
                           0,
                           "&^Hg",
                           "/test/path.hpp",
                           "0*1*/test/path.hpp*999999*3*0.000000*&^Hg*this is a test string");
  testStrings.emplace_back("123456789",
                           55555,
                           6543,
                           1,
                           0.0,
                           "tag1",
                           "really_long_obnoxious_path-with-lots^stuff.f90",
                           "55555*1*really_long_obnoxious_path-with-lots^stuff."
                           "f90*6543*1*0.000000*tag1*123456789");
  testStrings.emplace_back("                         asdf               ",
                           654987,
                           1,
                           2,
                           0,
                           "mytag",
                           "notimportantfilename",
                           "654987*1*notimportantfilename*1*2*0.000000*mytag*           "
                           "              asdf               ");
  testStrings.emplace_back("//comment",
                           158794,
                           9090,
                           4,
                           0,
                           "tag12",
                           "running out of ideas.cpp",
                           "158794*1*running out of ideas.cpp*9090*4*0.000000*tag12*//comment");
  testStrings.emplace_back("/* test string */",
                           12,
                           876543,
                           1,
                           0,
                           "tag",
                           "234234234.file",
                           "12*1*234234234.file*876543*1*0.000000*tag*/* test string */");
  testStrings.emplace_back(
    "~!@#$%^&*()_+{}|:\"<>?,./;'[]\\-='",
    12,
    654987,
    1,
    0,
    "tag",
    "filenameyepp",
    "12*1*filenameyepp*654987*1*0.000000*tag*~!@#$%^&*()_+{}|:\"<>?,./;'[]\\-='");
  return testStrings;
}

TEST(lumberjack_Message, getSet)
{
  std::vector<TestData> testData = getTestData();

  for(auto& td : testData)
  {
    axom::lumberjack::Message* m = new axom::lumberjack::Message(td.text,
                                                                 td.rank,
                                                                 td.fileName,
                                                                 td.lineNumber,
                                                                 td.level,
                                                                 td.creationTime,
                                                                 td.tag);

    EXPECT_EQ(m->text(), td.text);
    EXPECT_EQ(m->ranks().size(), (std::vector<int>::size_type)1);
    EXPECT_EQ(m->count(), 1);
    EXPECT_EQ(m->ranks()[0], td.rank);
    EXPECT_EQ(m->level(), td.level);
    EXPECT_EQ(m->tag(), td.tag);
    EXPECT_EQ(m->lineNumber(), td.lineNumber);
    EXPECT_EQ(m->fileName(), td.fileName);

    delete m;
  }
}

TEST(lumberjack_Message, getSetCaseConstCharToString)
{
  //Test that const char* will convert fine to string
  const char* textConstCharPointer = "Testing if const char* will convert fine.";
  std::string textString = "Testing if const char* will convert fine.";
  axom::lumberjack::Message m;
  m.text(textConstCharPointer);
  m.addRank(14, 5);

  EXPECT_EQ(m.text(), textString);
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
  EXPECT_EQ(m.count(), 1);
  EXPECT_EQ(m.ranks()[0], 14);
}

TEST(lumberjack_Message, getSetFillRankLimit)
{
  //Test case: one text, filled ranks to rank limit
  const int ranksLimit = 5;
  axom::lumberjack::Message m;
  m.text("Testing filling rank to rank limit.");
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    m.addRank(i + 1, ranksLimit);
  }

  EXPECT_EQ(m.text(), "Testing filling rank to rank limit.");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), ranksLimit);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i + 1);
  }
}

TEST(lumberjack_Message, getSetFillPastRankLimit)
{
  //Test case: one text, filled ranks to past rank limit
  const int ranksLimit = 5;
  axom::lumberjack::Message m;
  m.text("Test filling past rank limit.");
  for(int i = 0; i < (int)ranksLimit * 2; ++i)
  {
    m.addRank(i + 1, ranksLimit);
  }

  EXPECT_EQ(m.text(), "Test filling past rank limit.");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), (int)ranksLimit * 2);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i + 1);
  }
}

TEST(lumberjack_Message, getSetTestAddingVectorOfOneRank)
{
  //Test case: one text, fill ranks with vector of 1 rank
  const int ranksLimit = 5;
  std::vector<int> ranks;
  ranks.push_back(123);

  axom::lumberjack::Message m;
  m.text("Test adding vector of 1 rank.");
  m.addRanks(ranks, 1, ranksLimit);

  EXPECT_EQ(m.text(), "Test adding vector of 1 rank.");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
  EXPECT_EQ(m.count(), 1);
  EXPECT_EQ(m.ranks()[0], 123);
}

TEST(lumberjack_Message, getSetTestAddingVectorOfManyRanks)
{
  //Test case: one text, fill ranks with vector of many ranks
  const int ranksLimit = 5;
  const int ranksAdded = 10;
  std::vector<int> ranks;
  for(int i = 0; i < (int)ranksAdded; ++i)
  {
    ranks.push_back(i);
  }

  axom::lumberjack::Message m;
  m.text("Test adding vector of many ranks over limit.");
  m.addRanks(ranks, ranksAdded, ranksLimit);

  EXPECT_EQ(m.text(), "Test adding vector of many ranks over limit.");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), ranksAdded);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i);
  }
}

TEST(lumberjack_Message, getSetTestAddingExactNumberOfRanksLimit)
{
  //Test case: one text, fill ranks with vector of ranks don't go past
  // ranksLimit
  const int ranksLimit = 5;
  std::vector<int> ranks;
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    ranks.push_back(i + 1);
  }

  axom::lumberjack::Message m;
  m.text("Test adding vector of exactly ranksLimit of ranks.");
  m.addRanks(ranks, ranksLimit, ranksLimit);

  EXPECT_EQ(m.text(), "Test adding vector of exactly ranksLimit of ranks.");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), ranksLimit);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i + 1);
  }
}

TEST(lumberjack_Message, getSetAddSameRankMultipleTimes)
{
  //Test case: one text, add same rank multiple times to make sure ranks stays 1
  // but count increments
  const int ranksLimit = 5;
  axom::lumberjack::Message m;
  for(int i = 0; i < (int)ranksLimit * 3; ++i)
  {
    m.addRank(1, ranksLimit);
  }

  m.text("Test adding same rank over and over.");

  EXPECT_EQ(m.text(), "Test adding same rank over and over.");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
  EXPECT_EQ(m.count(), (int)ranksLimit * 3);
  EXPECT_EQ(m.ranks()[0], 1);
}

TEST(lumberjack_Message, getSetAddRanksWithOverRankLimit)
{
  //Test case: one text, fill ranks with vector of ranks that will go past
  // ranksLimit
  const int ranksLimit = 5;
  std::vector<int> ranks;
  for(int i = 0; i < (int)ranksLimit * 3; ++i)
  {
    ranks.push_back(1);
  }

  axom::lumberjack::Message m;
  m.text("This message is unimportant.");
  m.addRanks(ranks, ranksLimit * 3, ranksLimit);

  EXPECT_EQ(m.text(), "This message is unimportant.");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
  EXPECT_EQ(m.count(), (int)ranksLimit * 3);
  EXPECT_EQ(m.ranks()[0], 1);
}

TEST(lumberjack_Message, testConstructor01)
{
  //Test most basic case: one text, one rank, file name, line number
  axom::lumberjack::Message
    m("Testing the basic message constructor", 122, "foo.cpp", 154, 1, 0.0, "tag1");

  EXPECT_EQ(m.text(), "Testing the basic message constructor");
  EXPECT_EQ(m.fileName(), "foo.cpp");
  EXPECT_EQ(m.lineNumber(), 154);
  EXPECT_EQ(m.level(), 1);
  EXPECT_EQ(m.tag(), "tag1");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
  EXPECT_EQ(m.count(), 1);
  EXPECT_EQ(m.ranks()[0], 122);
}

TEST(lumberjack_Message, testConstructor02)
{
  //Test most basic case: one text, one rank, file name, line number
  const int ranksLimit = 5;
  std::vector<int> ranks;
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    ranks.push_back(i + 1);
  }

  axom::lumberjack::Message m("test message constructor with vector of ranks",
                              ranks,
                              ranksLimit,
                              ranksLimit,
                              "foo.cpp",
                              154,
                              2,
                              0.0,
                              "mytag");

  EXPECT_EQ(m.text(), "test message constructor with vector of ranks");
  EXPECT_EQ(m.fileName(), "foo.cpp");
  EXPECT_EQ(m.lineNumber(), 154);
  EXPECT_EQ(m.level(), 2);
  EXPECT_EQ(m.tag(), "mytag");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), ranksLimit);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i + 1);
  }
}

TEST(lumberjack_Message, stringOfRanks01)
{
  //Test most basic case: one rank
  const int ranksLimit = 5;
  axom::lumberjack::Message m;
  m.addRank(400, ranksLimit);

  EXPECT_EQ(m.text(), "");
  EXPECT_EQ(m.fileName(), "");
  EXPECT_EQ(m.lineNumber(), 0);
  EXPECT_EQ(m.level(), 0);
  EXPECT_EQ(m.tag(), "");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
  EXPECT_EQ(m.count(), 1);
  EXPECT_EQ(m.ranks()[0], 400);
  EXPECT_EQ(m.stringOfRanks(), "400");
}

TEST(lumberjack_Message, stringOfRanks02)
{
  //Test case: more than one rank
  const int ranksLimit = 5;
  std::vector<int> ranks;
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    ranks.push_back(i + 1);
  }

  axom::lumberjack::Message m;
  m.addRanks(ranks, ranksLimit, ranksLimit);

  EXPECT_EQ(m.text(), "");
  EXPECT_EQ(m.fileName(), "");
  EXPECT_EQ(m.lineNumber(), 0);
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), ranksLimit);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i + 1);
  }
  EXPECT_EQ(m.stringOfRanks(), "1,2,3,4,5...");
}

TEST(lumberjack_Message, stringOfRanks03)
{
  //Test case: full Message
  const int ranksLimit = 5;
  std::vector<int> ranks;
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    ranks.push_back(i * 2);
  }

  axom::lumberjack::Message
    m("Unimportant message", ranks, ranksLimit, ranksLimit, "test/foo.cpp", 987654321, 0, 0.0, "");

  EXPECT_EQ(m.text(), "Unimportant message");
  EXPECT_EQ(m.fileName(), "test/foo.cpp");
  EXPECT_EQ(m.lineNumber(), 987654321);
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), ranksLimit);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i * 2);
  }
  EXPECT_EQ(m.stringOfRanks(), "0,2,4,6,8...");
}

TEST(lumberjack_Message, pack01)
{
  //Test case: full Message
  const int ranksLimit = 5;
  std::vector<int> ranks;
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    ranks.push_back(i * 2);
  }

  axom::lumberjack::Message
    m("Unimportant message", ranks, ranksLimit, ranksLimit, "test/foo.cpp", 987654321, 0, 0.0, "");

  EXPECT_EQ(m.text(), "Unimportant message");
  EXPECT_EQ(m.fileName(), "test/foo.cpp");
  EXPECT_EQ(m.lineNumber(), 987654321);
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), ranksLimit);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i * 2);
  }

  std::string packedMessage = m.pack();
  EXPECT_EQ(packedMessage,
            "0,2,4,6,8*5*test/foo.cpp*987654321*0*" + std::to_string(m.creationTime()) +
              "**Unimportant message");
}

TEST(lumberjack_Message, unpack01)
{
  //Test case: full Message
  const int ranksLimit = 5;

  axom::lumberjack::Message m;
  m.unpack("0,2,4,6,8*15*test/foo.cpp*987654321*1*0*gtest*Unimportant message", ranksLimit);

  EXPECT_EQ(m.text(), "Unimportant message");
  EXPECT_EQ(m.fileName(), "test/foo.cpp");
  EXPECT_EQ(m.lineNumber(), 987654321);
  EXPECT_EQ(m.level(), 1);
  EXPECT_EQ(m.tag(), "gtest");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), 15);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i * 2);
  }
}

TEST(lumberjack_Message, unpack02)
{
  //Test case: semi full Message
  const int ranksLimit = 5;

  axom::lumberjack::Message m;
  m.unpack("0,2,4,6,8*15*test/foo.cpp*987654321*0*0**Unimportant message", ranksLimit);

  EXPECT_EQ(m.text(), "Unimportant message");
  EXPECT_EQ(m.fileName(), "test/foo.cpp");
  EXPECT_EQ(m.lineNumber(), 987654321);
  EXPECT_EQ(m.level(), 0);
  EXPECT_EQ(m.tag(), "");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), 15);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i * 2);
  }
}

TEST(lumberjack_Message, packEmptyMessage)
{
  //Test case: pack empty Message
  axom::lumberjack::Message m("", 1, "test/foo.cpp", 987654321, 0, 0.0, "tag");

  std::string packedMessage = m.pack();
  EXPECT_EQ(packedMessage,
            "1*1*test/foo.cpp*987654321*0*" + std::to_string(m.creationTime()) + "*tag*");
}

TEST(lumberjack_Message, packEmptyTag)
{
  //Test case: pack empty Tag
  axom::lumberjack::Message m("asdfMessageadsf", 1, "test/foo.cpp", 987654321, 0, 0.0, "");

  std::string packedMessage = m.pack();
  EXPECT_EQ(
    packedMessage,
    "1*1*test/foo.cpp*987654321*0*" + std::to_string(m.creationTime()) + "**asdfMessageadsf");
}

TEST(lumberjack_Message, packEmptyTagAndMessage)
{
  //Test case: pack empty Tag and Message
  axom::lumberjack::Message m("", 1, "test/foo.cpp", 987654321, 0, 0.0, "");

  std::string packedMessage = m.pack();
  EXPECT_EQ(packedMessage,
            "1*1*test/foo.cpp*987654321*0*" + std::to_string(m.creationTime()) + "**");
}

TEST(lumberjack_Message, unpackEmptyMessage)
{
  //Test case: unpack empty Message
  const int ranksLimit = 5;

  axom::lumberjack::Message m;
  m.unpack("0,2,4,6,8*15*test/foo.cpp*987654321*0*" + std::to_string(m.creationTime()) + "*tag*",
           ranksLimit);

  EXPECT_EQ(m.text(), "");
  EXPECT_EQ(m.fileName(), "test/foo.cpp");
  EXPECT_EQ(m.lineNumber(), 987654321);
  EXPECT_EQ(m.level(), 0);
  EXPECT_EQ(m.tag(), "tag");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), 15);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i * 2);
  }
}

TEST(lumberjack_Message, unpackEmptyTag)
{
  //Test case: unpack empty Tag
  const int ranksLimit = 5;

  axom::lumberjack::Message m;
  m.unpack("0,2,4,6,8*15*test/foo.cpp*987654321*0*0**123unimportant123", ranksLimit);

  EXPECT_EQ(m.text(), "123unimportant123");
  EXPECT_EQ(m.fileName(), "test/foo.cpp");
  EXPECT_EQ(m.lineNumber(), 987654321);
  EXPECT_EQ(m.level(), 0);
  EXPECT_EQ(m.tag(), "");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), 15);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i * 2);
  }
}

TEST(lumberjack_Message, unpackEmptyTagAndMessage)
{
  //Test case: unpack empty Tag and Message
  const int ranksLimit = 5;

  axom::lumberjack::Message m;
  m.unpack("0,2,4,6,8*15*test/foo.cpp*987654321*0*0**", ranksLimit);

  EXPECT_EQ(m.text(), "");
  EXPECT_EQ(m.fileName(), "test/foo.cpp");
  EXPECT_EQ(m.lineNumber(), 987654321);
  EXPECT_EQ(m.level(), 0);
  EXPECT_EQ(m.tag(), "");
  EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)ranksLimit);
  EXPECT_EQ(m.count(), 15);
  for(int i = 0; i < (int)ranksLimit; ++i)
  {
    EXPECT_EQ(m.ranks()[i], i * 2);
  }
}

TEST(lumberjack_Message, packMessagesIndividually)
{
  std::vector<axom::lumberjack::Message*> messages;
  std::vector<TestData> testData = getTestData();
  for(auto& td : testData)
  {
    axom::lumberjack::Message* m = new axom::lumberjack::Message(td.text,
                                                                 td.rank,
                                                                 td.fileName,
                                                                 td.lineNumber,
                                                                 td.level,
                                                                 td.creationTime,
                                                                 td.tag);
    messages.push_back(m);

    const char* packedMessage = axom::lumberjack::packMessages(messages);

    std::string answer = "1*" + std::to_string((int)td.packed.length()) + "*" + td.packed;

    EXPECT_EQ(packedMessage, answer);

    // cleanup
    delete[] packedMessage;
    delete m;
    messages.clear();
  }
}

TEST(lumberjack_Message, packMessages)
{
  std::vector<axom::lumberjack::Message*> messages;
  std::vector<TestData> testData = getTestData();
  std::string answer = "";
  for(auto& td : testData)
  {
    axom::lumberjack::Message* m = new axom::lumberjack::Message(td.text,
                                                                 td.rank,
                                                                 td.fileName,
                                                                 td.lineNumber,
                                                                 td.level,
                                                                 td.creationTime,
                                                                 td.tag);
    messages.push_back(m);

    answer += std::to_string((int)td.packed.length()) + "*" + td.packed;
  }

  const char* packedMessages = axom::lumberjack::packMessages(messages);

  answer = std::to_string((int)testData.size()) + "*" + answer;
  EXPECT_EQ(packedMessages, answer);

  // cleanup
  delete[] packedMessages;
  for(auto* _m : messages)
  {
    delete _m;
  }
  messages.clear();
}

TEST(lumberjack_Message, unpackMessagesIndividually)
{
  std::vector<axom::lumberjack::Message*> messages;
  std::vector<TestData> testData = getTestData();
  for(auto& td : testData)
  {
    std::string answer = "1*" + std::to_string((int)td.packed.length()) + "*" + td.packed;
    axom::lumberjack::unpackMessages(messages, answer.c_str(), 100);

    axom::lumberjack::Message* m = messages[0];
    EXPECT_EQ(m->text(), td.text);
    EXPECT_EQ(m->ranks()[0], td.rank);
    EXPECT_EQ(m->fileName(), td.fileName);
    EXPECT_EQ(m->lineNumber(), td.lineNumber);
    EXPECT_EQ(m->level(), td.level);
    EXPECT_EQ(m->creationTime(), td.creationTime);
    EXPECT_EQ(m->tag(), td.tag);

    delete m;
    messages.clear();
  }
}

TEST(lumberjack_Message, unpackMessages)
{
  std::vector<axom::lumberjack::Message*> messages;
  std::vector<TestData> testData = getTestData();
  std::string packedMessages = std::to_string((int)testData.size()) + "*";
  for(auto& td : testData)
  {
    packedMessages += std::to_string((int)td.packed.length()) + "*" + td.packed;
  }

  axom::lumberjack::unpackMessages(messages, packedMessages.c_str(), 100);

  for(int i = 0; i < (int)messages.size(); ++i)
  {
    axom::lumberjack::Message* m = messages[i];
    TestData td = testData[i];

    EXPECT_EQ(m->text(), td.text);
    EXPECT_EQ(m->ranks()[0], td.rank);
    EXPECT_EQ(m->fileName(), td.fileName);
    EXPECT_EQ(m->lineNumber(), td.lineNumber);
    EXPECT_EQ(m->level(), td.level);
    EXPECT_EQ(m->tag(), td.tag);

    delete m;
  }
  messages.clear();
}
