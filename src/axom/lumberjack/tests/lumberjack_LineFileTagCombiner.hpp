#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "axom/lumberjack/Message.hpp"
#include "axom/lumberjack/LineFileTagCombiner.hpp"

// Helper function to create a Message
inline axom::lumberjack::Message createMessage(const std::string& text,
                                               int rank,
                                               int rankCount,
                                               const std::string& fileName,
                                               int lineNumber,
                                               const std::string& tag)
{
  axom::lumberjack::Message message;
  message.text(text);
  message.addRank(rank, rankCount);
  message.fileName(fileName);
  message.lineNumber(lineNumber);
  message.tag(tag);
  return message;
}

// Helper function to verify Message attributes
inline void verifyMessage(const axom::lumberjack::Message& message,
                          const std::string& expectedText,
                          int expectedLineNumber,
                          const std::string& expectedFileName,
                          int expectedCount,
                          const std::vector<int>& expectedRanks,
                          const std::string& expectedTag)
{
  EXPECT_EQ(message.text().compare(expectedText), 0);
  EXPECT_EQ(message.lineNumber(), expectedLineNumber);
  EXPECT_EQ(message.fileName().compare(expectedFileName), 0);
  EXPECT_EQ(message.count(), expectedCount);
  EXPECT_EQ(message.ranks().size(), expectedRanks.size());
  for(size_t i = 0; i < expectedRanks.size(); ++i)
  {
    EXPECT_EQ(message.ranks()[i], expectedRanks[i]);
  }
  EXPECT_EQ(message.tag().compare(expectedTag), 0);
}

// Struct to hold test parameters
struct TestParams
{
  std::string caseName {""};
  int lineNumber1 {-1};
  int lineNumber2 {-1};
  std::string fileName1 {""};
  std::string fileName2 {""};
  std::string tag1 {""};
  std::string tag2 {""};
  bool shouldCombine {false};
};

// Parameterized Test Fixture
class LumberjackLineFileTagCombinerTest : public ::testing::TestWithParam<TestParams>
{ };

TEST_P(LumberjackLineFileTagCombinerTest, CombineMessages)
{
  // Arrange
  TestParams params = GetParam();
  std::string text1 = "This message does not matter because we do not filter by text";
  std::string text2 = "This message ALSO does not matter because we do not filter by text";

  axom::lumberjack::Message m1 =
    createMessage(text1, 13, 5, params.fileName1, params.lineNumber1, params.tag1);
  axom::lumberjack::Message m2 =
    createMessage(text2, 14, 5, params.fileName2, params.lineNumber2, params.tag2);

  axom::lumberjack::LineFileTagCombiner combiner;

  // Act
  bool shouldMessagesBeCombined = combiner.shouldMessagesBeCombined(m1, m2);
  if(shouldMessagesBeCombined)
  {
    combiner.combine(m1, m2, 5);
  }

  // Assert
  EXPECT_EQ(shouldMessagesBeCombined, params.shouldCombine);
  if(params.shouldCombine)
  {
    verifyMessage(m1, text1, params.lineNumber1, params.fileName1, 2, {13, 14}, params.tag1);
    verifyMessage(m2, text2, params.lineNumber2, params.fileName2, 1, {14}, params.tag2);
  }
  else
  {
    verifyMessage(m1, text1, params.lineNumber1, params.fileName1, 1, {13}, params.tag1);
    verifyMessage(m2, text2, params.lineNumber2, params.fileName2, 1, {14}, params.tag2);
  }
}

// Custom Name Generator for Parameterized Tests
inline std::string CustomNameGenerator(const ::testing::TestParamInfo<TestParams>& info)
{
  const TestParams& params = info.param;
  return params.caseName;
}

// Define test cases
INSTANTIATE_TEST_SUITE_P(
  LumberjackTests,
  LumberjackLineFileTagCombinerTest,
  ::testing::Values(
    // Positive case: Equal line numbers, filenames, and tags
    TestParams {"case1", 154, 154, "foo.cpp", "foo.cpp", "myTag", "myTag", true},

    // Negative case: Different line numbers, same filenames, same tags
    TestParams {"case2", 12000, 154, "foo.cpp", "foo.cpp", "myTag", "myTag", false},

    // Negative case: Same line numbers, different filenames, same tags
    TestParams {"case3", 154, 154, "foo.cpp", "bar.cpp", "myTag", "myTag", false},

    // Negative case: Same line numbers, same filenames, different tags
    TestParams {"case4", 154, 154, "foo.cpp", "foo.cpp", "myTag", "myOtherTag", false},

    // Negative case: Different line numbers, different filenames, same tags
    TestParams {"case5", 12000, 154, "foo.cpp", "bar.cpp", "myTag", "myTag", false},

    // Negative case: Different line numbers, same filenames, different tags
    TestParams {"case6", 12000, 154, "foo.cpp", "foo.cpp", "myTag", "myOtherTag", false},

    // Negative case: Same line numbers, different filenames, different tags
    TestParams {"case7", 154, 154, "foo.cpp", "bar.cpp", "myTag", "myOtherTag", false},

    // Negative case: Different line numbers, different filenames, different tags
    TestParams {"case8", 12000, 154, "foo.cpp", "bar.cpp", "myTag", "myOtherTag", false}),
  CustomNameGenerator  // Use the custom name generator
);