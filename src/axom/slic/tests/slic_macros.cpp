// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// axom includes
#include "axom/config.hpp"
#include "axom/core/utilities/FileUtilities.hpp"

// slic includes
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

// gtest includes
#include "gtest/gtest.h"  // for gtest macros

// C/C++ includes
#include <string>   // for C++ string
#include <sstream>  // for std::ostringstream

// namespace alias
namespace slic = axom::slic;

namespace axom
{
namespace slic
{
namespace internal
{
std::ostringstream test_stream;

bool is_stream_empty() { return test_stream.str().empty(); }

void clear()
{
  test_stream.clear();
  test_stream.str("");
  EXPECT_TRUE(is_stream_empty());
}

}  // namespace internal
}  // namespace slic
}  // namespace axom

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
void check_level(const std::string& msg, const std::string& expected_level)
{
  EXPECT_FALSE(msg.empty());

  // extract level
  size_t start = msg.find("[") + 1;
  size_t end = expected_level.length();
  std::string level = msg.substr(start, end);

  // check with expected
  EXPECT_FALSE(level.empty());
  EXPECT_EQ(level, expected_level);
}

//------------------------------------------------------------------------------
void check_msg(const std::string& msg, const std::string& expected_message)
{
  EXPECT_FALSE(msg.empty());

  // extract message
  size_t start = msg.find(";;") + 2;
  size_t end = expected_message.length();
  std::string mymsg = msg.substr(start, end);

  EXPECT_FALSE(mymsg.empty());
  EXPECT_EQ(mymsg, expected_message);
}

//------------------------------------------------------------------------------
void check_file(const std::string& msg)
{
  EXPECT_FALSE(msg.empty());

  const std::string EXPECTED_FILE = std::string(__FILE__);

  // extract filename
  size_t start = msg.find("@@") + 2;
  size_t end = EXPECTED_FILE.length();
  std::string myfile = msg.substr(start, end);

  EXPECT_FALSE(myfile.empty());
  EXPECT_EQ(myfile, std::string(__FILE__));
}

//------------------------------------------------------------------------------
void check_line(const std::string& msg, int expected_line)
{
  EXPECT_FALSE(msg.empty());

  // extract line
  size_t start = msg.rfind("@@") + 2;
  std::string l = msg.substr(start);
  const int line = std::stoi(l);

  // expected value is the line number of the SLIC message
  EXPECT_EQ(line, expected_line);
}

//------------------------------------------------------------------------------
void check_tag(const std::string& msg, const std::string& expected_tag)
{
  EXPECT_FALSE(msg.empty());

  // extract tag
  size_t start = msg.rfind("##") + 2;
  std::string tag = msg.substr(start, expected_tag.length());
  EXPECT_EQ(tag, expected_tag);
}

//------------------------------------------------------------------------------
int check_count(const std::string& msg, const std::string& expected_level)
{
  EXPECT_FALSE(msg.empty());

  int count = 0;
  for(size_t pos = msg.find(expected_level); pos != std::string::npos;
      pos = msg.find(expected_level, pos + expected_level.size()))
  {
    ++count;
  }
  return count;
}

//------------------------------------------------------------------------------
// Checks level, message, line number, and file location.
// Clears stream when finished.
void check_level_msg_line_file(const std::string& level, const std::string& message, int expected_line)
{
  EXPECT_FALSE(slic::internal::is_stream_empty());
  const std::string str = slic::internal::test_stream.str();

  check_level(str, level);
  check_msg(str, message);
  check_line(str, expected_line);
  check_file(str);
  slic::internal::clear();
}

// Convenience test macro that checks slic macro has logged a message
#define EXPECT_SLIC_LOG(macro_call, level, message)           \
  do                                                          \
  {                                                           \
    const int expected_line = __LINE__;                       \
    macro_call;                                               \
    check_level_msg_line_file(level, message, expected_line); \
  } while(false)

// Convenience test macro that checks SLIC_*_ONCE macro logs one message
#define EXPECT_SLIC_ONCE(macro_call, level, message)                     \
  do                                                                     \
  {                                                                      \
    const int expected_line = __LINE__;                                  \
    for(int i = 0; i < 2; i++)                                           \
    {                                                                    \
      macro_call;                                                        \
    }                                                                    \
    EXPECT_EQ(check_count(slic::internal::test_stream.str(), level), 1); \
    check_level_msg_line_file(level, message, expected_line);            \
  } while(false)

}  // end anonymous namespace

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(slic_macros, test_error_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  EXPECT_SLIC_LOG(SLIC_ERROR("test error message"), "ERROR", "test error message");

  SLIC_ERROR_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_LOG(SLIC_ERROR_IF(true, "this message is logged!"), "ERROR", "this message is logged!");
  // clang-format on

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_ERROR_ROOT_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  SLIC_ERROR_ROOT_IF(true, "this message is logged!");

  // clang-format off
  EXPECT_SLIC_LOG(SLIC_ERROR_ROOT_IF(true, "this message is logged!"), "ERROR", "this message is logged!");
  // clang-format on

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_ERROR_ROOT_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_ERROR_ROOT_IF(true, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_warning_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  EXPECT_SLIC_LOG(SLIC_WARNING("test warning message"), "WARNING", "test warning message");

  // Called once per call site
  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_ONCE(SLIC_WARNING_ONCE("test warning message once"), "WARNING", "test warning message once");
  // clang-format on

  // Two different call sites, will have two messages
  SLIC_WARNING_ONCE("test warning message #1");
  SLIC_WARNING_ONCE("test warning message #2");
  EXPECT_EQ(check_count(slic::internal::test_stream.str(), "WARNING"), 2);
  slic::internal::clear();

  SLIC_WARNING_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_WARNING_IF_ONCE(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // clang-format off
  EXPECT_SLIC_LOG(SLIC_WARNING_IF(true, "this message is logged!"), "WARNING", "this message is logged!");

  EXPECT_SLIC_ONCE(SLIC_WARNING_IF_ONCE(true, "this message is logged once!"), "WARNING", "this message is logged once!");
  // clang-format on

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_WARNING_ROOT_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_WARNING_ROOT_IF_ONCE(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  axom::slic::setIsRoot(true);
  // clang-format off
  EXPECT_SLIC_LOG(SLIC_WARNING_ROOT("this message is logged on root!"), "WARNING", "this message is logged on root!");

  EXPECT_SLIC_ONCE(SLIC_WARNING_ROOT_ONCE("this message is logged on root once!"), "WARNING", "this message is logged on root once!");

  // Check selective filter based on root == true
  EXPECT_SLIC_LOG(SLIC_WARNING_ROOT_IF(true, "this message is logged!"), "WARNING", "this message is logged!");

  EXPECT_SLIC_ONCE(SLIC_WARNING_ROOT_IF_ONCE(true, "this message is logged once!"), "WARNING", "this message is logged once!");
  // clang-format on

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_WARNING_ROOT_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_WARNING_ROOT_IF_ONCE(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_WARNING_ROOT_IF(true, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_WARNING_ROOT_IF_ONCE(true, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_info_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  EXPECT_SLIC_LOG(SLIC_INFO("test info message"), "INFO", "test info message");

  EXPECT_SLIC_ONCE(SLIC_INFO_ONCE("this info message once"), "INFO", "this info message once");

  SLIC_INFO_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_INFO_IF_ONCE(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  EXPECT_SLIC_LOG(SLIC_INFO_IF(true, "this message is logged!"), "INFO", "this message is logged!");

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_ONCE(SLIC_INFO_IF_ONCE(true, "this message is logged once!"), "INFO", "this message is logged once!");
  // clang-format on

  axom::slic::setIsRoot(true);
  // clang-format off
  EXPECT_SLIC_LOG(SLIC_INFO_ROOT("this message is logged on root!"), "INFO", "this message is logged on root!");

  EXPECT_SLIC_ONCE(SLIC_INFO_ROOT_ONCE("this message is logged on root once!"), "INFO", "this message is logged on root once!");
  // clang-format on

  // is root, but conditional is false -> no message
  SLIC_INFO_ROOT_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_INFO_ROOT_IF_ONCE(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_INFO_ROOT_IF(true, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_INFO_ROOT_IF_ONCE(true, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_debug_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
#ifdef AXOM_DEBUG
  EXPECT_SLIC_LOG(SLIC_DEBUG("test debug message"), "DEBUG", "test debug message");

  EXPECT_SLIC_ONCE(SLIC_DEBUG_ONCE("test debug message once"), "DEBUG", "test debug message once");

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_LOG(SLIC_DEBUG_IF(true, "this message is logged!"), "DEBUG", "this message is logged!");

  EXPECT_SLIC_ONCE(SLIC_DEBUG_IF_ONCE(true, "this message is logged once!"), "DEBUG", "this message is logged once!");
  // clang-format on

  axom::slic::setIsRoot(true);
  EXPECT_SLIC_LOG(SLIC_DEBUG_ROOT("this message is logged!"), "DEBUG", "this message is logged!");

  // clang-format off
  EXPECT_SLIC_ONCE(SLIC_DEBUG_ROOT_ONCE("this message is logged once!"), "DEBUG", "this message is logged once!");

  // Check selective filter based on root == true
  EXPECT_SLIC_LOG(SLIC_DEBUG_ROOT_IF(true, "this message is logged!"), "DEBUG", "this message is logged!");

  EXPECT_SLIC_ONCE(SLIC_DEBUG_ROOT_IF_ONCE(true, "this message is logged once!"), "DEBUG", "this message is logged once!");
  // clang-format on

#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined

  SLIC_DEBUG("test debug message");
  SLIC_DEBUG_ONCE("test debug message");

  SLIC_DEBUG_IF(true, "this message is logged!");
  SLIC_DEBUG_IF_ONCE(true, "this message is logged!");

  axom::slic::setIsRoot(true);
  SLIC_DEBUG_ROOT("this message is logged!");
  SLIC_DEBUG_ROOT_ONCE("this message is logged!");

  SLIC_DEBUG_ROOT_IF(true, "this message is logged!");
  SLIC_DEBUG_ROOT_IF_ONCE(true, "this message is logged!");

  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_DEBUG_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_DEBUG_IF_ONCE(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_DEBUG_ROOT_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_DEBUG_ROOT_IF_ONCE(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_DEBUG_ROOT_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_DEBUG_ROOT_IF_ONCE(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_DEBUG_ROOT_IF(true, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_DEBUG_ROOT_IF_ONCE(true, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_assert_macros)
{
  [[maybe_unused]] int expected_line_number;
  EXPECT_TRUE(slic::internal::is_stream_empty());

  constexpr int val = 42;
  SLIC_ASSERT(val < 0);
  expected_line_number = __LINE__ - 1;
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  check_level_msg_line_file("ERROR", "Failed Assert: val < 0", expected_line_number);
#else
  // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  AXOM_UNUSED_VAR(expected_line_number);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_ASSERT(val > 0);
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_ASSERT_MSG(val < 0, "val should be negative!");
  expected_line_number = __LINE__ - 1;
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  check_level_msg_line_file("ERROR",
                            "Failed Assert: val < 0\nval should be negative!",
                            expected_line_number);
#else
  // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  AXOM_UNUSED_VAR(expected_line_number);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_check_macros)
{
  [[maybe_unused]] int expected_line_number;
  EXPECT_TRUE(slic::internal::is_stream_empty());

  constexpr int val = 42;
  SLIC_CHECK(val < 0);
  expected_line_number = __LINE__ - 1;
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  check_level_msg_line_file("WARNING", "Failed Check: val < 0", expected_line_number);
#else
  // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  AXOM_UNUSED_VAR(expected_line_number);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_CHECK(val > 0);
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_CHECK_MSG(val < 0, "val should be negative!");
  expected_line_number = __LINE__ - 1;
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  check_level_msg_line_file("WARNING",
                            "Failed Check: val < 0\nval should be negative!",
                            expected_line_number);
#else
  // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  AXOM_UNUSED_VAR(expected_line_number);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_tagged_macros)
{
  int expected_line_number;

  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_INFO_TAGGED("test tagged info message", "myTag");
  expected_line_number = __LINE__ - 1;

  check_tag(slic::internal::test_stream.str(), "myTag");
  check_level_msg_line_file("INFO", "test tagged info message", expected_line_number);

  for(int i = 0; i < 2; i++)
  {
    SLIC_INFO_TAGGED_ONCE("test tagged info message once", "myTag");
  }
  expected_line_number = __LINE__ - 2;

  EXPECT_EQ(check_count(slic::internal::test_stream.str(), "INFO"), 1);
  check_tag(slic::internal::test_stream.str(), "myTag");
  check_level_msg_line_file("INFO", "test tagged info message once", expected_line_number);

  SLIC_INFO_TAGGED("this message should not be logged (no tag given)!", "");
  SLIC_INFO_TAGGED_ONCE("this message should not be logged (no tag given)!", "");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_INFO_TAGGED("this message should not be logged (tag DNE)!", "tag404");
  SLIC_INFO_TAGGED_ONCE("this message should not be logged (tag DNE)!", "tag404");
  EXPECT_TRUE(slic::internal::is_stream_empty());
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_if_once_macros)
{
  // Check that message is logged when condition is satisfied only once
  for(int i = 0; i < 3; i++)
  {
    SLIC_INFO_IF_ONCE(i > 0, i << "th message is logged!");
  }
  int expected_line_number = __LINE__ - 2;
  EXPECT_EQ(check_count(slic::internal::test_stream.str(), "INFO"), 1);
  check_level_msg_line_file("INFO", "1th message is logged", expected_line_number);

  axom::slic::setIsRoot(true);
  for(int i = 0; i < 3; i++)
  {
    SLIC_INFO_ROOT_IF_ONCE(i > 0, i << "th message is logged!");
  }
  expected_line_number = __LINE__ - 2;
  EXPECT_EQ(check_count(slic::internal::test_stream.str(), "INFO"), 1);
  check_level_msg_line_file("INFO", "1th message is logged", expected_line_number);

  // Two call-sites have a single message each
  for(int i = 0; i < 3; i++)
  {
    SLIC_INFO_IF_ONCE(i == 0, "message 1 logs " << i);
    SLIC_INFO_IF_ONCE(i > 0, "message 2 logs " << i);
  }
  int msg_1_line = __LINE__ - 2;
  int msg_2_line = __LINE__ - 3;

  EXPECT_FALSE(slic::internal::is_stream_empty());
  const std::string str = slic::internal::test_stream.str();
  EXPECT_EQ(check_count(str, "INFO"), 2);
  check_level(str, "INFO");
  check_msg(str, "message 1 logs 0");
  check_msg(str.substr(str.size() / 2), "message 2 logs 1");
  check_line(str, msg_1_line);
  check_line(str, msg_2_line);
  check_file(str);
  slic::internal::clear();
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_macros_file_output)
{
  int expected_line_number;

  std::string msgfmt = "<MESSAGE>";

  // GenericOutputStream(std::string stream) and
  // GenericOutputStream(std::string stream, std::string format) constructors
  // do not create a a file until macros called, then flushed
  std::string no_fmt = "file_no_fmt.txt";
  std::string with_fmt = "file_with_fmt.txt";

  slic::addStreamToAllMsgLevels(new slic::GenericOutputStream(no_fmt));
  slic::addStreamToAllMsgLevels(new slic::GenericOutputStream(with_fmt, msgfmt));

  EXPECT_FALSE(axom::utilities::filesystem::pathExists(no_fmt));
  EXPECT_FALSE(axom::utilities::filesystem::pathExists(with_fmt));

  // streams flushed with no buffered messages, no files created
  slic::flushStreams();

  EXPECT_FALSE(axom::utilities::filesystem::pathExists(no_fmt));
  EXPECT_FALSE(axom::utilities::filesystem::pathExists(with_fmt));

  // message is buffered but not yet flushed, no files created
  SLIC_INFO("Test");
  expected_line_number = __LINE__ - 1;

  EXPECT_FALSE(axom::utilities::filesystem::pathExists(no_fmt));
  EXPECT_FALSE(axom::utilities::filesystem::pathExists(with_fmt));

  // message has been buffered and now flushed, files are created
  slic::flushStreams();

  EXPECT_TRUE(axom::utilities::filesystem::pathExists(no_fmt));
  EXPECT_TRUE(axom::utilities::filesystem::pathExists(with_fmt));

  // Verify file contents
  std::ifstream no_fmt_contents(no_fmt);
  std::stringstream no_fmt_buffer;
  no_fmt_buffer << no_fmt_contents.rdbuf();

  std::string no_fmt_expected;
  no_fmt_expected += "*****\n[INFO]\n\n Test \n\n ";
  no_fmt_expected += __FILE__;
  no_fmt_expected += "\n" + std::to_string(expected_line_number);
  no_fmt_expected += "\n****\n";

  EXPECT_EQ(no_fmt_buffer.str(), no_fmt_expected);

  std::ifstream with_fmt_contents(with_fmt);
  std::stringstream with_fmt_buffer;
  with_fmt_buffer << with_fmt_contents.rdbuf();

  EXPECT_EQ(with_fmt_buffer.str(), "Test");

  no_fmt_contents.close();
  with_fmt_contents.close();

  // Closes open file streams associated with Slic streams when destructors
  // called during slic::finalize().
  // Windows _unlink file deletion fails if file is still in use.
  slic::finalize();

  // Cleanup generated files
  EXPECT_EQ(axom::utilities::filesystem::removeFile(no_fmt), 0);
  EXPECT_EQ(axom::utilities::filesystem::removeFile(with_fmt), 0);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // initialize slic
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError(); /* disable abort for testing purposes */

  std::string msgfmt = "[<LEVEL>]:;;<MESSAGE>;;\n@@<FILE>\n@@<LINE>";

  slic::addStreamToAllMsgLevels(new slic::GenericOutputStream(&slic::internal::test_stream, msgfmt));

  std::string msgtagfmt = "[<LEVEL>]:;;<MESSAGE>;;\n##<TAG>\n@@<FILE>\n@@<LINE>";
  slic::addStreamToTag(new slic::GenericOutputStream(&slic::internal::test_stream, msgtagfmt),
                       "myTag");

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();

  // Finalize slic
  slic::finalize();

  return result;
}
