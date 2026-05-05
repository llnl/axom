// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/FileUtilities.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/slic/streams/SynchronizedStream.hpp"
#include "axom/slic/streams/LumberjackStream.hpp"

#include "gtest/gtest.h"

#include <string>
#include <sstream>

#include "mpi.h"

// namespace alias
namespace slic = axom::slic;

namespace axom
{
namespace slic
{
namespace internal
{
// Stream for message levels
std::ostringstream test_stream;

// Stream for tagged streams
std::ostringstream test_tag_stream;

bool is_stream_empty() { return test_stream.str().empty(); }

bool is_tag_stream_empty() { return test_tag_stream.str().empty(); }

bool are_all_streams_empty() { return (is_stream_empty() && is_tag_stream_empty()); }

void clear_streams()
{
  test_stream.clear();
  test_stream.str("");

  test_tag_stream.clear();
  test_tag_stream.str("");

  EXPECT_TRUE(are_all_streams_empty());
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

  // message is output only once (combining successful for LumberjackStream)
  int count = 0;
  size_t pos = 0;

  while((pos = msg.find(expected_message, pos)) != std::string::npos)
  {
    count++;
    pos += expected_message.length();
  }

  EXPECT_EQ(count, 1);

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
// For SynchronizedStream where each message has one associated rank
void check_rank(const std::string& msg, int expected_rank)
{
  EXPECT_FALSE(msg.empty());

  // extract rank
  size_t start = msg.rfind("$$") + 2;
  size_t end = msg.rfind("&&");
  std::string ranks_string = msg.substr(start, end - start);

  EXPECT_TRUE(ranks_string.find(std::to_string(expected_rank)) != std::string::npos);
}

//------------------------------------------------------------------------------
// For LumberjackStream where message from multiple ranks are at the root
void check_ranks(const std::string& msg, int expected_ranks)
{
  // Check all ranks from [0, expected_ranks) are in message
  for(int i = 0; i < expected_ranks; i++)
  {
    check_rank(msg, i);
  }
}

//------------------------------------------------------------------------------
void check_rank_count(const std::string& msg, const std::string& stream_type, int expected_rank_count)
{
  EXPECT_FALSE(msg.empty());

  // Always 1 for SynchronizedStream
  if(stream_type == "Synchronized")
  {
    expected_rank_count = 1;
  }

  // extract rank
  size_t start = msg.rfind("&&") + 2;
  std::string rc = msg.substr(start);
  const int rank_count = std::stoi(rc);
  EXPECT_EQ(rank_count, expected_rank_count);
}

//------------------------------------------------------------------------------
// Use level to determine number of messages - used for SLIC_*_ONCE macros
int check_msg_count(const std::string& msg, const std::string& expected_level)
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
// Checks message logged on all ranks
void check_all_ranks(const std::string& stream_type,
                     const std::string& level,
                     const std::string& message,
                     int expected_line,
                     int rank,
                     int nranks)
{
  if(stream_type == "Synchronized" || (stream_type == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    const std::string str = slic::internal::test_stream.str();

    SCOPED_TRACE(std::string("SLIC trace (all ranks): Failed line was ") +
                 std::to_string(expected_line));

    check_level(str, level);
    check_msg(str, message);
    check_line(str, expected_line);
    if(stream_type == "Synchronized")
    {
      check_rank(str, rank);
    }
    else
    {
      check_ranks(str, nranks);
    }
    check_rank_count(str, stream_type, nranks);

    check_file(str);
  }
}

// Convenience test macro that calls given slic macro and checks
// message logged on all ranks
#define EXPECT_SLIC_LOG_ALL_RANKS(macro_call, level, message)                 \
  do                                                                          \
  {                                                                           \
    const int expected_line = __LINE__;                                       \
    macro_call;                                                               \
    slic::flushStreams();                                                     \
    check_all_ranks(GetParam(), level, message, expected_line, rank, nranks); \
    slic::internal::clear_streams();                                          \
  } while(false)

//------------------------------------------------------------------------------
// Checks message logged only on root
void check_root(const std::string& stream_type,
                const std::string& level,
                const std::string& message,
                int expected_line,
                int rank)
{
  if(rank == 0)
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    const std::string str = slic::internal::test_stream.str();

    SCOPED_TRACE(std::string("SLIC trace (root): Failed line was ") + std::to_string(expected_line));

    check_level(str, level);
    check_msg(str, message);
    check_line(str, expected_line);
    check_rank(str, rank);

    // Only one rank has logged a message
    check_rank_count(str, stream_type, 1);

    check_file(str);
  }
  else
  {
    EXPECT_TRUE(slic::internal::are_all_streams_empty());
  }
}

// Convenience test macro that calls given slic macro and checks
// message logged only on root
#define EXPECT_SLIC_LOG_ROOT(macro_call, level, message)         \
  do                                                             \
  {                                                              \
    axom::slic::setIsRoot(rank == 0);                            \
    const int expected_line = __LINE__;                          \
    macro_call;                                                  \
    slic::flushStreams();                                        \
    check_root(GetParam(), level, message, expected_line, rank); \
    slic::internal::clear_streams();                             \
  } while(false)

//------------------------------------------------------------------------------
// Checks message logged only on even ranks
void check_even(const std::string& stream_type,
                const std::string& level,
                const std::string& message,
                int expected_line,
                int rank,
                int nranks)
{
  if(((rank % 2) == 0 && stream_type == "Synchronized") || (rank == 0 && stream_type == "Lumberjack"))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    const std::string str = slic::internal::test_stream.str();

    SCOPED_TRACE(std::string("SLIC trace (even): Failed line was ") + std::to_string(expected_line));

    check_level(str, level);
    check_msg(str, message);
    check_line(str, expected_line);
    if(stream_type == "Synchronized")
    {
      check_rank(str, rank);
    }
    else
    {
      for(int i = 0; i < nranks; i += 2)
      {
        check_rank(str, i);
      }
    }
    check_rank_count(str, stream_type, (nranks / 2) + (nranks % 2));
    check_file(str);
  }
  else
  {
    EXPECT_TRUE(slic::internal::are_all_streams_empty());
  }
}

// Convenience test macro that calls given slic macro and checks
// message logged only on even ranks
#define EXPECT_SLIC_LOG_EVEN(macro_call, level, message)                 \
  do                                                                     \
  {                                                                      \
    const int expected_line = __LINE__;                                  \
    macro_call;                                                          \
    slic::flushStreams();                                                \
    check_even(GetParam(), level, message, expected_line, rank, nranks); \
    slic::internal::clear_streams();                                     \
  } while(false)

//------------------------------------------------------------------------------
bool has_aborted = false;
void custom_abort_function() { has_aborted = true; }

void reset_state()
{
  axom::slic::internal::clear_streams();
  has_aborted = false;
  AXOM_UNUSED_VAR(has_aborted);
}
}  // end anonymous namespace

//------------------------------------------------------------------------------
// TEST FIXTURE CLASS
//------------------------------------------------------------------------------
class SlicMacrosParallel : public ::testing::TestWithParam<std::string>
{
public:
  void SetUp() override
  {
    stream_type = GetParam();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    // initialize slic
    slic::initialize();
    slic::setLoggingMsgLevel(slic::message::Debug);
    slic::disableAbortOnError(); /* disable abort for testing purposes */

    if(stream_type == "Lumberjack")
    {
      slic::addStreamToAllMsgLevels(
        new slic::LumberjackStream(&slic::internal::test_stream, MPI_COMM_WORLD, RLIMIT, msgfmt));

      slic::addStreamToTag(
        new slic::LumberjackStream(&slic::internal::test_tag_stream, MPI_COMM_WORLD, RLIMIT, msgtagfmt),
        "myTag");
    }

    if(stream_type == "Synchronized")
    {
      slic::addStreamToAllMsgLevels(
        new slic::SynchronizedStream(&slic::internal::test_stream, MPI_COMM_WORLD, msgfmt));

      slic::addStreamToTag(
        new slic::SynchronizedStream(&slic::internal::test_tag_stream, MPI_COMM_WORLD, msgtagfmt),
        "myTag");
    }
  }

  void TearDown() override { slic::finalize(); }

  std::string stream_type;

  int rank;
  int nranks;
  const int RLIMIT = 8;

  const std::string msgfmt =
    "[<LEVEL>]:;;<MESSAGE>;;\n@@<FILE>\n@@<LINE>\n$$<RANK>\n&&<RANK_COUNT>";

  const std::string msgtagfmt =
    "[<LEVEL>]:;;<MESSAGE>;;\n##<TAG>\n@@<FILE>\n@@<LINE>\n$$<RANK>\n&&<RANK_"
    "COUNT>";
};

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_error_macros)
{
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_ERROR("test error message"), "ERROR", "test error message");

  SLIC_ERROR_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_ERROR_IF(true, "this message is logged!"), "ERROR", "this message is logged!");
  // clang-format on

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_ERROR_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_ERROR_ROOT_IF(true, "this message is logged!"), "ERROR", "this message is logged!");
  // clang-format on

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_ERROR_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_ERROR_ROOT_IF(true, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // clang-format off
  // Check for message only on root
  EXPECT_SLIC_LOG_ROOT(SLIC_ERROR_ROOT_IF(true, "this message is logged!"), "ERROR", "this message is logged!");

  // Check for message on every even rank only
  EXPECT_SLIC_LOG_EVEN(SLIC_ERROR_IF((rank % 2) == 0, "this message is logged!"), "ERROR", "this message is logged!");
  // clang-format on
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_warning_macros)
{
  int expected_line_number;

  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_WARNING("test warning message"), "WARNING", "test warning message");

  SLIC_WARNING_IF(false, "this message should not be logged!");
  SLIC_WARNING_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Called once per call site; since this test is ran for parameters
  // Synchronized and Lumberjack, a separate call site is needed for each parameter,
  // otherwise ONCE macro logs nothing when its Lumberjack's turn to run.
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_WARNING_ONCE("test warning message " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_WARNING_ONCE("test warning message " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "WARNING"), 1);
  }
  check_all_ranks(GetParam(), "WARNING", "test warning message 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_WARNING_IF(true, "this message is logged!"), "WARNING", "this message is logged!");
  // clang-format on

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_WARNING_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_WARNING_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "WARNING"), 1);
  }
  check_all_ranks(GetParam(), "WARNING", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_WARNING_ROOT_IF(false, "this message should not be logged!");
  SLIC_WARNING_ROOT_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_WARNING_ROOT_IF(true, "this message is logged!"), "WARNING", "this message is logged!");
  // clang-format on

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_WARNING_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_WARNING_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "WARNING"), 1);
  }
  check_all_ranks(GetParam(), "WARNING", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_WARNING_ROOT_IF(false, "this message should not be logged!");
  SLIC_WARNING_ROOT_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_WARNING_ROOT_IF(true, "this message should not be logged!");
  SLIC_WARNING_ROOT_IF_ONCE(true, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // clang-format off
  // Check for message only on root
  EXPECT_SLIC_LOG_ROOT(SLIC_WARNING_ROOT_IF(true, "this message is logged!"), "WARNING", "this message is logged!");

  axom::slic::setIsRoot(rank == 0);
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_WARNING_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_WARNING_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(rank == 0)
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "WARNING"), 1);
  }
  check_root(GetParam(), "WARNING", "this message is logged 0", expected_line_number, rank);
  slic::internal::clear_streams();

  // Check for message on every even rank only
  EXPECT_SLIC_LOG_EVEN(SLIC_WARNING_IF((rank % 2) == 0, "this message is logged!"), "WARNING", "this message is logged!");
  // clang-format on

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_WARNING_IF_ONCE((rank % 2) == 0, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_WARNING_IF_ONCE((rank % 2) == 0, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(((rank % 2) == 0 && stream_type == "Synchronized") || (rank == 0 && stream_type == "Lumberjack"))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "WARNING"), 1);
  }
  check_even(GetParam(), "WARNING", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_info_macros)
{
  int expected_line_number;
  int expected_tag_number;

  EXPECT_TRUE(slic::internal::are_all_streams_empty());
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_INFO("test info message"), "INFO", "test info message");

  // Called once per call site; since this test is ran for parameters
  // Synchronized and Lumberjack, a separate call site is needed for each parameter,
  // otherwise ONCE macro logs nothing when its Lumberjack's turn to run.
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_INFO_ONCE("test info message " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_INFO_ONCE("test info message " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "INFO"), 1);
  }
  check_all_ranks(GetParam(), "INFO", "test info message 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  SLIC_INFO_TAGGED("test tagged info message", "myTag");
  expected_tag_number = __LINE__ - 1;
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_tag_stream_empty());
    check_level(slic::internal::test_tag_stream.str(), "INFO");
    check_msg(slic::internal::test_tag_stream.str(), "test tagged info message");
    check_file(slic::internal::test_tag_stream.str());
    check_line(slic::internal::test_tag_stream.str(), expected_tag_number);
    check_tag(slic::internal::test_tag_stream.str(), "myTag");
    if(GetParam() == "Synchronized")
    {
      check_rank(slic::internal::test_tag_stream.str(), rank);
    }
    else
    {
      check_ranks(slic::internal::test_tag_stream.str(), nranks);
    }
    check_rank_count(slic::internal::test_tag_stream.str(), GetParam(), nranks);
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();

  SLIC_INFO("test info message only for normal message-level stream");
  expected_line_number = __LINE__ - 1;
  SLIC_INFO_TAGGED("test tagged info message only for tagged stream", "myTag");
  expected_tag_number = __LINE__ - 1;
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "INFO");
    check_msg(slic::internal::test_stream.str(),
              "test info message only for normal message-level stream");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), expected_line_number);
    if(GetParam() == "Synchronized")
    {
      check_rank(slic::internal::test_stream.str(), rank);
    }
    else
    {
      check_ranks(slic::internal::test_stream.str(), nranks);
    }
    check_rank_count(slic::internal::test_stream.str(), GetParam(), nranks);

    EXPECT_FALSE(slic::internal::is_tag_stream_empty());
    check_level(slic::internal::test_tag_stream.str(), "INFO");
    check_msg(slic::internal::test_tag_stream.str(),
              "test tagged info message only for tagged stream");
    check_file(slic::internal::test_tag_stream.str());
    check_line(slic::internal::test_tag_stream.str(), expected_tag_number);
    check_tag(slic::internal::test_tag_stream.str(), "myTag");
    if(GetParam() == "Synchronized")
    {
      check_rank(slic::internal::test_tag_stream.str(), rank);
    }
    else
    {
      check_ranks(slic::internal::test_tag_stream.str(), nranks);
    }
    check_rank_count(slic::internal::test_tag_stream.str(), GetParam(), nranks);
  }
  slic::internal::clear_streams();

  SLIC_INFO_TAGGED("this message should not be logged (no tag given)!", "");
  SLIC_INFO_TAGGED_ONCE("this message should not be logged (no tag given)!", "");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  SLIC_INFO_TAGGED("this message should not be logged (tag DNE)!", "tag404");
  SLIC_INFO_TAGGED_ONCE("this message should not be logged (tag DNE)!", "tag404");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  SLIC_INFO_IF(false, "this message should not be logged!");
  SLIC_INFO_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_INFO_IF(true, "this message is logged!"), "INFO", "this message is logged!");
  // clang-format on

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_INFO_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_INFO_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "INFO"), 1);
  }
  check_all_ranks(GetParam(), "INFO", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_INFO_ROOT_IF(false, "this message should not be logged!");
  SLIC_INFO_ROOT_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_INFO_ROOT_IF(true, "this message is logged!"), "INFO", "this message is logged!");
  // clang-format on

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_INFO_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_INFO_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "INFO"), 1);
  }
  check_all_ranks(GetParam(), "INFO", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_INFO_ROOT_IF(false, "this message should not be logged!");
  SLIC_INFO_ROOT_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_INFO_ROOT_IF(true, "this message should not be logged!");
  SLIC_INFO_ROOT_IF_ONCE(true, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // clang-format off
  // Check for message only on root
  EXPECT_SLIC_LOG_ROOT(SLIC_INFO_ROOT_IF(true, "this message is logged!"), "INFO", "this message is logged!");

  axom::slic::setIsRoot(rank == 0);
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_INFO_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_INFO_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(rank == 0)
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "INFO"), 1);
  }
  check_root(GetParam(), "INFO", "this message is logged 0", expected_line_number, rank);
  slic::internal::clear_streams();

  // Check for message on every even rank only
  EXPECT_SLIC_LOG_EVEN(SLIC_INFO_IF((rank % 2) == 0, "this message is logged!"), "INFO", "this message is logged!");
  // clang-format on

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_INFO_IF_ONCE((rank % 2) == 0, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_INFO_IF_ONCE((rank % 2) == 0, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(((rank % 2) == 0 && stream_type == "Synchronized") || (rank == 0 && stream_type == "Lumberjack"))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "INFO"), 1);
  }
  check_even(GetParam(), "INFO", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_debug_macros)
{
  int expected_line_number;

  EXPECT_TRUE(slic::internal::are_all_streams_empty());

#ifdef AXOM_DEBUG

  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_DEBUG("test debug message"), "DEBUG", "test debug message");

  // Called once per call site; since this test is ran for parameters
  // Synchronized and Lumberjack, a separate call site is needed for each parameter,
  // otherwise ONCE macro logs nothing when its Lumberjack's turn to run.
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_DEBUG_ONCE("test debug message " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_DEBUG_ONCE("test debug message " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "DEBUG"), 1);
  }
  check_all_ranks(GetParam(), "DEBUG", "test debug message 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_DEBUG_IF(true, "this message is logged!"), "DEBUG", "this message is logged!");

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_DEBUG_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_DEBUG_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "DEBUG"), 1);
  }
  check_all_ranks(GetParam(), "DEBUG", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_DEBUG_ROOT_IF(true, "this message is logged!"), "DEBUG", "this message is logged!");

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_DEBUG_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_DEBUG_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "DEBUG"), 1);
  }
  check_all_ranks(GetParam(), "DEBUG", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // Check for message only on root
  EXPECT_SLIC_LOG_ROOT(SLIC_DEBUG_ROOT_IF(true, "this message is logged!"), "DEBUG", "this message is logged!");

  axom::slic::setIsRoot(rank == 0);
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_DEBUG_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_DEBUG_ROOT_IF_ONCE(true, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(rank == 0)
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "DEBUG"), 1);
  }
  check_root(GetParam(), "DEBUG", "this message is logged 0", expected_line_number, rank);
  slic::internal::clear_streams();

  // Check for message on every even rank only
  EXPECT_SLIC_LOG_EVEN(SLIC_DEBUG_IF((rank % 2) == 0, "this message is logged!"), "DEBUG", "this message is logged!");
  // clang-format on

  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_DEBUG_IF_ONCE((rank % 2) == 0, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_DEBUG_IF_ONCE((rank % 2) == 0, "this message is logged " << i);
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(((rank % 2) == 0 && stream_type == "Synchronized") || (rank == 0 && stream_type == "Lumberjack"))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "DEBUG"), 1);
  }
  check_even(GetParam(), "DEBUG", "this message is logged 0", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(expected_line_number);

  SLIC_DEBUG("test debug message");
  SLIC_DEBUG_ONCE("test debug message");

  SLIC_DEBUG_IF(true, "this message is logged!");
  SLIC_DEBUG_IF_ONCE(true, "this message is logged!");

  axom::slic::setIsRoot(true);
  SLIC_DEBUG_ROOT_IF(true, "this message is logged!");
  SLIC_DEBUG_ROOT_IF_ONCE(true, "this message is logged!");

  SLIC_DEBUG_IF((rank % 2) == 0, "this message is logged!");
  SLIC_DEBUG_IF_ONCE((rank % 2) == 0, "this message is logged!");

  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());
#endif

  SLIC_DEBUG_IF(false, "this message should not be logged!");
  SLIC_DEBUG_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_DEBUG_ROOT_IF(false, "this message should not be logged!");
  SLIC_DEBUG_ROOT_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_DEBUG_ROOT_IF(false, "this message should not be logged!");
  SLIC_DEBUG_ROOT_IF_ONCE(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_DEBUG_ROOT_IF(true, "this message should not be logged!");
  SLIC_DEBUG_ROOT_IF_ONCE(true, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());
}

TEST_P(SlicMacrosParallel, test_abort_error_macros)
{
  int expected_line_number;

  const int NUM_ABORT_STATES = 2;

  slic::enableAbortOnError(); /* enable abort for testing purposes */
  slic::setAbortFunction(custom_abort_function);

  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Test for each rank
  for(int i = 0; i < nranks; i++)
  {
    // Test abort enabled/disabled
    for(int j = 0; j < NUM_ABORT_STATES; j++)
    {
      if(j == 0)
      {
        slic::enableAbortOnError();
      }
      else
      {
        slic::disableAbortOnError();
      }

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

      bool abort_enabled = slic::isAbortOnErrorsEnabled();

      axom::slic::setIsRoot(rank == i);

      if(rank == i)
      {
        SLIC_ERROR("SLIC_ERROR message is logged!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(), "SLIC_ERROR message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        int val = rank == i ? 42 : -42;
        SLIC_ERROR_IF(val == 42, "SLIC_ERROR_IF message is logged!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(), "SLIC_ERROR_IF message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        SLIC_ERROR_ROOT("SLIC_ERROR_ROOT message is logged!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(), "SLIC_ERROR_ROOT message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        SLIC_ERROR_ROOT_IF(val == 42, "SLIC_ERROR_ROOT_IF message is logged!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(), "SLIC_ERROR_ROOT_IF message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        SLIC_ASSERT(val < 0);
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(), "Failed Assert: val < 0");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        SLIC_ASSERT_MSG(val < 0, "val should be negative!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(),
                  "Failed Assert: val < 0\nval should be negative!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();
      }
      axom::slic::setIsRoot(true);
#else
      // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined

      // Quiet warning about has_aborted and reset_state never being referenced
      AXOM_UNUSED_VAR(has_aborted);
      AXOM_UNUSED_VAR(expected_line_number);
      reset_state();

      EXPECT_TRUE(slic::internal::are_all_streams_empty());
#endif

    }  //end NUM_ABORT_STATES loop
  }  // end nranks loop

  slic::disableAbortOnError(); /* disable abort for testing purposes */
  slic::setAbortFunction(axom::utilities::processAbort);
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_abort_warning_macros)
{
  int expected_line_number;

  const int NUM_ABORT_STATES = 2;

  slic::enableAbortOnWarning(); /* enable abort for testing purposes */
  slic::setAbortFunction(custom_abort_function);

  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  // Test for each rank
  for(int i = 0; i < nranks; i++)
  {
    // Test abort enabled/disabled
    for(int j = 0; j < NUM_ABORT_STATES; j++)
    {
      if(j == 0)
      {
        slic::enableAbortOnWarning();
      }
      else
      {
        slic::disableAbortOnWarning();
      }

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

      bool abort_enabled = slic::isAbortOnWarningsEnabled();

      axom::slic::setIsRoot(rank == i);

      if(rank == i)
      {
        SLIC_WARNING("SLIC_WARNING message is logged!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(), "SLIC_WARNING message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        int val = rank == i ? 42 : -42;
        SLIC_WARNING_IF(val == 42, "SLIC_WARNING_IF message is logged!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(), "SLIC_WARNING_IF message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        SLIC_WARNING_ROOT("SLIC_WARNING_ROOT message is logged!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(), "SLIC_WARNING_ROOT message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        SLIC_WARNING_ROOT_IF(val == 42, "SLIC_WARNING_ROOT_IF msg logged!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(), "SLIC_WARNING_ROOT_IF msg logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        SLIC_CHECK(val < 0);
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(), "Failed Check: val < 0");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();

        SLIC_CHECK_MSG(val < 0, "val should be negative!");
        expected_line_number = __LINE__ - 1;
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::are_all_streams_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(),
                  "Failed Check: val < 0\nval should be negative!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), expected_line_number);
        check_rank(slic::internal::test_stream.str(), rank);
        check_rank_count(slic::internal::test_stream.str(), GetParam(), 1);
        reset_state();
      }

      axom::slic::setIsRoot(true);

#else
      AXOM_UNUSED_VAR(expected_line_number);
      // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
      EXPECT_TRUE(slic::internal::are_all_streams_empty());
#endif

    }  //end NUM_ABORT_STATES loop
  }  // end nranks loop

  slic::disableAbortOnWarning(); /* disable abort for testing purposes */
  slic::setAbortFunction(axom::utilities::processAbort);
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_assert_macros)
{
  slic::internal::clear_streams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());
  constexpr int val = 42;

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_ASSERT(val < 0), "ERROR", "Failed Assert: val < 0");

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_ASSERT_MSG(val < 0, "val should be negative!"), "ERROR", "Failed Assert: val < 0\nval should be negative!");
  // clang-format on

#else
  // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);

  SLIC_ASSERT(val < 0);

  SLIC_ASSERT_MSG(val < 0, "val should be negative!");

  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());
#endif

  SLIC_ASSERT(val > 0);
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());
}

// ------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_check_macros)
{
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  constexpr int val = 42;

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_CHECK(val < 0), "WARNING", "Failed Check: val < 0");

  // Single line - Placement of ")" matters for __LINE__ for slic call and checking
  // clang-format off
  EXPECT_SLIC_LOG_ALL_RANKS(SLIC_CHECK_MSG(val < 0, "val should be negative!"), "WARNING", "Failed Check: val < 0\nval should be negative!");
  // clang-format on

#else
  // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);

  SLIC_CHECK(val < 0);

  SLIC_CHECK_MSG(val < 0, "val should be negative!");

  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());
#endif

  SLIC_CHECK(val > 0);
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::are_all_streams_empty());
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_if_once_macros)
{
  int expected_line_number;

  // Check that message is logged when condition is satisfied only once
  // Called once per call site; since this test is ran for parameters
  // Synchronized and Lumberjack, a separate call site is needed for each parameter,
  // otherwise ONCE macro logs nothing when its Lumberjack's turn to run.
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_INFO_IF_ONCE(i > 0, i << "th message is logged!");
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_INFO_IF_ONCE(i > 0, i << "th message is logged!");
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "INFO"), 1);
  }
  check_all_ranks(GetParam(), "INFO", "1th message is logged", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  axom::slic::setIsRoot(true);
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_INFO_ROOT_IF_ONCE(i > 0, i << "th message is logged!");
      expected_line_number = __LINE__ - 1;
    }
    else
    {
      SLIC_INFO_ROOT_IF_ONCE(i > 0, i << "th message is logged!");
      expected_line_number = __LINE__ - 1;
    }
  }
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    EXPECT_EQ(check_msg_count(slic::internal::test_stream.str(), "INFO"), 1);
  }
  check_all_ranks(GetParam(), "INFO", "1th message is logged", expected_line_number, rank, nranks);
  slic::internal::clear_streams();

  // Two call-sites have a single message each
  int msg_1_line;
  int msg_2_line;
  for(int i = 0; i < 3; i++)
  {
    if(GetParam() == "Synchronized")
    {
      SLIC_INFO_IF_ONCE(i == 0, "message 1 logs " << i);
      msg_1_line = __LINE__ - 1;
      SLIC_INFO_IF_ONCE(i > 0, "message 2 logs " << i);
      msg_2_line = __LINE__ - 1;
    }
    else
    {
      SLIC_INFO_IF_ONCE(i == 0, "message 1 logs " << i);
      msg_1_line = __LINE__ - 1;
      SLIC_INFO_IF_ONCE(i > 0, "message 2 logs " << i);
      msg_2_line = __LINE__ - 1;
    }
  }
  slic::flushStreams();

  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::are_all_streams_empty());
    const std::string str = slic::internal::test_stream.str();
    EXPECT_EQ(check_msg_count(str, "INFO"), 2);

    check_level(str, "INFO");

    std::string msg_1 = str.substr(0, str.size() / 2);
    std::string msg_2 = str.substr(str.size() / 2);

    check_level(msg_1, "INFO");
    check_level(msg_2, "INFO");
    check_msg(msg_1, "message 1 logs 0");
    check_msg(msg_2, "message 2 logs 1");
    check_line(msg_1, msg_1_line);
    check_line(msg_2, msg_2_line);
    check_file(msg_1);
    check_file(msg_2);
    if(GetParam() == "Synchronized")
    {
      check_rank(msg_1, rank);
      check_rank(msg_2, rank);
    }
    else
    {
      check_ranks(msg_1, nranks);
      check_ranks(msg_2, nranks);
    }
    check_rank_count(msg_1, GetParam(), nranks);
    check_rank_count(msg_2, GetParam(), nranks);
    check_file(msg_1);
    check_file(msg_2);
  }
  slic::internal::clear_streams();
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_macros_file_output)
{
  EXPECT_TRUE(slic::internal::are_all_streams_empty());

  std::string no_fmt;
  std::string with_fmt;

  if(GetParam() == "Synchronized")
  {
    // SynchronizedStream(std::string stream, MPI_Comm comm) and
    // SynchronizedStream(std::string stream, MPI_Comm comm, std::string format)
    // constructors do not create a a file until macros called, then flushed

    no_fmt = "file_ss_rank_" + std::to_string(rank) + "_no_fmt.txt";
    with_fmt = "file_ss_rank_" + std::to_string(rank) + "_with_fmt.txt";

    slic::addStreamToAllMsgLevels(new slic::SynchronizedStream(no_fmt, MPI_COMM_WORLD));

    slic::addStreamToAllMsgLevels(new slic::SynchronizedStream(with_fmt, MPI_COMM_WORLD, msgfmt));
  }

  else
  {
    // LumberjackStream(std::string stream, MPI_Comm comm, int ranksLimit) and
    // LumberjackStream(std::string stream, MPI_Comm comm, int ranksLimit,
    //                  std::string format)
    // constructors do not create a a file until macros called, then flushed

    no_fmt = "file_lj_rank_" + std::to_string(rank) + "_no_fmt.txt";
    with_fmt = "file_lj_rank_" + std::to_string(rank) + "_with_fmt.txt";

    slic::addStreamToAllMsgLevels(new slic::LumberjackStream(no_fmt, MPI_COMM_WORLD, RLIMIT));

    slic::addStreamToAllMsgLevels(
      new slic::LumberjackStream(with_fmt, MPI_COMM_WORLD, RLIMIT, msgfmt));
  }

  EXPECT_FALSE(axom::utilities::filesystem::pathExists(no_fmt));
  EXPECT_FALSE(axom::utilities::filesystem::pathExists(with_fmt));

  // streams flushed with no buffered messages, no files created
  slic::flushStreams();

  EXPECT_FALSE(axom::utilities::filesystem::pathExists(no_fmt));
  EXPECT_FALSE(axom::utilities::filesystem::pathExists(with_fmt));

  // message is buffered but not yet flushed, no files created
  SLIC_INFO("Test");
  int expected_line_number_flush = __LINE__ - 1;

  EXPECT_FALSE(axom::utilities::filesystem::pathExists(no_fmt));
  EXPECT_FALSE(axom::utilities::filesystem::pathExists(with_fmt));

  // message has been buffered and now flushed, files are created
  slic::flushStreams();

  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_TRUE(axom::utilities::filesystem::pathExists(no_fmt));
    EXPECT_TRUE(axom::utilities::filesystem::pathExists(with_fmt));

    // Verify file contents
    std::ifstream no_fmt_contents(no_fmt);
    std::stringstream no_fmt_buffer;
    no_fmt_buffer << no_fmt_contents.rdbuf();
    no_fmt_contents.close();

    std::string no_fmt_expected;
    no_fmt_expected += "*****\n[INFO]\n\n Test \n\n ";
    no_fmt_expected += __FILE__;
    no_fmt_expected += "\n" + std::to_string(expected_line_number_flush);
    no_fmt_expected += "\n****\n";

    EXPECT_EQ(no_fmt_buffer.str(), no_fmt_expected);

    std::ifstream with_fmt_contents(with_fmt);
    std::stringstream with_fmt_buffer;
    with_fmt_buffer << with_fmt_contents.rdbuf();
    with_fmt_contents.close();

    check_level(with_fmt_buffer.str(), "INFO");
    check_msg(with_fmt_buffer.str(), "Test");
    check_file(with_fmt_buffer.str());
    check_line(with_fmt_buffer.str(), expected_line_number_flush);
    if(GetParam() == "Synchronized")
    {
      check_rank(with_fmt_buffer.str(), rank);
    }
    else
    {
      check_ranks(with_fmt_buffer.str(), nranks);
    }
    check_rank_count(with_fmt_buffer.str(), GetParam(), nranks);
  }

  else
  {
    // Expect non-output Lumberjack ranks to not create files
    EXPECT_FALSE(axom::utilities::filesystem::pathExists(no_fmt));
    EXPECT_FALSE(axom::utilities::filesystem::pathExists(with_fmt));
  }

  // In case of an abort, outputLocal() is called.
  // Expect non-output Lumberjack ranks to create files if possible
  // (cannot guarantee all ranks will output before non-collective MPI Abort)
  SLIC_INFO("Test outputLocalMessages()");
  int expected_line_number_local = __LINE__ - 1;
  slic::outputLocalMessages();

  EXPECT_TRUE(axom::utilities::filesystem::pathExists(no_fmt));
  EXPECT_TRUE(axom::utilities::filesystem::pathExists(with_fmt));

  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    // Verify file contents
    std::ifstream no_fmt_output(no_fmt);
    std::stringstream no_fmt_out_buf;
    no_fmt_out_buf << no_fmt_output.rdbuf();
    no_fmt_output.close();

    std::string no_fmt_output_expected;
    no_fmt_output_expected += "*****\n[INFO]\n\n Test \n\n ";
    no_fmt_output_expected += __FILE__;
    no_fmt_output_expected += "\n" + std::to_string(expected_line_number_flush);
    no_fmt_output_expected += "\n****\n";
    no_fmt_output_expected += "*****\n[INFO]\n\n Test outputLocalMessages() \n\n ";
    no_fmt_output_expected += __FILE__;
    no_fmt_output_expected += "\n" + std::to_string(expected_line_number_local);
    no_fmt_output_expected += "\n****\n";

    EXPECT_EQ(no_fmt_out_buf.str(), no_fmt_output_expected);

    std::ifstream with_fmt_output(with_fmt);
    std::stringstream with_fmt_out_buf;
    with_fmt_out_buf << with_fmt_output.rdbuf();
    with_fmt_output.close();

    // Check the next message, "Test outputLocalMessages()"
    size_t next_index = with_fmt_out_buf.str().rfind("[INFO]");
    std::string next_msg = with_fmt_out_buf.str().substr(next_index);
    check_level(next_msg, "INFO");
    check_msg(next_msg, "Test outputLocalMessages()");
    check_file(next_msg);
    check_line(next_msg, (expected_line_number_local));

    // For outputLocalMessages(), only current rank and rank count of 1
    // output for LumberjackStreams. Behaves like SynchronizedStream.
    check_rank(next_msg, rank);
    check_rank_count(next_msg, GetParam(), 1);
  }

  else
  {
    // Verify file contents
    std::ifstream no_fmt_output(no_fmt);
    std::stringstream no_fmt_out_buf;
    no_fmt_out_buf << no_fmt_output.rdbuf();
    no_fmt_output.close();

    std::string no_fmt_output_expected;
    no_fmt_output_expected += "*****\n[INFO]\n\n Test outputLocalMessages() \n\n ";
    no_fmt_output_expected += __FILE__;
    no_fmt_output_expected += "\n" + std::to_string(expected_line_number_local);
    no_fmt_output_expected += "\n****\n";

    EXPECT_EQ(no_fmt_out_buf.str(), no_fmt_output_expected);

    std::ifstream with_fmt_output(with_fmt);
    std::stringstream with_fmt_out_buf;
    with_fmt_out_buf << with_fmt_output.rdbuf();
    with_fmt_output.close();

    check_level(with_fmt_out_buf.str(), "INFO");
    check_msg(with_fmt_out_buf.str(), "Test outputLocalMessages()");
    check_file(with_fmt_out_buf.str());
    check_line(with_fmt_out_buf.str(), expected_line_number_local);

    // For outputLocalMessages(), only current rank and rank count of 1
    // output for LumberjackStreams. Behaves like SynchronizedStream.
    check_rank(with_fmt_out_buf.str(), rank);
    check_rank_count(with_fmt_out_buf.str(), GetParam(), 1);
  }

  // Closes open file streams associated with Slic streams when destructors
  // called during slic::finalize().
  // Windows _unlink file deletion fails if file is still in use.
  slic::finalize();

  // Cleanup generated files
  EXPECT_EQ(axom::utilities::filesystem::removeFile(no_fmt), 0);
  EXPECT_EQ(axom::utilities::filesystem::removeFile(with_fmt), 0);

  // Clear out ostringstreams for other tests
  slic::flushStreams();
  slic::internal::clear_streams();
}

//------------------------------------------------------------------------------
const std::string parallel_streams[] = {"Synchronized", "Lumberjack"};
INSTANTIATE_TEST_SUITE_P(core_memory_management,
                         SlicMacrosParallel,
                         ::testing::ValuesIn(parallel_streams));

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
