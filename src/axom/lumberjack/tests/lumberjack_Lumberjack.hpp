// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/lumberjack/Lumberjack.hpp"

#include "axom/lumberjack/Communicator.hpp"
#include "axom/lumberjack/Message.hpp"

#include <stdlib.h>
#include <time.h>

class TestCommunicator : public axom::lumberjack::Communicator
{
public:
  void initialize(MPI_Comm comm, int ranksLimit)
  {
    m_mpiComm = comm;
    m_ranksLimit = ranksLimit;
    m_isOutputNode = true;
    srand(time(nullptr));
  }

  void finalize() { }

  int rank() { return rand() % (m_ranksLimit * 4); }

  void ranksLimit(int value) { m_ranksLimit = value; }

  int ranksLimit() { return m_ranksLimit; }

  int numPushesToFlush() { return 1; }

  void push(const char* /* packedMessagesToBeSent */,
            std::vector<const char*>& /* receivedPackedMessages */)
  { }

  bool isOutputNode() { return m_isOutputNode; }

  void outputNode(bool value) { m_isOutputNode = value; }

private:
  MPI_Comm m_mpiComm;
  int m_ranksLimit;
  bool m_isOutputNode;
};

TEST(lumberjack_Lumberjack, combineMessagesPushOnce01)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");

  lumberjack.pushMessagesOnce();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), 1);
  EXPECT_EQ(messages[0]->text(), "Should be combined.");
  EXPECT_EQ(messages[0]->count(), 6);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessagesPushOnce02)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");

  lumberjack.pushMessagesOnce();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), 2);
  EXPECT_EQ(messages[0]->text(), "");
  EXPECT_EQ(messages[0]->count(), 1);
  EXPECT_EQ(messages[1]->text(), "Should be combined.");
  EXPECT_EQ(messages[1]->count(), 5);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessagesPushOnceEmpty)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("");

  lumberjack.pushMessagesOnce();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), 1);
  EXPECT_EQ(messages[0]->text(), "");
  EXPECT_EQ(messages[0]->count(), 1);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessagesPushOnceEmptyNonOutputNode)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  communicator.outputNode(false);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("");

  lumberjack.pushMessagesOnce();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), 0);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessagesPushOnceNothingNonOutputNode)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  communicator.outputNode(false);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.pushMessagesOnce();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), 0);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessages01)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");

  lumberjack.pushMessagesFully();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), 1);
  EXPECT_EQ(messages[0]->text(), "Should be combined.");
  EXPECT_EQ(messages[0]->count(), 6);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessages02)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("Should be combined 1.");
  lumberjack.queueMessage("Should be combined 1.");
  lumberjack.queueMessage("Should be combined 1.");
  lumberjack.queueMessage("Should be combined 2.");
  lumberjack.queueMessage("Should be combined 2.");
  lumberjack.queueMessage("Should be combined 2.");

  lumberjack.pushMessagesFully();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), 2);
  EXPECT_EQ(messages[0]->text(), "Should be combined 1.");
  EXPECT_EQ(messages[0]->count(), 3);
  EXPECT_EQ(messages[1]->text(), "Should be combined 2.");
  EXPECT_EQ(messages[1]->count(), 3);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessages03)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("Should not be combined 1.");
  lumberjack.queueMessage("Should not be combined 2.");
  lumberjack.queueMessage("Should not be combined 3.");
  lumberjack.queueMessage("Should not be combined 4.");
  lumberjack.queueMessage("Should not be combined 5.");
  lumberjack.queueMessage("Should not be combined 6.");

  lumberjack.pushMessagesFully();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), 6);
  for(int i = 0; i < 6; ++i)
  {
    std::string s = "Should not be combined " + std::to_string(i + 1) + ".";
    EXPECT_EQ(messages[i]->text(), s);
    EXPECT_EQ(messages[i]->count(), 1);
  }

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMixedMessages01)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("Should not be combined 1.");
  lumberjack.queueMessage("Should not be combined 2.");
  lumberjack.queueMessage("");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");

  lumberjack.pushMessagesFully();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  // Check total messages size
  EXPECT_EQ((int)messages.size(), 4);
  for(int i = 0; i < 2; ++i)
  {
    std::string s = "Should not be combined " + std::to_string(i + 1) + ".";
    EXPECT_EQ(messages[i]->text(), s);
    EXPECT_EQ(messages[i]->count(), 1);
  }

  EXPECT_EQ(messages[2]->text(), "");
  EXPECT_EQ(messages[2]->count(), 1);

  EXPECT_EQ(messages[3]->text(), "Should be combined.");
  EXPECT_EQ(messages[3]->count(), 3);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMixedMessages02)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("Should not be combined 1.");
  lumberjack.queueMessage("Should not be combined 2.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("");

  lumberjack.pushMessagesFully();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  // Check total messages size
  EXPECT_EQ((int)messages.size(), 4);
  for(int i = 0; i < 2; ++i)
  {
    std::string s = "Should not be combined " + std::to_string(i + 1) + ".";
    EXPECT_EQ(messages[i]->text(), s);
    EXPECT_EQ(messages[i]->count(), 1);
  }

  EXPECT_EQ(messages[2]->text(), "Should be combined.");
  EXPECT_EQ(messages[2]->count(), 3);

  EXPECT_EQ(messages[3]->text(), "");
  EXPECT_EQ(messages[3]->count(), 1);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMixedMessages03)
{
  int ranksLimit = 5;
  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  lumberjack.queueMessage("");
  lumberjack.queueMessage("Should not be combined 2.");
  lumberjack.queueMessage("Should not be combined 3.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");
  lumberjack.queueMessage("Should be combined.");

  lumberjack.pushMessagesFully();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  // Check total messages size
  EXPECT_EQ((int)messages.size(), 4);

  EXPECT_EQ(messages[0]->text(), "");
  EXPECT_EQ(messages[0]->count(), 1);

  for(int i = 1; i < 3; ++i)
  {
    std::string s = "Should not be combined " + std::to_string(i + 1) + ".";
    EXPECT_EQ(messages[i]->text(), s);
    EXPECT_EQ(messages[i]->count(), 1);
  }

  EXPECT_EQ(messages[3]->text(), "Should be combined.");
  EXPECT_EQ(messages[3]->count(), 3);

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessagesManyMessages)
{
  int ranksLimit = 5;
  const int loopCount = 10000;

  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  for(int i = 0; i < loopCount; ++i)
  {
    std::string s = "Should not be combined " + std::to_string(i) + ".";
    lumberjack.queueMessage(s);
  }

  lumberjack.pushMessagesFully();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), loopCount);
  for(int i = 0; i < loopCount; ++i)
  {
    std::string s = "Should not be combined " + std::to_string(i) + ".";
    EXPECT_EQ(messages[i]->text(), s);
    EXPECT_EQ(messages[i]->count(), 1);
  }

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, combineMessagesLargeMessages)
{
  int ranksLimit = 5;
  const int loopCount = 10;
  const int padSize = 1000;
  std::string padding = "";
  for(int j = 0; j < padSize; ++j)
  {
    padding += "0";
  }

  TestCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  for(int i = 0; i < loopCount; ++i)
  {
    std::string s = std::to_string(i) + ":" + padding;
    lumberjack.queueMessage(s);
  }

  lumberjack.pushMessagesFully();

  std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

  EXPECT_EQ((int)messages.size(), loopCount);
  for(int i = 0; i < loopCount; ++i)
  {
    std::string s = std::to_string(i) + ":" + padding;
    EXPECT_EQ(messages[i]->text(), s);
    EXPECT_EQ(messages[i]->count(), 1);
  }

  lumberjack.finalize();
  communicator.finalize();
}

TEST(lumberjack_Lumberjack, setNonOwnedCommunicator)
{
  int ranksLimit = 5;
  auto communicator1 = new TestCommunicator();
  communicator1->initialize(MPI_COMM_NULL, ranksLimit);

  auto communicator2 = new TestCommunicator();
  communicator2->initialize(MPI_COMM_NULL, ranksLimit);

  axom::lumberjack::Lumberjack lumberjack;

  lumberjack.initialize(communicator1, ranksLimit);

  EXPECT_EQ(communicator1, lumberjack.getCommunicator());
  EXPECT_EQ(lumberjack.isCommunicatorOwned(), false);

  lumberjack.setCommunicator(communicator2, false);

  /* communicator1 should still be valid after set
     because it's not owned by Lumberjack*/
  EXPECT_NE(communicator1, nullptr);
  EXPECT_EQ(communicator2, lumberjack.getCommunicator());

  lumberjack.finalize();
  communicator1->finalize();
  communicator2->finalize();
  delete communicator1;
  delete communicator2;
}

TEST(lumberjack_Lumberjack, setOwnedCommunicator)
{
  int ranksLimit = 5;
  auto communicator = new TestCommunicator();
  communicator->initialize(MPI_COMM_NULL, ranksLimit);

  axom::lumberjack::Lumberjack lumberjack;

  lumberjack.initialize(new TestCommunicator(), ranksLimit, true);

  EXPECT_NE(lumberjack.getCommunicator(), nullptr);
  EXPECT_EQ(lumberjack.isCommunicatorOwned(), true);

  lumberjack.setCommunicator(communicator, true);

  EXPECT_EQ(lumberjack.getCommunicator(), communicator);
  EXPECT_EQ(lumberjack.isCommunicatorOwned(), true);

  lumberjack.finalize();
}

TEST(lumberjack_Lumberjack, setOwnedAndNonOwnedCommunicator)
{
  int ranksLimit = 5;
  auto communicator = new TestCommunicator();
  communicator->initialize(MPI_COMM_NULL, ranksLimit);

  axom::lumberjack::Lumberjack lumberjack;

  lumberjack.initialize(new TestCommunicator(), ranksLimit, true);

  EXPECT_NE(lumberjack.getCommunicator(), nullptr);
  EXPECT_EQ(lumberjack.isCommunicatorOwned(), true);

  lumberjack.setCommunicator(communicator, false);

  EXPECT_EQ(lumberjack.getCommunicator(), communicator);
  EXPECT_EQ(lumberjack.isCommunicatorOwned(), false);

  lumberjack.finalize();
  communicator->finalize();
  delete communicator;
}
