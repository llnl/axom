// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/lumberjack/Lumberjack.hpp"

#include <ctime>
#include <iostream>
#include <string>
#include <vector>

class DummyCommunicator : public axom::lumberjack::Communicator
{
public:
  void initialize(MPI_Comm comm, int ranksLimit) override
  {
    m_mpiComm = comm;
    m_ranksLimit = ranksLimit;
    m_isOutputNode = true;
    m_startTime = 0.0;
  }

  void finalize() override { }

  MPI_Comm comm() override { return m_mpiComm; }

  int rank() override { return 0; }

  void ranksLimit(int value) override { m_ranksLimit = value; }

  int ranksLimit() override { return m_ranksLimit; }

  int numPushesToFlush() override { return 1; }

  void push(const char* /* packedMessagesToBeSent */,
            std::vector<const char*>& /* receivedPackedMessages */) override
  { }

  bool isOutputNode() override { return m_isOutputNode; }

  double startTime() override { return m_startTime; }

private:
  MPI_Comm m_mpiComm {MPI_COMM_NULL};
  int m_ranksLimit {0};
  bool m_isOutputNode {true};
  double m_startTime {0.0};
};

class TagOnlyTextCombiner : public axom::lumberjack::Combiner
{
public:
  explicit TagOnlyTextCombiner(std::string tag)
    : m_tag(std::move(tag))
    , m_id(std::string("TagOnlyTextCombiner:") + m_tag)
  { }

  const std::string id() override { return m_id; }

  // NOTE: This method intentionally does not use the `override` keyword so this
  // file can be compiled against older Axom revisions that do not yet declare
  // `Combiner::isMessageCandidateForCombiner()`. In those builds, this method
  // is simply ignored and the benchmark demonstrates the slow path.
  bool isMessageCandidateForCombiner(const axom::lumberjack::Message& message)
  {
    return message.tag() == m_tag;
  }

  bool shouldMessagesBeCombined(const axom::lumberjack::Message& leftMessage,
                                const axom::lumberjack::Message& rightMessage) override
  {
    if(leftMessage.tag() != m_tag || rightMessage.tag() != m_tag)
    {
      return false;
    }

    return leftMessage.text() == rightMessage.text();
  }

  void combine(axom::lumberjack::Message& combined,
               const axom::lumberjack::Message& combinee,
               const int ranksLimit) override
  {
    combined.addRanks(combinee.ranks(), combinee.count(), ranksLimit);
    if(combinee.creationTime() < combined.creationTime())
    {
      combined.creationTime(combinee.creationTime());
    }
  }

private:
  std::string m_tag;
  std::string m_id;
};

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  const int ranksLimit = 5;

  int numPassThroughMessages = 250000;
  int numCombinableMessages = 5000;

  if(argc >= 2)
  {
    numPassThroughMessages = std::stoi(argv[1]);
  }
  if(argc >= 3)
  {
    numCombinableMessages = std::stoi(argv[2]);
  }

  DummyCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);

  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  // Use a single combiner that is only interested in a specific tag. Messages
  // with other tags are "pass through" and should not pay quadratic duplicate
  // checking cost in Lumberjack::combineMessages().
  lumberjack.clearCombiners();
  lumberjack.addCombiner(new TagOnlyTextCombiner("combine_tag"));

  for(int i = 0; i < numPassThroughMessages; ++i)
  {
    lumberjack.queueMessage("pass-through message", "", -1, 0, static_cast<double>(i), "nocombine");
  }

  for(int i = 0; i < numCombinableMessages; ++i)
  {
    lumberjack.queueMessage("combinable message", "", -1, 0, static_cast<double>(i), "combine_tag");
  }

  std::clock_t begin = clock();
  lumberjack.pushMessagesOnce();
  std::clock_t end = clock();

  std::cout << "Elapsed time to push+combine (ms): "
            << ((double)(end - begin) * 1000) / CLOCKS_PER_SEC << std::endl;
  std::cout << "Final message count: " << lumberjack.getMessages().size() << std::endl;

  lumberjack.finalize();
  communicator.finalize();

  return 0;
}
