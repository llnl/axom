// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/lumberjack/Lumberjack.hpp"

#include <ctime>
#include <iostream>

class DummyCommunicator : public axom::lumberjack::Communicator
{
public:
  void initialize(MPI_Comm comm, int ranksLimit)
  {
    m_mpiComm = comm;
    m_ranksLimit = ranksLimit;
    m_isOutputNode = true;
    srand(time(nullptr));
    m_startTime = 0.0;
  }

  void finalize() { }

  MPI_Comm comm() { return m_mpiComm; }

  int rank() { return 0; }

  void ranksLimit(int value) { m_ranksLimit = value; }

  int ranksLimit() { return m_ranksLimit; }

  int numPushesToFlush() { return 1; }

  void push(const char* /* packedMessagesToBeSent */,
            std::vector<const char*>& /* receivedPackedMessages */)
  { }

  bool isOutputNode() { return m_isOutputNode; }

  void outputNode(bool value) { m_isOutputNode = value; }

  double startTime() { return m_startTime; }

private:
  MPI_Comm m_mpiComm;
  int m_ranksLimit;
  bool m_isOutputNode;
  double m_startTime;
};

// Inherits from TextTagCombiner - only messages with tag "tagA" are candidates
class CombinerTagA : public axom::lumberjack::TextTagCombiner
{
public:
  const std::string id() { return m_id; }

  bool isMessageCandidateForCombiner(const axom::lumberjack::Message& m)
  {
    return m.tag() == "tagA";
  }

private:
  const std::string m_id = "TagACombiner";
};

//------------------------------------------------------------------------------
int main()
{
  int ranksLimit = 5;
  DummyCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  // Remove default combiner, use custom combiner
  lumberjack.removeCombiner("TextTagCombiner");
  lumberjack.addCombiner(new CombinerTagA);

  // Combinable message
  lumberjack.queueMessage("Should be combined.", "", 0, 0, 0.0, "tagA");

  // Unique, uncombinable messages (different text)
  for(int i = 0; i < 10000; i++)
  {
    lumberjack.queueMessage("Unique message " + std::to_string(i) + " is not combinable ",
                            "",
                            0,
                            0,
                            1.0,
                            "ignoreTag");
  }
  lumberjack.queueMessage("Should be combined.", "", 0, 0, 0.0, "tagA");

  std::clock_t begin = clock();
  lumberjack.pushMessagesOnce();
  std::clock_t end = clock();

  std::cout << "After combining, there are " << lumberjack.getMessages().size() << " messages."
            << std::endl;

  lumberjack.finalize();
  communicator.finalize();

  std::cout << "Elapsed time to push messages (ms): "
            << ((double)(end - begin) * 1000) / CLOCKS_PER_SEC << std::endl;

  return 0;
}
