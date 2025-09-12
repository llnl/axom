// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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

//------------------------------------------------------------------------------
int main()
{
  int ranksLimit = 5;
  DummyCommunicator communicator;
  communicator.initialize(MPI_COMM_NULL, ranksLimit);
  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  // Remove default combiner (no combiners now)
  lumberjack.removeCombiner("TextTagCreationTimeCombiner");

  for(int i = 0; i < 100000; i++)
  {
    lumberjack.queueMessage("Should not be combined.", static_cast<double>(i));
  }

  std::clock_t begin = clock();
  lumberjack.pushMessagesOnce();
  std::clock_t end = clock();

  lumberjack.finalize();
  communicator.finalize();

  std::cout << "Elapsed time to push messages (ms): "
            << ((double)(end - begin) * 1000) / CLOCKS_PER_SEC << std::endl;

  return 0;
}
