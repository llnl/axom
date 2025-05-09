// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <cstring>
#include <string>
#include <vector>

#include "mpi.h"

#include "axom/lumberjack/RootCommunicator.hpp"

TEST(lumberjack_RootCommunicator, basic)
{
  MPI_Barrier(MPI_COMM_WORLD);

  const int ranksLimit = 5;
  axom::lumberjack::RootCommunicator c;
  c.initialize(MPI_COMM_WORLD, ranksLimit);

  int commRank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  int commSize = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  // Check rank
  EXPECT_EQ(c.rank(), commRank);

  // Check if we are an output node
  if(commRank != 0)
  {
    EXPECT_EQ(c.isOutputNode(), false);
  }
  else
  {
    EXPECT_EQ(c.isOutputNode(), true);
  }

  // Ranks limit get/set
  EXPECT_EQ(c.ranksLimit(), ranksLimit);
  const int newRanksLimit = 98876;
  c.ranksLimit(newRanksLimit);
  EXPECT_EQ(c.ranksLimit(), newRanksLimit);

  // Check number of pushes
  EXPECT_EQ(c.numPushesToFlush(), 1);

  // Push
  std::string s = std::to_string(commRank);
  std::string origS = std::to_string(commRank);

  const char* packedMessage = s.c_str();
  std::vector<const char*> receivedPackedMessages;
  c.push(packedMessage, receivedPackedMessages);
  EXPECT_EQ(strcmp(packedMessage, origS.c_str()), 0);

  if(commRank != 0)
  {
    EXPECT_EQ((int)receivedPackedMessages.size(), 0);
  }
  else
  {
    const int numMessagesToReceive = commSize - 1;
    EXPECT_EQ((int)receivedPackedMessages.size(), numMessagesToReceive);
    for(int i = 1; i <= numMessagesToReceive; ++i)
    {
      std::string currMessage = std::to_string(i);
      bool found = false;
      for(auto& rm : receivedPackedMessages)
      {
        if(strcmp(rm, currMessage.c_str()) == 0)
        {
          found = true;
        }
      }
      if(!found)
      {
        std::cout << "Error: Message not received:" << currMessage << std::endl;
      }
      EXPECT_EQ(found, true);
    }
  }

  c.finalize();

  MPI_Barrier(MPI_COMM_WORLD);

  // cleanup allocated memory from received messages
  for(auto& rm : receivedPackedMessages)
  {
    delete[] rm;
  }
  receivedPackedMessages.clear();
}

TEST(lumberjack_RootCommunicator, pushNothing)
{
  MPI_Barrier(MPI_COMM_WORLD);

  int commRank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  // Initialize
  axom::lumberjack::RootCommunicator c;
  c.initialize(MPI_COMM_WORLD, 5);
  std::vector<const char*> receivedPackedMessages;

  // Call push when there is empty message
  std::string emptyString = "";
  std::string origEmptyString = "";
  const char* packedMessage = emptyString.c_str();
  c.push(packedMessage, receivedPackedMessages);
  EXPECT_EQ(strcmp(packedMessage, origEmptyString.c_str()), 0);
  EXPECT_EQ((int)receivedPackedMessages.size(), 0);

  // Call push with a nullptr
  packedMessage = nullptr;
  c.push(packedMessage, receivedPackedMessages);
  EXPECT_EQ(packedMessage, nullptr);  // Message should still be nullptr
  EXPECT_EQ((int)receivedPackedMessages.size(), 0);

  // Finalize
  c.finalize();

  MPI_Barrier(MPI_COMM_WORLD);
}
