// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <cstring>
#include <string>
#include <vector>

#include "mpi.h"

#include "axom/lumberjack/BinaryTreeCommunicator.hpp"

#include "axom/core/utilities/Utilities.hpp"

TEST(lumberjack_BinaryCommunicator, basic)
{
  MPI_Barrier(MPI_COMM_WORLD);

  const int ranksLimit = 5;
  axom::lumberjack::BinaryTreeCommunicator c;
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
  EXPECT_EQ(c.numPushesToFlush(), axom::utilities::log2(commSize));

  // Push
  std::string s = std::to_string(commRank);
  std::string origS = std::to_string(commRank);

  const char* packedMessage = s.c_str();
  std::vector<const char*> receivedPackedMessages;
  c.push(packedMessage, receivedPackedMessages);
  EXPECT_EQ(strcmp(packedMessage, origS.c_str()), 0);  // Message shouldn't be
                                                       // altered

  // Calculate how many children/messages you should have at this rank
  int numChildren = 0;
  int leftChild = (commRank * 2) + 1;
  int rightChild = (commRank * 2) + 2;
  if(leftChild < commSize)
  {
    numChildren++;
  }
  if(rightChild < commSize)
  {
    numChildren++;
  }
  EXPECT_EQ((int)receivedPackedMessages.size(), numChildren);

  // Verify the messages received were correct
  std::string currMessage = "";
  bool found = false;

  if(numChildren == 2)
  {
    currMessage = std::to_string(rightChild);
    found = false;
    for(auto& rm : receivedPackedMessages)
    {
      if(strcmp(rm, currMessage.c_str()) == 0)
      {
        found = true;
      }
    }
    EXPECT_EQ(found, true) << "Error: Rank: " << commRank
                           << ": Message not received: " << currMessage << std::endl;
  }

  if(numChildren > 0)
  {
    currMessage = std::to_string(leftChild);
    found = false;
    for(auto& rm : receivedPackedMessages)
    {
      if(strcmp(rm, currMessage.c_str()) == 0)
      {
        found = true;
      }
    }
    EXPECT_EQ(found, true) << "Error: Rank: " << commRank
                           << ": Message not received: " << currMessage << std::endl;
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

TEST(lumberjack_BinaryCommunicator, pushNothing)
{
  MPI_Barrier(MPI_COMM_WORLD);

  int commRank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  // Initialize
  axom::lumberjack::BinaryTreeCommunicator c;
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
