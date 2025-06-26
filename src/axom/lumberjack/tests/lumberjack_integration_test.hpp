// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/lumberjack/Lumberjack.hpp"
#include "axom/lumberjack/BinaryTreeCommunicator.hpp"
#include "axom/lumberjack/Message.hpp"

#include <mpi.h>
#include <iomanip>

TEST(lumberjack_Lumberjack, sortMessages)
{
  MPI_Barrier(MPI_COMM_WORLD);

  int commSize;
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  int commRank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  const int ranksLimit = 5;
  axom::lumberjack::BinaryTreeCommunicator communicator;
  communicator.initialize(MPI_COMM_WORLD, ranksLimit);

  axom::lumberjack::Lumberjack lumberjack;
  lumberjack.initialize(&communicator, ranksLimit);

  if (commRank == 0) {
    lumberjack.queueMessage("entering some interactive input command");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::string output_msg = "outputting results from the interactive input command command on rank " + std::to_string(commRank);
  lumberjack.queueMessage(output_msg);

  lumberjack.pushMessagesFully();

  if (commRank == 0) {
    std::vector<axom::lumberjack::Message*> messages = lumberjack.getMessages();

    EXPECT_EQ((int)messages.size(), commSize+1);

    EXPECT_EQ(messages[0]->text(), "entering some interactive input command");

    EXPECT_TRUE(messages[0]->creationTime() <= messages[1]->creationTime());

    for (int i = 1; i < commSize; i++) {
      EXPECT_TRUE(messages[i]->creationTime() <= messages[i+1]->creationTime());
    }
  }


  lumberjack.finalize();
  communicator.finalize();

  MPI_Barrier(MPI_COMM_WORLD);
}