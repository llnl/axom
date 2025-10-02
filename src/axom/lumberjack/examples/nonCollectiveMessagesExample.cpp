// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/lumberjack.hpp"
#include "axom/lumberjack/NonCollectiveRootCommunicator.hpp"

#include <mpi.h>
#include <unistd.h>
#include <iostream>

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  /*
    This test illustrates the use of a non-collective root communicator 
    to send/receive messages when a runtime error is thrown.
  */
  // Initialize MPI and get rank and comm size
  MPI_Init(&argc, &argv);

  int commRank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  int commSize = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  // Determine how many ranks we want to individually track per message
  int ranksLimit = commSize / 2;

  axom::lumberjack::NonCollectiveRootCommunicator communicator;
  communicator.initialize(MPI_COMM_WORLD, ranksLimit);

  // Initialize lumberjack
  axom::lumberjack::Lumberjack lj;
  lj.initialize(&communicator, ranksLimit);

  if(commRank > 1)
  {
    /* sleep simulates work being done by other ranks 
       that prevents them from sending/receiving messages */
    sleep(5u);
  }

  if(commRank == 1)
  {
    lj.queueMessage("This is a sample error message sent from rank 1");

    lj.pushMessagesOnce();

    // sleep for a bit to allow some time for the message to be sent
    sleep(3u);

    /* This is typically the point where rank 1 calls
       MPI_abort and exits. It does not wait for other ranks 
       because they may be hanging or doing other work. 
       To simulate the program exiting, we will put a 
       print statement below. */
    std::cout << "--> simulated program termination" << std::endl;
  }

  if(commRank == 0)
  {
    /* Call push to determine whether error messages have 
       arrived from any rank.  In this example, we expect an
       error message to arrive from rank 1.  In applications, 
       this check should be done at a regular interval to see if 
       any rank has reported an error.
    */
    lj.pushMessagesOnce();
    std::vector<axom::lumberjack::Message*> messages = lj.getMessages();
    for(int i = 0; i < (int)(messages.size()); ++i)
    {
      // Output a single Message at a time to screen
      std::cout << "(" << messages[i]->stringOfRanks() << ") " << messages[i]->count() << " '"
                << messages[i]->text() << "'" << std::endl;
    }
    lj.clearMessages();
  }

  // Finalize lumberjack
  lj.finalize();
  // Finalize the lumberjack communicator
  communicator.finalize();
  // Finalize MPI
  MPI_Finalize();

  return 0;
}
