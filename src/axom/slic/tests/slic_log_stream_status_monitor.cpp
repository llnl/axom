// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/lumberjack/Lumberjack.hpp"
#include "axom/lumberjack/RootCommunicator.hpp"
#include "axom/slic/core/LogStreamStatusMonitorMPI.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"
#include "axom/slic/streams/LumberjackStream.hpp"
#include "axom/slic/streams/SynchronizedStream.hpp"
#include "gtest/gtest.h"

#include <mpi.h>

// namespace alias
namespace slic = axom::slic;

TEST(SlicLogStreamMonitorTest, test_has_pending_messages)
{

  std::ostringstream test_stream;

  auto generic_stream = slic::GenericOutputStream(&test_stream);

  axom::slic::LogStreamStatusMonitor logStreamStatusMonitor;
  axom::slic::LogStreamStatusMonitorMPI logStreamStatusMonitorMPI;

  logStreamStatusMonitor.addStream(&generic_stream);
  logStreamStatusMonitorMPI.addStream(&generic_stream);

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);
  EXPECT_EQ(logStreamStatusMonitorMPI.hasPendingMessages(), false);

  generic_stream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  /* hasPendingMessages will not probe for unflushed messages in 
     GenericOutputStream because it relies on ofstream for flushing rather than MPI 
   */
  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);
  EXPECT_EQ(logStreamStatusMonitorMPI.hasPendingMessages(), false);

  generic_stream.flush();

  auto synchronized_stream = slic::SynchronizedStream(&test_stream, MPI_COMM_WORLD);

  logStreamStatusMonitorMPI.addStream(&synchronized_stream);

  EXPECT_EQ(logStreamStatusMonitorMPI.hasPendingMessages(), false);

  synchronized_stream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  EXPECT_EQ(logStreamStatusMonitorMPI.hasPendingMessages(), true);

  synchronized_stream.flush();

  auto lj_stream = slic::LumberjackStream(&test_stream, MPI_COMM_WORLD, 1);

  logStreamStatusMonitorMPI.addStream(&lj_stream);

  EXPECT_EQ(logStreamStatusMonitorMPI.hasPendingMessages(), false);

  lj_stream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  EXPECT_EQ(logStreamStatusMonitorMPI.hasPendingMessages(), true);

  lj_stream.flush();

  logStreamStatusMonitor.finalize();
  logStreamStatusMonitorMPI.finalize();

}

//------------------------------------------------------------------------------
TEST(SlicLogStreamMonitorTest, test_add_streams_multiple_comms)
{

  MPI_Group world_group, group;
  MPI_Comm comm = MPI_COMM_NULL;

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  // Split ranks into even and odd
  std::vector<int> ranks;
  if (size == 1)
  {
    ranks.push_back(0);
  }
  else
  {
    for(int i = 0; i < size; ++i)
    {
      if(i % 2 == 0 && rank % 2 == 0)
      {
        ranks.push_back(i);
      }

      if (i % 2 != 0 && rank % 2 != 0)
      {
        ranks.push_back(i);
      }
    }
  }

  MPI_Group_incl(world_group, ranks.size(), ranks.data(), &group);
  MPI_Comm_create(MPI_COMM_WORLD, group, &comm);

  MPI_Group_free(&world_group);
  MPI_Group_free(&group);

  std::ostringstream test_stream;

  auto lj_comm = axom::lumberjack::RootCommunicator();

  auto lj = axom::lumberjack::Lumberjack();
  lj_comm.initialize(comm, 1);
  lj.initialize(&lj_comm, 1);

  auto ljstream = slic::LumberjackStream(&test_stream, &lj);

  axom::slic::LogStreamStatusMonitorMPI logStreamStatusMonitor;

  if(rank % 2 == 0)
  {
    logStreamStatusMonitor.addStream(&ljstream);
  }

  MPI_Barrier(comm);

  /* test that no messages should exist */
  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  if(rank % 2 == 0)
  {
    ljstream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);
  }

  /*
    test that pending messages should only exist on even ranks
    when logging messages from a stream whose communicator contains
    only even ranks
   */
  if(rank % 2 == 0)
  {
    EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), true);
  }
  else
  {
    EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);
  }

  if(rank % 2 != 0)
  {
    logStreamStatusMonitor.addStream(&ljstream);
  }
  
  MPI_Barrier(comm);

  /*
    test that pending messages should only exist on even ranks
    when logging messages from a stream whose communicator contains
    only even ranks
   */
  if(rank % 2 == 0)
  {
    EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), true);
  }
  else
  {
    EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);
  }

  ljstream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  MPI_Barrier(comm);

  /*
    Test that pending messages exist on all ranks
   */
  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), true);

  ljstream.flush();

  /*
    Test that no pending messages exist on all ranks
   */
  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  if(comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&comm);
  }

  logStreamStatusMonitor.finalize();

}

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