// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/lumberjack/Lumberjack.hpp"
#include "axom/lumberjack/RootCommunicator.hpp"
#include "axom/slic/core/LogStreamStatusMonitor.hpp"
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

  logStreamStatusMonitor.addStream(&generic_stream);

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  generic_stream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  /* hasPendingMessages will not probe for unflushed messages in 
     GenericOutputStream because it relies on ofstream for flushing rather than MPI 
   */
  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  generic_stream.flush();

  auto synchronized_stream = slic::SynchronizedStream(&test_stream, MPI_COMM_WORLD);

  logStreamStatusMonitor.addStream(&synchronized_stream);

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  synchronized_stream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), true);

  synchronized_stream.flush();

  auto lj_stream = slic::LumberjackStream(&test_stream, MPI_COMM_WORLD, 1);

  logStreamStatusMonitor.addStream(&lj_stream);

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  lj_stream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), true);

  lj_stream.flush();

}

//------------------------------------------------------------------------------
TEST(SlicLogStreamMonitorTest, test_add_streams_different_comms)
{

  /*
    This test checks hasPendingMessages when different 
    ranks have different MPI communicators.
  */

  MPI_Group world_group, group;
  MPI_Comm comm = MPI_COMM_NULL;

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  bool rank_is_even = (rank % 2 == 0);

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
      if((i % 2 == 0) == rank_is_even)
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

  axom::slic::LogStreamStatusMonitor logStreamStatusMonitor;

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  logStreamStatusMonitor.addStream(&ljstream);

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  ljstream.append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), true);

  ljstream.flush();

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  lj.finalize();
  lj_comm.finalize();

  if(comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&comm);
  }

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
