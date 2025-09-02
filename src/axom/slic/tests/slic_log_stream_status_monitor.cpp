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

enum class StreamType
{
    Generic,
    Synchronized,
    Lumberjack
};

class LogStreamStatusMonitorParamTest : public ::testing::TestWithParam<StreamType>
{
protected:
    std::ostringstream test_stream;
    std::unique_ptr<axom::slic::LogStream> stream;
    axom::slic::LogStreamStatusMonitor logStreamStatusMonitor;

    void SetUp() override
    {
      switch(GetParam())
      {
        case StreamType::Generic:
          stream = std::make_unique<axom::slic::GenericOutputStream>(&test_stream);
          break;
        case StreamType::Synchronized:
          stream = std::make_unique<axom::slic::SynchronizedStream>(&test_stream, MPI_COMM_WORLD);
          break;
        case StreamType::Lumberjack:
          stream = std::make_unique<axom::slic::LumberjackStream>(&test_stream, MPI_COMM_WORLD, 1);
          break;
      }
    }
};

//------------------------------------------------------------------------------
TEST_P(LogStreamStatusMonitorParamTest, test_has_pending_messages)
{
  logStreamStatusMonitor.addStream(stream.get());

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  stream->append(axom::slic::message::Debug, "test message", "test tag", "test file name", 1, false, false);

  if(stream->canHavePendingMessages() == true)
  {
    EXPECT_TRUE(logStreamStatusMonitor.hasPendingMessages());
  }
  else
  {
    EXPECT_FALSE(logStreamStatusMonitor.hasPendingMessages());
  }

  stream->flush();

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

}

INSTANTIATE_TEST_SUITE_P(
    StreamTypes,
    LogStreamStatusMonitorParamTest,
    ::testing::Values(
        StreamType::Generic,
        StreamType::Synchronized,
        StreamType::Lumberjack
    )
);

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

  auto ljstream = axom::slic::LumberjackStream(&test_stream, &lj);

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
