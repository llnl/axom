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

  if(stream->isUsingMPI() == true)
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

INSTANTIATE_TEST_SUITE_P(StreamTypes,
                         LogStreamStatusMonitorParamTest,
                         ::testing::Values(StreamType::Generic,
                                           StreamType::Synchronized,
                                           StreamType::Lumberjack));

//------------------------------------------------------------------------------
TEST(SlicLogStreamMonitorTest, test_add_streams_different_comms)
{
  /*
    This test checks hasPendingMessages when different 
    ranks have different MPI communicators.
  */

  MPI_Group world_group;
  MPI_Comm comm0 = MPI_COMM_NULL, comm1 = MPI_COMM_NULL;

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  const bool rank_is_even = (rank % 2 == 0);

  int color = (rank_is_even) ? 0 : 1;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm0);

  // To create communicator asymmetry, rank 0 is included in both even and odd communicators.
  if(rank == 0)
  {
    color = 1;
  }

  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm1);

  std::ostringstream test_stream;

  axom::slic::LogStreamStatusMonitor logStreamStatusMonitor;

  auto ljstream0 = axom::slic::LumberjackStream(&test_stream, comm0, 1);
  logStreamStatusMonitor.addStream(&ljstream0);

  auto ljstream1 = axom::slic::LumberjackStream(&test_stream, comm1, 1);
  logStreamStatusMonitor.addStream(&ljstream1);

  if(rank_is_even)
  {
    ljstream0.append(axom::slic::message::Debug,
                     "test message",
                     "test tag",
                     "test file name",
                     1,
                     false,
                     false);
  }

  /*
    All ranks should have pending messages because logStreamStatusMonitor
    calls MPI_Allreduce on all MPI communicators
  */
  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), true);

  ljstream0.flush();
  ljstream1.flush();

  EXPECT_EQ(logStreamStatusMonitor.hasPendingMessages(), false);

  if(comm0 != MPI_COMM_NULL)
  {
    MPI_Comm_free(&comm0);
  }

  if(comm1 != MPI_COMM_NULL)
  {
    MPI_Comm_free(&comm1);
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
