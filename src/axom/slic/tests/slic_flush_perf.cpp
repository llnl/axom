// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/interface/slic.hpp"

#include "axom/slic/streams/LumberjackStream.hpp"
#include "axom/slic/streams/SynchronizedStream.hpp"
#include "gtest/gtest.h"

#include <mpi.h>

// namespace alias
namespace slic = axom::slic;

//------------------------------------------------------------------------------
TEST(SlicPerfTest, test_many_flushes)
{
  // initialize slic
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError(); /* disable abort for testing purposes */

  // Stream for message levels
  std::ostringstream test_stream;

  // Stream for tagged streams
  std::ostringstream test_tag_stream;

  slic::addStreamToAllMsgLevels(
    new slic::LumberjackStream(&test_stream, MPI_COMM_WORLD, 1, "<MESSAGE>\n"));

  slic::addStreamToTag(
    new slic::LumberjackStream(&test_tag_stream, MPI_COMM_WORLD, 1, "<MESSAGE>\n"),
    "myTag");

  slic::addStreamToAllMsgLevels(
    new slic::SynchronizedStream(&test_stream, MPI_COMM_WORLD, "<MESSAGE>\n"));

  slic::addStreamToTag(new slic::SynchronizedStream(&test_tag_stream, MPI_COMM_WORLD, "<MESSAGE>\n"),
                       "myTag");

  for(int i = 0; i < 10000; i++)
  {
    slic::flushStreams();
  }

  slic::finalize();
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