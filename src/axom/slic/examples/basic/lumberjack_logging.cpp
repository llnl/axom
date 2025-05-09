// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// C/C++ includes
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// axom includes
#include "axom/core/Types.hpp"
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/LumberjackStream.hpp"

// MPI
#include <mpi.h>

using namespace axom;

#define CYCLELIMIT 5
#define RANKSLIMIT 5

#define N 20

slic::message::Level getRandomEvent(const int start, const int end)
{
  return (static_cast<slic::message::Level>(std::rand() % (end - start) + start));
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Initialize MPI
  MPI_Init(&argc, &argv);
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Initialize SLIC
  constexpr const char* format = R"(
<MESSAGE>
\t<TIMESTAMP>
\tLEVEL=<LEVEL>
\tRANKS=<RANK>
\tRANK_COUNT=<RANK_COUNT>
\tFILE=<FILE>
\tLINE=<LINE>
)";
  slic::initialize();

  // Set SLIC logging level and Lumberjack Logging stream
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError();
  slic::LumberjackStream* ljStream =
    new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RANKSLIMIT, format);
  slic::addStreamToAllMsgLevels(ljStream);

  // Queue messages
  int cycleCount = 0;
  for(int i = 0; i < N; ++i)
  {
    std::ostringstream oss;
    oss << "message " << i << "/" << N - 1;

    slic::logMessage(getRandomEvent(0, slic::message::Num_Levels), oss.str(), __FILE__, __LINE__);

    ++cycleCount;
    if(cycleCount > CYCLELIMIT)
    {
      // Incrementally push messages through the log stream
      slic::pushStreams();
      cycleCount = 0;
    }
  }

  // Fully flush system of messages
  slic::flushStreams();

  // Shutdown SLIC which in turn shutsdown Lumberjack
  slic::finalize();

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
