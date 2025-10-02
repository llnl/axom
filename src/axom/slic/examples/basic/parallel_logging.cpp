// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// C/C++ includes
#include <cstdlib>  // for rand()
#include <sstream>  // for ostringstream

// Logging includes
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/SynchronizedStream.hpp"

// MPI
#include <mpi.h>

using namespace axom;

#define N 20

slic::message::Level getRandomEvent(const int start, const int end)
{
  return (static_cast<slic::message::Level>(std::rand() % (end - start) + start));
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // STEP 0: initialize MPI & logging environment
  MPI_Init(&argc, &argv);

  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string format = std::string("[<RANK>]: <MESSAGE>\n") + std::string("\t<TIMESTAMP>\n") +
    std::string("\tLEVEL=<LEVEL>\n") + std::string("\tFILE=<FILE>\n") +
    std::string("\tLINE=<LINE>\n");

  slic::initialize();

  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError();
  slic::addStreamToAllMsgLevels(new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, format));

  // STEP 3: loop N times and generate a random logging event
  for(int i = 0; i < N; ++i)
  {
    std::ostringstream oss;
    oss << "message " << i << "/" << N - 1;

    slic::logMessage(getRandomEvent(0, slic::message::Num_Levels), oss.str(), __FILE__, __LINE__);

    // Flush every 5 cycles
    if((i % 5) == 0)
    {
      slic::flushStreams();

    }  // END if
  }

  // STEP 4: shutdown logging environment
  slic::finalize();

  // STEP 5: Finalize MPI
  MPI_Finalize();

  return 0;
}
