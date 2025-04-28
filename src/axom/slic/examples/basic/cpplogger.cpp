// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// C/C++ includes
#include <fstream>  // for file stream
#include <cstdlib>  // for rand()

#include "axom/core.hpp"

// Logging includes
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

using namespace axom;
#define N 10

slic::message::Level getRandomEvent(const int start, const int end)
{
  return (static_cast<slic::message::Level>(std::rand() % (end - start) + start));
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  //----------------------------------------------------------------------------
  // STEP 0: Initialize logger
  //----------------------------------------------------------------------------
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError();

  //----------------------------------------------------------------------------
  // STEP 1: Create log streams
  //----------------------------------------------------------------------------

  // setup log stream for ERROR and WARNING messages
  std::ofstream hspStream;
  hspStream.open("HSP.dat");

  std::string hsp_format = std::string("***********************************\n") +
    std::string("* <TIMESTAMP>\n\n") + std::string("* LEVEL=<LEVEL>\n") +
    std::string("* MESSAGE=<MESSAGE>\n") + std::string("* FILE=<FILE>\n") +
    std::string("* LINE=<LINE>\n") + std::string("***********************************\n");

  slic::LogStream* hspLogStream = new slic::GenericOutputStream(&hspStream, hsp_format);

  // setup log stream for ALL messages, including ERROR and WARNING
  std::string console_format = std::string("[<LEVEL>]: <MESSAGE>\n");
  slic::LogStream* console = new slic::GenericOutputStream(&std::cerr, console_format);

  //----------------------------------------------------------------------------
  // STEP 2: add streams to logger
  //----------------------------------------------------------------------------
  slic::addStreamToMsgLevel(hspLogStream, slic::message::Error);
  slic::addStreamToMsgLevel(hspLogStream, slic::message::Warning);

  slic::addStreamToAllMsgLevels(console);

  //----------------------------------------------------------------------------
  // STEP 3: Loop N times and generate random logging events
  //----------------------------------------------------------------------------
  for(int i = 0; i < N; ++i)
  {
    slic::logMessage(getRandomEvent(0, slic::message::Num_Levels),
                     "a random message",
                     __FILE__,
                     __LINE__);
  }

  //----------------------------------------------------------------------------
  // STEP 4: Loop N times and generate random logging events
  //----------------------------------------------------------------------------

  // finalize logging & close HSP file
  slic::finalize();
  hspStream.close();

  return 0;
}
