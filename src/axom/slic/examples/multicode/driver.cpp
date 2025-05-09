// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// C/C++ includes
#include <fstream>  // for filestream
#include <cstdlib>  // for rand()
#include <sstream>  // for string stream

#include "physicsA.hpp"
#include "physicsB.hpp"
#include "axom/config.hpp"
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

#define N 10

using namespace axom;

std::ofstream hsp;

//-----------------------------------------------------------------------------
int getRandInt(const int start, const int end) { return (std::rand() % (end - start) + start); }

//-----------------------------------------------------------------------------
slic::message::Level getRandomLevel()
{
  return (static_cast<slic::message::Level>(getRandInt(0, slic::message::Num_Levels)));
}

//-----------------------------------------------------------------------------
void driver_init()
{
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError();

  std::string hsp_format = std::string("[<LEVEL>]: <MESSAGE>\n") + std::string("\t FILE:<FILE>\n") +
    std::string("\t LINE:<LINE>\n");

  // setup main hsp output
  hsp.open("hsp.log");
  slic::LogStream* ls = new slic::GenericOutputStream(&hsp, hsp_format);
  slic::addStreamToAllMsgLevels(ls);

  std::string console_format = std::string("[<LEVEL>]: <MESSAGE>\n");

  // setup console output
  slic::LogStream* console = new slic::GenericOutputStream(&std::cout, console_format);
  slic::addStreamToMsgLevel(console, slic::message::Error);
  slic::addStreamToMsgLevel(console, slic::message::Warning);
  slic::addStreamToMsgLevel(console, slic::message::Info);
}

//-----------------------------------------------------------------------------
void driver_finalize()
{
  slic::finalize();
  hsp.close();
}

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  driver_init();
  physicsA::init();
  physicsB::init();

  std::ostringstream oss;
  oss.clear();

  for(int t = 0; t < N; ++t)
  {
    // STEP 0: Execute hydro
    slic::logMessage(slic::message::Info, "== Executing physicsA == \n");
    physicsA::timestep(t, getRandInt(2, N));

    slic::message::Level random = getRandomLevel();
    oss.str("");
    oss << "t=" << t << " a random message!";
    slic::logMessage(random, oss.str(), __FILE__, __LINE__);

    slic::logMessage(slic::message::Info, "== Done Executing physicsA ==\n");

    // STEP 1: Execute
    slic::logMessage(slic::message::Info, "== Executing physicsB ==\n");
    physicsB::timestep(t, getRandInt(2, N));

    slic::logMessage(slic::message::Info, "== Done Executing physicsB ==\n");
  }

  physicsB::finalize();
  physicsA::finalize();
  driver_finalize();
}
