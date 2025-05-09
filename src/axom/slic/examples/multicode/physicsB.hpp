// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PHYSICSB_HPP_
#define PHYSICSB_HPP_

// SLIC includes
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

using namespace axom;

namespace physicsB
{
std::ofstream physicsB_log;

inline int getRandInt(const int start, const int end)
{
  return (std::rand() % (end - start) + start);
}

inline slic::message::Level getRandomLevel()
{
  return (static_cast<slic::message::Level>(getRandInt(0, slic::message::Num_Levels)));
}
//------------------------------------------------------------------------------
void init()
{
  std::string current_logger = slic::getActiveLoggerName();

  std::string physicsB_format =
    std::string("***************************************************\n") +
    std::string("<TIMESTAMP>\n") +
    std::string("***************************************************\n") +
    std::string("[<LEVEL>]: <MESSAGE>\n") + std::string("\t FILE:<FILE>\n") +
    std::string("\t LINE:<LINE>\n");

  physicsB_log.open("physicsB.log");
  slic::LogStream* ls = new slic::GenericOutputStream(&physicsB_log, physicsB_format);

  slic::createLogger("physicsB", slic::inherit::errors_and_warnings);
  slic::activateLogger("physicsB");
  slic::disableAbortOnError();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::addStreamToAllMsgLevels(ls);

  slic::activateLogger(current_logger);
}

//------------------------------------------------------------------------------
void timestep(int step, int n)
{
  std::string current_logger = slic::getActiveLoggerName();
  slic::activateLogger("physicsB");

  std::ostringstream oss;
  oss << "n=" << n << " physicsB cycles";
  slic::logMessage(slic::message::Info, oss.str(), __FILE__, __LINE__);

  for(int i = 0; i < n; ++i)
  {
    slic::message::Level random = getRandomLevel();

    oss.str("");
    oss << "cycle=" << step << " subcycle=" << i << " a random message!";
    slic::logMessage(random, oss.str(), __FILE__, __LINE__);
  }

  slic::activateLogger(current_logger);
}

//------------------------------------------------------------------------------
inline void finalize() { physicsB_log.close(); }

} /* namespace physicsB */

#endif /* PHYSICSB_HPP_ */
