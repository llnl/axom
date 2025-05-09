// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// C/C++ includes
#include <iostream>  // For std::cerr

// Logging includes
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

using namespace axom;

void customAbortFunction()
{
  // finalize logging if needed
  if(slic::isInitialized())
  {
    slic::finalize();
  }
  std::cerr << "Custom abort function called!\n";
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

  //----------------------------------------------------------------------------
  // STEP 1: Create log streams
  //----------------------------------------------------------------------------

  // setup log stream for ALL messages, including FATAL, ERROR and WARNING
  std::string console_format = std::string("[<LEVEL>]: <MESSAGE>\n");
  slic::LogStream* console = new slic::GenericOutputStream(&std::cerr, console_format);

  //----------------------------------------------------------------------------
  // STEP 2: add streams to logger
  //----------------------------------------------------------------------------

  slic::addStreamToAllMsgLevels(console);

  //----------------------------------------------------------------------------
  // STEP 3: Register a custom abort function
  //----------------------------------------------------------------------------
  slic::setAbortFunction(customAbortFunction);

  //----------------------------------------------------------------------------
  // STEP 4: Trigger the abort function by raising an error
  //----------------------------------------------------------------------------
  SLIC_ERROR("Here is an error message!");

  return 0;
}
