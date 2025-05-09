// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// SPHINX_SLIC_BASIC_EXAMPLE_BEGIN

// SPHINX_SLIC_INCLUDES_BEGIN
// Slic includes
#include "axom/slic.hpp"
// SPHINX_SLIC_INCLUDES_END

using namespace axom;

//------------------------------------------------------------------------------
int main(int AXOM_UNUSED_PARAM(argc), char** AXOM_UNUSED_PARAM(argv))
{
  // SPHINX_SLIC_INIT_BEGIN

  slic::initialize();

  // SPHINX_SLIC_INIT_END

  slic::disableAbortOnError();

  // SPHINX_SLIC_FORMAT_MSG_BEGIN

  std::string format = std::string("<TIMESTAMP>\n") +
    std::string("[ <LEVEL> <TAG>]: <MESSAGE> \n") + std::string("FILE=<FILE>\n") +
    std::string("LINE=<LINE>\n\n");

  // SPHINX_SLIC_FORMAT_MSG_END

  // SPHINX_SLIC_SET_SEVERITY_BEGIN

  slic::setLoggingMsgLevel(slic::message::Debug);

  // SPHINX_SLIC_SET_SEVERITY_END

  // SPHINX_SLIC_SET_STREAM_BEGIN
  slic::addStreamToAllMsgLevels(new slic::GenericOutputStream(&std::cout, format));

  // SPHINX_SLIC_SET_STREAM_END

  // SPHINX_SLIC_SET_TAGGED_STREAM_BEGIN
  slic::addStreamToTag(new slic::GenericOutputStream(&std::cout, format), "myTag");

  // SPHINX_SLIC_SET_TAGGED_STREAM_END

  // SPHINX_SLIC_LOG_MESSAGES_BEGIN

  SLIC_DEBUG("Here is a debug message!");
  SLIC_INFO("Here is an info mesage!");
  SLIC_WARNING("Here is a warning!");
  SLIC_ERROR("Here is an error message!");
  SLIC_INFO_TAGGED("Here is a message for tagged streams with tag 'myTag'!", "myTag");

  // SPHINX_SLIC_LOG_MESSAGES_END

  // SPHINX_SLIC_FINALIZE_BEGIN

  slic::finalize();

  // SPHINX_SLIC_FINALIZE_END

  return 0;
}

// SPHINX_SLIC_BASIC_EXAMPLE_END
