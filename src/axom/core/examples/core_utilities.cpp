// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_utilities.cpp
 *  \brief This example code is a demonstration of the Axom Core utilites.
 */

/* This example code contains snippets used in the Core Sphinx documentation.
  * They begin and end with comments such as
  *
  * timer_start
  * timer_end
  *
  * each prepended with an underscore.
  */

#ifdef WIN32
  #include "windows.h"
void sleep(int numSeconds)
{
  int numMilliSecs = numSeconds * 1000;
  Sleep(numMilliSecs);
}
#else
  #include <unistd.h>  // for sleep()
#endif

// C/C++ includes
#include <iostream>
#include <vector>

// Axom includes
// _header_start
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/core/utilities/Timer.hpp"
#include "axom/core/utilities/About.hpp"
// _header_end

// _fs_string_start
void demoFileSystemAndString(const char* argv0)
{
  using namespace axom::utilities;

  // Get the current directory
  std::string cwd = filesystem::getCWD();

  // Split it on file separator.
#if WIN32
  const char pathsep = '\\';
#else
  const char pathsep = '/';
#endif
  std::vector<std::string> cmp = string::split(cwd, pathsep);

  // Count how many start with "ax" or end with "exe"
  // (we could also use std::count_if)
  int matchcount = 0;
  std::string prefix {"ax"};
  std::string suffix {"exe"};
  const int N = static_cast<int>(cmp.size());
  for(int i = 0; i < N; ++i)
  {
    if(string::startsWith(cmp[i], prefix) || string::endsWith(cmp[i], suffix))
    {
      matchcount += 1;
    }
  }
  std::cout << "Found " << matchcount << " path components starting with " << prefix
            << " or ending with " << suffix << "." << std::endl;

  // Append "hello.txt"
  std::string hellofile = filesystem::joinPath(cwd, "hello.txt", std::string {pathsep});

  // Does this exist?
  std::cout << "The file \"hello.txt\" ";
  if(filesystem::pathExists(hellofile))
  {
    std::cout << "exists ";
  }
  else
  {
    std::cout << "DOES NOT exist ";
  }
  std::cout << "in the current working directory." << std::endl;

  // Does argv0 exist?
  if(filesystem::pathExists(argv0))
  {
    std::cout << argv0 << " exists ";
  }
  else
  {
    std::cout << argv0 << " DOES NOT exist ";
  }
  std::cout << "in the current working directory." << std::endl;

  // sleep for a second
  sleep(1);
}
// _fs_string_end

int main(int argc, char** argv)
{
  // _about_start
  std::cout << "Here is a message telling you about Axom." << std::endl;
  axom::about();

  std::cout << "The version string '" << axom::getVersion() << "' is part of the previous message, "
            << std::endl
            << " and is also available separately." << std::endl;
  // _about_end

  // _timer_start
  axom::utilities::Timer t;

  t.start();

  if(argc == 1)
  {
    std::cerr << "Error: no path given on command line" << std::endl;
    return 1;
  }
  else
  {
    demoFileSystemAndString(argv[0]);
  }

  t.stop();

  std::cout << "The tests took " << t.elapsedTimeInMilliSec() << " ms." << std::endl;
  // _timer_end

  return 0;
}
