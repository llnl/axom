// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: check_axom_configuration.cpp
///
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/fmt.hpp"

#include <iostream>

int main()
{
  std::cout << "Checking properties of installed Axom library. \n";

  // Check Axom version
  std::cout << axom::fmt::format("Version: {}", axom::getVersion()) << "\n\n";

  // Print Axom about
  axom::about();

  return 0;
}
