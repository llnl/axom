// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: example.cpp
///
//-----------------------------------------------------------------------------

#include "axom/core.hpp"
#include "axom/fmt.hpp"

#include <iostream>

int main()
{
  // Using fmt library exported by axom
  std::cout << axom::fmt::format("Example of using an installed version of Axom {}",
                                 axom::getVersion())
            << std::endl
            << std::endl;

  // Uses installed axom library
  axom::about();

  return 0;
}
