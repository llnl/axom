// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * This file tests the ability to disable warnings about unused local typedefs
 * on all supported compilers using the AXOM_DISABLE_UNUSED_LOCAL_TYPEDEF
 * build variable.
 */

#include <iostream>

int main()
{
  // NOLINTBEGIN
  typedef int IntT;
  // NOLINTEND

  std::cout << "I have defined type IntT, but am not using it." << std::endl;
  return 0;
}
