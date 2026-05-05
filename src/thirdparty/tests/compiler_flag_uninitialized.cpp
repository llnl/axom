// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * This file tests the ability to disable warnings about uninitialized variables
 * on all supported compilers using the AXOM_DISABLE_UNINITIALIZED_WARNINGS
 * build variable.
 */

#include <iostream>
#include <cstdlib>

int main()
{
  int* result;  // Note: variable not allocated or initialized

  if(rand() % 2 == 0)
  {
    *result = 5;
  }

  return 0;
}
