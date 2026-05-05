// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * This file tests the ability to disable warnings about unused parameters
 * on all supported compilers using the AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS
 * build variable.
 */

#include <iostream>

void foo(int param) { std::cout << "Hello " << std::endl; }

int main()
{
  foo(2);
  return 0;
}
