// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * This file tests the ability to disable warnings about unknown [openmp]
 * pragmas on all supported compilers using the
 * AXOM_DISABLE_OMP_PRAGMA_WARNINGS build variable.
 * It should be compiled without the openmp flag.
 */

#include <iostream>

int main()
{
  const int SZ = 100;
  int arr[SZ];

#pragma omp parallel for
  for(int i = 0; i < SZ; ++i)
  {
    arr[i] = i;
  }

  std::cout << "Value of array element 0 is " << arr[0] << std::endl;

  return 0;
}
