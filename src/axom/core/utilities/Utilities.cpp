// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file    Utilities.cpp
 *
 * \brief   Implementation file for utility functions.
 *
 */

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include <cstdlib>  // for exit, EXIT_SUCCESS, EXIT_FAILURE

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

namespace axom
{
namespace utilities
{
void processAbort()
{
#ifndef AXOM_USE_MPI
  abort();
#else
  int mpi = 0;
  MPI_Initialized(&mpi);
  if(mpi)
  {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  abort();
#endif
}

int binomialCoefficient(int n, int k)
{
  if(k > n || k < 0)  // check if out-of-bounds
  {
    return 0;
  }
  if(k == n || k == 0)  // early return
  {
    return 1;
  }
  if(k > n - k)  // exploit symmetry to reduce work
  {
    k = n - k;
  }

  int val = 1;
  for(int i = 1; i <= k; ++i)
  {
    val *= (n - k + i);
    val /= i;
  }
  return val;
}

}  // end namespace utilities
}  // end namespace axom
