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

namespace axom::utilities
{

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

}  // end namespace axom::utilities
