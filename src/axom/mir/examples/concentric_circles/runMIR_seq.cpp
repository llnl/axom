// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

int runMIR_seq(int dimension,
               const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result)
{
  int retval = 0;
  if(dimension == 3)
  {
    retval = runMIR<axom::SEQ_EXEC, 3>(mesh, options, result);
  }
  else
  {
    retval = runMIR<axom::SEQ_EXEC, 2>(mesh, options, result);
  }
  return retval;
}
