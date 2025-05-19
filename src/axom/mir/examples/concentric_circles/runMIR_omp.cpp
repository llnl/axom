// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_OPENMP)
int runMIR_omp(int dimension,
               const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result)
{
  int retval = 0;
  if(dimension == 3)
  {
    retval = runMIR<axom::OMP_EXEC, 3>(mesh, options, result);
  }
  else
  {
    retval = runMIR<axom::OMP_EXEC, 2>(mesh, options, result);
  }
  return retval;
}
#else
int runMIR_omp(int AXOM_UNUSED_PARAM(dimension),
               const conduit::Node &AXOM_UNUSED_PARAM(mesh),
               const conduit::Node &AXOM_UNUSED_PARAM(options),
               conduit::Node &AXOM_UNUSED_PARAM(result))
{
  return 0;
}
#endif
