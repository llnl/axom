// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_HIP)
int runMIR_hip(int dimension,
               const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result)
{
  constexpr int HIP_BLOCK_SIZE = 64;
  using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
  int retval = 0;
  if(dimension == 3)
  {
    retval = runMIR<hip_exec, 3>(mesh, options, result);
  }
  else
  {
    retval = runMIR<hip_exec, 2>(mesh, options, result);
  }
  return retval;
}
#else
int runMIR_hip(int AXOM_UNUSED_PARAM(dimension),
               const conduit::Node &AXOM_UNUSED_PARAM(mesh),
               const conduit::Node &AXOM_UNUSED_PARAM(options),
               conduit::Node &AXOM_UNUSED_PARAM(result))
{
  return 0;
}
#endif
