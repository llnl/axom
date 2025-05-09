// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_OPENMP)
int runMIR_omp_quad(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result)
{
  return runMIR_quad<axom::OMP_EXEC>(mesh, options, result);
}
#else
int runMIR_omp_quad(const conduit::Node &AXOM_UNUSED_PARAM(mesh),
                    const conduit::Node &AXOM_UNUSED_PARAM(options),
                    conduit::Node &AXOM_UNUSED_PARAM(result))
{
  return 0;
}
#endif
