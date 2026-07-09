// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/utilities/Abort.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

#include <cstdlib>

namespace axom::utilities
{
[[noreturn]] void processAbort()
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
}  // end namespace axom::utilities
