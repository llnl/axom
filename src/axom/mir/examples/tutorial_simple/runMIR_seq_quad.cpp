// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

int runMIR_seq_quad(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result)
{
  return runMIR_quad<axom::SEQ_EXEC>(mesh, options, result);
}
