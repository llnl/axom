// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

int runMIR_seq_hex(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result)
{
  return runMIR_hex<axom::SEQ_EXEC>(mesh, options, result);
}
