// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_HIP)
// Prototypes
int runMIR_hip_tri(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_hip_quad(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_hip_hex(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);

int runMIR_hip(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result)
{
  std::string shape = mesh["topologies/mesh/elements/shape"].as_string();
  int retval = 0;
  if(shape == "tri")
    retval = runMIR_hip_tri(mesh, options, result);
  else if(shape == "quad")
    retval = runMIR_hip_quad(mesh, options, result);
  else if(shape == "hex")
    retval = runMIR_hip_hex(mesh, options, result);
  return retval;
}
#else
int runMIR_hip(const conduit::Node &AXOM_UNUSED_PARAM(mesh),
               const conduit::Node &AXOM_UNUSED_PARAM(options),
               conduit::Node &AXOM_UNUSED_PARAM(result))
{
  return 0;
}
#endif
