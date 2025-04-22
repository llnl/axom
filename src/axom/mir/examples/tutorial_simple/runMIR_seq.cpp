// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

// Prototypes
int runMIR_seq_tri(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_seq_quad(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);
int runMIR_seq_hex(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result);

int runMIR_seq(const conduit::Node &mesh, const conduit::Node &options, conduit::Node &result)
{
  std::string shape = mesh["topologies/mesh/elements/shape"].as_string();
  int retval = 0;
  if(shape == "tri")
    retval = runMIR_seq_tri(mesh, options, result);
  else if(shape == "quad")
    retval = runMIR_seq_quad(mesh, options, result);
  else if(shape == "hex")
    retval = runMIR_seq_hex(mesh, options, result);
  return retval;
}
