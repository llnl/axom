// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "MIRApplication.hpp"

int main(int argc, char **argv)
{
  MIRApplication app;
  int retval = app.initialize(argc, argv);
  if(retval == 0)
  {
    retval = app.execute();
  }

  return retval;
}
