// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "HMApplication.hpp"

int main(int argc, char **argv)
{
  HMApplication app;
  int retval = app.initialize(argc, argv);
  if(retval == 0)
  {
    retval = app.execute();
  }

  return retval;
}
