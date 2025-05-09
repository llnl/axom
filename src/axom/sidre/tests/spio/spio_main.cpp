// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "mpi.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

namespace
{
const std::string PROTOCOL = axom::sidre::Group::getDefaultIOProtocol();
const std::string ROOT_EXT = ".root";
}  // namespace

#include "spio_basic.hpp"
#include "spio_parallel.hpp"
#include "spio_serial.hpp"

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}
