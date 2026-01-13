// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**************************************************************************
 *************************************************************************/

#include "mpi.h"

#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::Group;
using axom::sidre::IOManager;
using namespace axom::utilities;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  axom::slic::SimpleLogger logger;

  SLIC_ERROR_IF(
    argc != 3,
    "Missing command line arguments. \n\t" << "Usage: spio_IOWrite <num_files> <base_file_name>");

  size_t num_files = 0;
  std::string file_base;
  if(argc == 3)
  {
    num_files = static_cast<size_t>(atoi(argv[1]));
    file_base = argv[2];
  }
  else
  {
    return 0;
  }

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  DataStore* ds = new DataStore();
  SLIC_ASSERT(ds);
  Group* root = ds->getRoot();

  Group* flds = root->createGroup("fields");
  Group* flds2 = root->createGroup("fields2");

  Group* ga = flds->createGroup("a");
  Group* gb = flds2->createGroup("b");
  ga->createViewScalar<int>("i0", my_rank + 101);
  gb->createViewScalar<int>("i1", 4 * my_rank * my_rank + 404);

  if(my_rank == 0)
  {
    std::string dir;
    filesystem::getDirName(dir, file_base);
    filesystem::makeDirsForPath(dir);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  IOManager writer(MPI_COMM_WORLD);
  writer.write(root, num_files, file_base, "sidre_hdf5");

  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank == 0)
  {
    Group* extra = root->createGroup("extra");
    extra->createViewScalar<double>("dval", 1.1);
    Group* child = extra->createGroup("child");
    child->createViewScalar<int>("ival", 7);
    child->createViewString("word0", "hello");
    child->createViewString("word1", "world");

    std::string root_name = file_base + ".root";
    writer.writeGroupToRootFile(extra, root_name);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  delete ds;

  MPI_Finalize();

  return 0;
}
