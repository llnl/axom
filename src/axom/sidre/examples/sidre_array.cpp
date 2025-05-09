// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core/Types.hpp"   // for Axom types
#include "axom/core/Macros.hpp"  // for Axom macros

#include "axom/sidre.hpp"  // for sidre
#include "axom/slic.hpp"   // for logging with slic

// MPI includes
#include <mpi.h>

// C/C++ includes
#include <iostream>

// aliases
namespace sidre = axom::sidre;

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm problem_comm = MPI_COMM_WORLD;

  int nranks = -1;
  int myrank = -1;
  MPI_Comm_rank(problem_comm, &myrank);
  MPI_Comm_size(problem_comm, &nranks);

  axom::slic::SimpleLogger logger;

  // STEP 0: create the data store
  sidre::DataStore* dataStore1 = new sidre::DataStore();
  sidre::Group* root1 = dataStore1->getRoot();

  // STEP 1: create an array with some data
  constexpr axom::IndexType NUM_NODES = 10;
  constexpr axom::IndexType DIMENSION = 4;
  sidre::MCArray<int> nodes_1(root1->createView("nodes_1/data"), NUM_NODES, DIMENSION);

  int value = 0;
  for(axom::IndexType i = 0; i < NUM_NODES; ++i)
  {
    for(axom::IndexType j = 0; j < DIMENSION; ++j)
    {
      nodes_1(i, j) = value;
      ++value;
    }  // END for all components
  }    // END for all nodes

  // DEBUG
  SLIC_INFO("Here is the array data in DataStore_1:\n");
  root1->print();
  std::cout << std::endl;
  // END DEBUG

  // STEP 2: save the array data in to a file
  sidre::IOManager sidre_io(problem_comm);
#if defined(AXOM_USE_HDF5)
  sidre_io.write(root1, nranks, "sidre_array_mesh", "sidre_hdf5");
#else
  sidre_io.write(root1, nranks, "sidre_array_mesh", "sidre_conduit_json");
#endif

  // STEP 3: read the data from the file into a new DataStore
  sidre::DataStore* dataStore2 = new sidre::DataStore();
  sidre::Group* root2 = dataStore2->getRoot();
  sidre_io.read(root2, "sidre_array_mesh.root");

  // DEBUG
  SLIC_INFO("Here is the array data in DataStore_2:\n");
  root2->print();
  std::cout << std::endl;
  // END DEBUG

  sidre::MCArray<int> nodes_2(root2->getView("nodes_1/data"));
  SLIC_ASSERT(nodes_2.shape()[0] == NUM_NODES);
  SLIC_ASSERT(nodes_2.shape()[1] == DIMENSION);

  // STEP 4: ensure the data is correct
  int expected_value = 0;
  for(axom::IndexType i = 0; i < NUM_NODES; ++i)
  {
    for(axom::IndexType j = 0; j < DIMENSION; ++j)
    {
      SLIC_ASSERT(nodes_2(i, j) == expected_value);
      ++expected_value;
    }  // END for all components
  }    // END for all nodes

  // STEP 5: delete the datastores
  delete dataStore2;
  dataStore2 = nullptr;

  delete dataStore1;
  dataStore1 = nullptr;

  MPI_Finalize();
  return 0;
}
