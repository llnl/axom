// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"

// MPI includes
#include <mpi.h>

// C/C++ includes
#include <iostream>

// aliases
namespace sidre = axom::sidre;

int main(int argc, char** argv)
{
  //
  // Basic MPI setup with global communicator
  //
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);
  const int nranks = mpi_raii_wrapper.num_ranks();

  axom::slic::SimpleLogger logger;

  // STEP 0: create a data store instance, with one child group of root
  sidre::DataStore* dataStore1 = new sidre::DataStore();
  sidre::Group* root1 = dataStore1->getRoot();

  // STEP 1: Create a sidre::MCArray with data held in a Sidre View.
  //         A sidre::MCArray is similar to axom::Array, except that
  //         the array data can be accessed using a Sidre View.
  constexpr axom::IndexType NUM_NODES = 10;
  constexpr axom::IndexType DIMENSION = 4;
  sidre::MCArray<int> nodes_1(root1->createView("nodes_1/data"), NUM_NODES, DIMENSION);

  // STEP 1a: Initialize array data values
  int value = 0;
  for(axom::IndexType i = 0; i < NUM_NODES; ++i)
  {
    for(axom::IndexType j = 0; j < DIMENSION; ++j)
    {
      nodes_1(i, j) = value;
      ++value;
    }
  }

  // DEBUG
  // STEP 1b: Print contents of data store showing that array data is in View
  SLIC_INFO("Here is the array data in DataStore_1:\n");
  root1->print();
  std::cout << std::endl;
  // END DEBUG

  // STEP 2: Write out the data store contents to a file.
  sidre::IOManager sidre_io(MPI_COMM_WORLD);
  const auto protocol = sidre::Group::getDefaultIOProtocol();
  sidre_io.write(root1, nranks, "sidre_array_mesh", protocol);

  // STEP 3: Create a new data store object and read the data from the file into it.
  sidre::DataStore* dataStore2 = new sidre::DataStore();
  sidre::Group* root2 = dataStore2->getRoot();
  sidre_io.read(root2, "sidre_array_mesh.root");

  // DEBUG
  // STEP 3a: Print contents of data store showing that array data was read into a View
  SLIC_INFO("Here is the array data in DataStore_2:\n");
  root2->print();
  std::cout << std::endl;
  // END DEBUG

  // STEP 3b: Create a new sidre::MCArray with data held in a View in the
  //          new data store. Verify that the shape and size of the data is
  //          what it is expected to be.
  sidre::MCArray<int> nodes_2(root2->getView("nodes_1/data"));
  SLIC_ASSERT(nodes_2.shape()[0] == NUM_NODES);
  SLIC_ASSERT(nodes_2.shape()[1] == DIMENSION);

  // STEP 4: Check the values are correct when accessed through the sidre::MCArray
  [[maybe_unused]] int expected_value = 0;
  for(axom::IndexType i = 0; i < NUM_NODES; ++i)
  {
    for(axom::IndexType j = 0; j < DIMENSION; ++j)
    {
      SLIC_ASSERT(nodes_2(i, j) == expected_value);
      ++expected_value;
    }
  }

  // STEP 5: Delete both data stores (and associated data).
  delete dataStore2;
  dataStore2 = nullptr;

  delete dataStore1;
  dataStore1 = nullptr;

  return 0;
}
