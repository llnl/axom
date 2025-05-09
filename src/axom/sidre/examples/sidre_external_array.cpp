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

#if defined(AXOM_USE_HDF5)
//------------------------------------------------------------------------------
void sidre_write(MPI_Comm comm,
                 const std::string& file,
                 int* data,
                 axom::IndexType numTuples,
                 axom::IndexType numComponents)
{
  SLIC_ASSERT(comm != MPI_COMM_NULL);
  SLIC_ASSERT(data != nullptr);
  SLIC_ASSERT(!file.empty());

  int nranks = -1;
  MPI_Comm_size(comm, &nranks);

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  sidre::View* view = root->createView("data");

  sidre::IndexType shape[2];
  shape[0] = numTuples;
  shape[1] = numComponents;

  view->setExternalDataPtr(sidre::INT32_ID, 2, shape, data);

  // DEBUG
  SLIC_INFO("Here is the data that is begin dumped:");
  root->print();
  std::cout << std::endl;
  // DEBUG

  // STEP 2: save the array data in to a file
  sidre::IOManager sidre_io(comm);
  sidre_io.write(root, nranks, file, "sidre_hdf5");
}

//------------------------------------------------------------------------------
void sidre_read(MPI_Comm comm,
                const std::string& file,
                int*& data,
                axom::IndexType& numTuples,
                axom::IndexType& numComponents)
{
  SLIC_ASSERT(comm != MPI_COMM_NULL);
  SLIC_ASSERT(data == nullptr);
  SLIC_ASSERT(!file.empty());

  int nranks = -1;
  MPI_Comm_size(comm, &nranks);

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  sidre::IOManager sidre_io(comm);
  sidre_io.read(root, file);

  SLIC_ASSERT(root->hasChildView("data"));
  sidre::View* view = root->getView("data");
  SLIC_ASSERT(view->isDescribed());
  SLIC_ASSERT(view->isExternal());

  // get the array shape information
  sidre::IndexType shape[2];
  view->getShape(2, shape);
  numTuples = shape[0];
  numComponents = shape[1];

  axom::IndexType nelems = view->getNumElements();
  SLIC_ASSERT(nelems == (numTuples * numComponents));

  // allocate external data
  data = axom::allocate<int>(nelems);
  SLIC_ASSERT(data != nullptr);

  // set external data for the view
  view->setExternalDataPtr(data);

  // load the external data
  sidre_io.loadExternalData(root, file);

  // DEBUG
  SLIC_INFO("Here is the data that was read back:");
  root->print();
  std::cout << std::endl;
  // DEBUG
}
#endif
//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
#if defined(AXOM_USE_HDF5)
  MPI_Comm problem_comm = MPI_COMM_WORLD;

  axom::slic::SimpleLogger logger;

  // STEP 0: create some data
  constexpr axom::IndexType NUM_NODES = 10;
  constexpr axom::IndexType DIMENSION = 4;
  constexpr axom::IndexType NSIZE = NUM_NODES * DIMENSION;

  int* data = axom::allocate<int>(NSIZE);
  SLIC_ASSERT(data != nullptr);

  for(int i = 0; i < NSIZE; ++i)
  {
    data[i] = (i + 1) * 10;
  }

  // STEP 1: dump the data to a file using sidre
  SLIC_INFO("Writting data...");
  sidre_write(problem_comm, "sidre_external_array_mesh", data, NUM_NODES, DIMENSION);
  SLIC_INFO("[DONE]");

  // STEP 2: read the data from a file using sidre
  int* data2 = nullptr;
  axom::IndexType ntuples = -1;
  axom::IndexType ncomp = -1;

  SLIC_INFO("Reading data...");
  sidre_read(problem_comm, "sidre_external_array_mesh.root", data2, ntuples, ncomp);
  SLIC_INFO("[DONE]");

  // STEP 3: check the data
  SLIC_ASSERT(ntuples == NUM_NODES);
  SLIC_ASSERT(ncomp == DIMENSION);
  for(int i = 0; i < NSIZE; ++i)
  {
    SLIC_ASSERT(data[i] == data2[i]);
  }

  // STEP 4: deallocate
  axom::deallocate(data);
  axom::deallocate(data2);
#endif

  MPI_Finalize();
  return 0;
}
