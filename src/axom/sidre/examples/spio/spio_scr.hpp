// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"

#include "mpi.h"

#ifdef AXOM_USE_SCR
  #include "scr.h"
#endif

#include "axom/sidre.hpp"
#include "conduit_relay.hpp"

using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::Group;
using axom::sidre::IOManager;
using axom::sidre::View;

#ifdef AXOM_USE_SCR
//------------------------------------------------------------------------------
TEST(spio_scr, spio_scr_writeread)
{
  SCR_Init();

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int num_output = num_ranks;

  /*
   * Create a DataStore and give it a small hierarchy of groups and views.
   *
   * The views are filled with repeatable nonsense data that will vary based
   * on rank.
   */
  DataStore* ds = new DataStore();

  Group* root = ds->getRoot();

  Group* flds = root->createGroup("fields");
  Group* flds2 = root->createGroup("fields2");

  Group* ga = flds->createGroup("a");
  Group* gb = flds2->createGroup("b");
  ga->createViewScalar<int>("i0", 101 * my_rank);
  gb->createView("i1")->allocate(DataType::c_int(10));
  int* i1_vals = gb->getView("i1")->getData();

  for(int i = 0; i < 10; i++)
  {
    i1_vals[i] = (i + 10) * (404 - my_rank - i);
  }

  SCR_Start_output("out_spio_parallel_write_read", SCR_FLAG_OUTPUT);

  /*
   * Contents of the DataStore written to files with IOManager.
   */
  int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD, true);

  writer.write(root, num_files, "out_spio_parallel_write_read", "sidre_hdf5");

  MPI_Barrier(MPI_COMM_WORLD);

  int valid = 1;
  int complete_rc = SCR_Complete_output(valid);
  EXPECT_EQ(complete_rc, SCR_SUCCESS);

  /*
   * SCR write operation is complete. Must finalize and initialize
   * before reading.
   */
  SCR_Finalize();
  SCR_Init();

  std::string root_name = "out_spio_parallel_write_read.root";

  /*
   * Create another DataStore that holds nothing but the root group.
   */
  DataStore* ds2 = new DataStore();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD, true);

  reader.read(ds2->getRoot(), root_name, false);

  /*
   * Verify that the contents of ds2 match those written from ds.
   */
  EXPECT_TRUE(ds2->getRoot()->isEquivalentTo(root));

  int testvalue = ds->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();
  int testvalue2 = ds2->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();

  EXPECT_EQ(testvalue, testvalue2);

  View* view_i1_orig = ds->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");
  View* view_i1_restored = ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");

  int num_elems = view_i1_orig->getNumElements();
  EXPECT_EQ(view_i1_restored->getNumElements(), num_elems);

  int* i1_orig = view_i1_orig->getData();
  int* i1_restored = view_i1_restored->getData();

  for(int i = 0; i < num_elems; ++i)
  {
    EXPECT_EQ(i1_orig[i], i1_restored[i]);
  }

  delete ds;
  delete ds2;

  SCR_Finalize();
}
#endif
