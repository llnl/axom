// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* An excerpt from this test file is used in the Sidre Sphinx documentation,
 * denoted by the comment strings
 *
 * parallel_io_headers_start
 * parallel_io_headers_end
 * parallel_io_save_start
 * parallel_io_save_end
 * parallel_io_load_start
 * parallel_io_load_end
 *
 * prepended with an underscore.
 */

#include "gtest/gtest.h"

// _parallel_io_headers_start
#include "axom/config.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"
#ifdef AXOM_USE_HDF5
  #include "conduit_relay_io_hdf5.hpp"
#endif

#include "axom/sidre.hpp"
#include "axom/fmt.hpp"

#include "mpi.h"
// _parallel_io_headers_end

using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::Group;
using axom::sidre::IOManager;
using axom::sidre::View;

namespace
{
/** Returns the number of output files for spio  */
int numOutputFiles(int numRanks)
{
#ifdef AXOM_USE_HDF5
  return std::max(numRanks / 2, 1);
#else
  return numRanks;
#endif
}

}  // namespace

//------------------------------------------------------------------------------

TEST(spio_parallel, parallel_writeread)
{
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  const int num_output = numOutputFiles(num_ranks);

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

  // Add a scalar view under group "a"
  Group* ga = flds->createGroup("a");
  ga->createViewScalar<int>("i0", 101 * my_rank);

  // Add an array view under group "b"
  Group* gb = flds2->createGroup("b");
  gb->createView("i1")->allocate(DataType::c_int(10));
  int* i1_vals = gb->getView("i1")->getData();

  for(int i = 0; i < 10; i++)
  {
    i1_vals[i] = (i + 10) * (404 - my_rank - i);
  }

  // Add an array view of size 0 under group "c"
  Group* gc = flds2->createGroup("c");
  gc->createView("i2")->allocate(DataType::c_int(0));

  // The namespace qualifying IOManager is not necessary for compilation
  // (because are already using axom::sidre::IOManager), but we've kept it
  // because the following code is used in the documentation as an example.
  // _parallel_io_save_start
  /*
   * Contents of the DataStore written to files with IOManager.
   */
  int num_files = num_output;
  axom::sidre::IOManager writer(MPI_COMM_WORLD);

  const std::string file_name = "out_spio_parallel_write_read";

  writer.write(root, num_files, file_name, PROTOCOL, "tree_%07d");

  std::string root_name = file_name + ROOT_EXT;
  // _parallel_io_save_end

  /*
   * Extra stuff to exercise writeGroupToRootFile
   * Note: This is only valid for the 'sidre_hdf5' protocol
   */
  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank == 0 && PROTOCOL == "sidre_hdf5")
  {
    DataStore* dsextra = new DataStore();
    Group* extra = dsextra->getRoot()->createGroup("extra");
    extra->createViewScalar<double>("dval", 1.1);
    Group* child = extra->createGroup("child");
    child->createViewScalar<int>("ival", 7);
    child->createViewString("word0", "hello");
    child->createViewString("word1", "world");

    writer.writeGroupToRootFile(extra, root_name);

    Group* path_test = dsextra->getRoot()->createGroup("path_test");

    path_test->createViewScalar<int>("path_val", 9);
    path_test->createViewString("word2", "again");

    writer.writeGroupToRootFileAtPath(path_test, root_name, "extra/child");

    View* view_test = dsextra->getRoot()->createViewString("word3", "new_view");

    writer.writeViewToRootFileAtPath(view_test, root_name, "extra/child/path_test");

    delete dsextra;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /*
   * Read the root file on rank 1, unless this is a serial run.
   */
  if((my_rank == 1 || num_ranks == 1) && PROTOCOL == "sidre_hdf5")
  {
    conduit::Node root_node;
    conduit::relay::io::load(root_name, "hdf5", root_node);

    std::string pattern = root_node["tree_pattern"].as_string();

    const conduit::Node& n = root_node["extra"];

    EXPECT_EQ(pattern, "tree_%07d");

    double dval = n["dval"].to_double();

    EXPECT_TRUE(dval > 1.0000009 && dval < 1.1000001);

    EXPECT_EQ(n["child"]["ival"].to_int(), 7);
    EXPECT_EQ(n["child"]["word0"].as_string(), "hello");
    EXPECT_EQ(n["child"]["word1"].as_string(), "world");
    EXPECT_EQ(n["child"]["path_test"]["path_val"].to_int(), 9);
    EXPECT_EQ(n["child"]["path_test"]["word2"].as_string(), "again");
    EXPECT_EQ(n["child"]["path_test"]["word3"].as_string(), "new_view");
  }

  // _parallel_io_load_start
  /*
   * Create another DataStore that holds nothing but the root group.
   */
  DataStore* ds2 = new DataStore();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), root_name, PROTOCOL);
  // _parallel_io_load_end

  /*
   * Verify that the contents of ds2 match those written from ds.
   */
  EXPECT_TRUE(ds2->getRoot()->isEquivalentTo(root));

  // check group "a", which should contain a view with a single scalar
  {
    int testvalue = ds->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();
    int testvalue2 = ds2->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();

    EXPECT_EQ(testvalue, testvalue2);
  }

  // check group "b", which should contain a view with an array of ints
  {
    View* view_i1_orig = ds->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");
    View* view_i1_restored = ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");

    int num_elems = view_i1_orig->getNumElements();
    EXPECT_EQ(view_i1_restored->getNumElements(), num_elems);
    if(view_i1_restored->getNumElements() == num_elems)
    {
      int* i1_orig = view_i1_orig->getData();
      int* i1_restored = view_i1_restored->getData();

      for(int i = 0; i < num_elems; ++i)
      {
        EXPECT_EQ(i1_orig[i], i1_restored[i]);
      }
    }
  }

  // check group "c", which should contain a view with an array of size 0
  {
    View* view_i2_orig = ds->getRoot()->getGroup("fields2")->getGroup("c")->getView("i2");
    View* view_i2_restored = ds2->getRoot()->getGroup("fields2")->getGroup("c")->getView("i2");

    int num_elems_orig = view_i2_orig->getNumElements();
    int num_elems_rest = view_i2_restored->getNumElements();
    EXPECT_EQ(0, num_elems_orig);
    EXPECT_EQ(0, num_elems_rest);
  }

  delete ds;
  delete ds2;
}

//----------------------------------------------------------------------
TEST(spio_parallel, write_read_write)
{
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  const int num_files = numOutputFiles(num_ranks);
  const std::string filename = axom::fmt::format("out_spio_WRW_{}", num_ranks);

  // Initialize a datastore and dump to disk
  DataStore* ds = new DataStore();
  ds->getRoot()->createViewScalar("grp/i", 2);
  ds->getRoot()->createViewScalar("grp/f", 3.0);
  IOManager writer_a(MPI_COMM_WORLD);
  writer_a.write(ds->getRoot(), num_files, filename, PROTOCOL);

  // Create another DataStore to read into.
  DataStore ds_r;
  IOManager reader(MPI_COMM_WORLD);
  reader.read(ds_r.getRoot(), filename + ROOT_EXT);

  // Dump this datastore to disk.
  // Regression for sidre_hdf5 protocol:
  // This used to produce the following HDF5 error:
  //  HDF5-DIAG: Error detected in HDF5 (1.8.16) thread 0:
  //    #000: H5F.c line 522 in H5Fcreate(): unable to create file
  //      major: File accessibility
  //      minor: Unable to open file
  //    #001: H5Fint.c line 1024 in H5F_open(): unable to truncate a file which
  // is already open
  //      major: File accessibility
  //      minor: Unable to open file
  IOManager writer_b(MPI_COMM_WORLD);
  writer_b.write(ds_r.getRoot(), num_files, filename, PROTOCOL);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(spio_parallel, external_writeread)
{
  if(PROTOCOL != "sidre_hdf5")
  {
    SUCCEED() << "Loading external data in spio only currently supported "
              << " for 'sidre_hdf5' protocol";
    return;
  }

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  const int num_output = numOutputFiles(num_ranks);

  const int nvals = 10;
  int orig_vals1[nvals], orig_vals2[nvals];
  for(int i = 0; i < 10; i++)
  {
    orig_vals1[i] = (i + 10) * (404 - my_rank - i);
    orig_vals2[i] = (i + 10) * (404 - my_rank - i) + 20;
  }

  /*
   * Create a DataStore and give it a small hierarchy of groups and views.
   *
   * The views are filled with repeatable nonsense data that will vary based
   * on rank.
   */
  DataStore* ds1 = new DataStore();

  Group* root1 = ds1->getRoot();

  Group* flds = root1->createGroup("fields");
  Group* flds2 = root1->createGroup("fields2");

  Group* ga = flds->createGroup("a");
  Group* gb = flds2->createGroup("b");
  ga->createView("external_array", axom::sidre::INT_ID, nvals, orig_vals1);
  gb->createView("external_undescribed")->setExternalDataPtr(orig_vals2);

  /*
   * Contents of the DataStore written to files with IOManager.
   */
  const int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD);

  const std::string file_name = "out_spio_external_write_read";

  writer.write(root1, num_files, file_name, PROTOCOL);

  /*
   * Create another DataStore than holds nothing but the root group.
   */
  DataStore* ds2 = new DataStore();
  Group* root2 = ds2->getRoot();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD);

  reader.read(root2, file_name + ROOT_EXT);

  int restored_vals1[nvals], restored_vals2[nvals];
  for(int i = 0; i < nvals; ++i)
  {
    restored_vals1[i] = -1;
    restored_vals2[i] = -1;
  }

  View* view1 = root2->getView("fields/a/external_array");
  view1->setExternalDataPtr(restored_vals1);

  View* view2 = root2->getView("fields2/b/external_undescribed");
  view2->setExternalDataPtr(restored_vals2);

  reader.loadExternalData(root2, file_name + ROOT_EXT);

  enum SpioTestResult
  {
    SPIO_TEST_SUCCESS = 0,
    HIERARCHY_ERROR = 1 << 0,
    EXT_ARRAY_ERROR = 1 << 1,
    EXT_UNDESC_ERROR = 1 << 2
  };
  int result = SPIO_TEST_SUCCESS;

  /*
   * Verify that the contents of ds2 match those written from ds.
   */
  EXPECT_TRUE(ds2->getRoot()->isEquivalentTo(root1));
  if(!ds2->getRoot()->isEquivalentTo(root1))
  {
    result |= HIERARCHY_ERROR;
  }
  SLIC_WARNING_IF(result & HIERARCHY_ERROR, "Tree layouts don't match");

  EXPECT_EQ(view1->getNumElements(), nvals);
  if(view1->getNumElements() != nvals)
  {
    result |= EXT_ARRAY_ERROR;
  }
  else
  {
    for(int i = 0; i < nvals; ++i)
    {
      EXPECT_EQ(orig_vals1[i], restored_vals1[i]);
      if(orig_vals1[i] != restored_vals1[i])
      {
        result |= EXT_ARRAY_ERROR;
        break;
      }
    }
  }
  SLIC_WARNING_IF(result & EXT_ARRAY_ERROR, "External_array was not correctly loaded");

  /*
   * external_undescribed was not written to disk (since it is undescribed)
   * make sure it was not read in.
   */
  for(int i = 0; i < nvals; ++i)
  {
    EXPECT_EQ(-1, restored_vals2[i]);
    if(-1 != restored_vals2[i])
    {
      result |= EXT_UNDESC_ERROR;
      break;
    }
  }
  SLIC_WARNING_IF(result & EXT_UNDESC_ERROR, "External_undescribed data was modified.");

  delete ds1;
  delete ds2;
}

//------------------------------------------------------------------------------
TEST(spio_parallel, external_piecemeal_writeread)
{
  if(PROTOCOL != "sidre_hdf5")
  {
    SUCCEED() << "Loading external data in spio only currently supported "
              << " for 'sidre_hdf5' protocol";
    return;
  }

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  const int num_output = numOutputFiles(num_ranks);

  const int nvals = 10;
  int orig_vals1[nvals], orig_vals2[nvals];
  for(int i = 0; i < 10; i++)
  {
    orig_vals1[i] = (i + 10) * (404 - my_rank - i);
    orig_vals2[i] = (i + 10) * (404 - my_rank - i) + 20;
  }

  /*
   * Create a DataStore and give it a small hierarchy of groups and views.
   *
   * The views are filled with repeatable nonsense data that will vary based
   * on rank.
   */
  DataStore* ds1 = new DataStore();

  Group* root1 = ds1->getRoot();

  Group* flds = root1->createGroup("testdata/fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  ga->createView("vals1", axom::sidre::INT_ID, nvals, orig_vals1);
  gb->createView("vals2", axom::sidre::INT_ID, nvals, orig_vals2);

  /*
   * Contents of the DataStore written to files with IOManager.
   */
  const int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD);

  std::string file_name = "out_spio_external_piecemeal_write_read";

  writer.write(root1, num_files, file_name, PROTOCOL);

  /*
   * Create another DataStore to be the destination for input from file
   */
  DataStore* ds2 = new DataStore();
  Group* root2 = ds2->getRoot();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD);

  reader.read(root2, file_name + ROOT_EXT);

  // pollute values so we know they are changed
  int restored_vals1[nvals], restored_vals2[nvals];
  for(int i = 0; i < nvals; ++i)
  {
    restored_vals1[i] = -1;
    restored_vals2[i] = -1;
  }

  Group* group_a = root2->getGroup("testdata/fields/a");
  Group* group_b = root2->getGroup("testdata/fields/b");

  View* view1 = group_a->getView("vals1");
  view1->setExternalDataPtr(restored_vals1);

  View* view2 = group_b->getView("vals2");
  view2->setExternalDataPtr(restored_vals2);

  reader.loadExternalData(root2, group_a, file_name + ROOT_EXT);

  enum SpioTestResult
  {
    SPIO_TEST_SUCCESS = 0,
    HIERARCHY_ERROR = 1 << 0,
    EXT_ARRAY_ERROR = 1 << 1,
    EXT_UNDESC_ERROR = 1 << 2
  };
  int result = SPIO_TEST_SUCCESS;

  /*
   * Verify that the contents of view1 are restored.
   */
  EXPECT_EQ(view1->getNumElements(), nvals);
  if(view1->getNumElements() != nvals)
  {
    result |= EXT_ARRAY_ERROR;
  }
  else
  {
    for(int i = 0; i < nvals; ++i)
    {
      EXPECT_EQ(orig_vals1[i], restored_vals1[i]);
      if(orig_vals1[i] != restored_vals1[i])
      {
        result |= EXT_ARRAY_ERROR;
        break;
      }
    }
  }
  SLIC_WARNING_IF(result & EXT_ARRAY_ERROR, "External_array was not correctly loaded");

  result = SPIO_TEST_SUCCESS;

  /*
   * Load view within group_b
   */
  reader.loadExternalData(root2, group_b, file_name + ROOT_EXT);

  /*
   * Verify that the contents of view2 are restored.
   */
  EXPECT_EQ(view2->getNumElements(), nvals);
  if(view2->getNumElements() != nvals)
  {
    result |= EXT_ARRAY_ERROR;
  }
  else
  {
    for(int i = 0; i < nvals; ++i)
    {
      EXPECT_EQ(orig_vals2[i], restored_vals2[i]);
      if(orig_vals2[i] != restored_vals2[i])
      {
        result |= EXT_ARRAY_ERROR;
        break;
      }
    }
  }
  SLIC_WARNING_IF(result & EXT_ARRAY_ERROR, "External_array was not correctly loaded");

  /*
   * Create another DataStore to test another way of reading input
   */
  DataStore* ds3 = new DataStore();
  Group* root3 = ds3->getRoot();

  reader.read(root3, file_name + ROOT_EXT);

  /*
   * Reset the restored_vals arrays to incorrect values.
   */
  for(int i = 0; i < nvals; ++i)
  {
    restored_vals1[i] = -1;
    restored_vals2[i] = -1;
  }

  view1 = root3->getView("testdata/fields/a/vals1");
  view1->setExternalDataPtr(restored_vals1);

  view2 = root3->getView("testdata/fields/b/vals2");
  view2->setExternalDataPtr(restored_vals2);

  /*
   * This uses the API for piecewise loading of external data but calls it
   * to load into the root group.  When loading into the root group, it is
   * unnecessary to use this API, but here we test that it works nonetheless.
   */
  reader.loadExternalData(root3, root3, file_name + ROOT_EXT);

  result = SPIO_TEST_SUCCESS;

  /*
   * Verify that the contents of view1 are restored.
   */
  EXPECT_EQ(view1->getNumElements(), nvals);
  if(view1->getNumElements() != nvals)
  {
    result |= EXT_ARRAY_ERROR;
  }
  else
  {
    for(int i = 0; i < nvals; ++i)
    {
      EXPECT_EQ(orig_vals1[i], restored_vals1[i]);
      if(orig_vals1[i] != restored_vals1[i])
      {
        result |= EXT_ARRAY_ERROR;
        break;
      }
    }
  }
  SLIC_WARNING_IF(result & EXT_ARRAY_ERROR, "External_array was not correctly loaded");

  result = SPIO_TEST_SUCCESS;

  /*
   * Verify that the contents of view2 are restored.
   */
  EXPECT_EQ(view2->getNumElements(), nvals);
  if(view2->getNumElements() != nvals)
  {
    result |= EXT_ARRAY_ERROR;
  }
  else
  {
    for(int i = 0; i < nvals; ++i)
    {
      EXPECT_EQ(orig_vals2[i], restored_vals2[i]);
      if(orig_vals2[i] != restored_vals2[i])
      {
        result |= EXT_ARRAY_ERROR;
        break;
      }
    }
  }
  SLIC_WARNING_IF(result & EXT_ARRAY_ERROR, "External_array was not correctly loaded");

  /*
   * Create another DataStore to test another way of reading input
   */
  DataStore* ds4 = new DataStore();

  /*
   * Write the ds1 data to a different file name, starting at a subgroup,
   * not at the root.
   */
  file_name = "out_spio_external_partial_subgroup";
  writer.write(root1->getGroup("testdata"), num_files, file_name, PROTOCOL);

  Group* root4 = ds4->getRoot();
  Group* testdata = root4->createGroup("testdata");
  reader.read(testdata, file_name + ROOT_EXT);

  /*
   * Reset the restored_vals arrays to incorrect values.
   */
  for(int i = 0; i < nvals; ++i)
  {
    restored_vals1[i] = -1;
    restored_vals2[i] = -1;
  }

  group_a = root4->getGroup("testdata/fields/a");
  group_b = root4->getGroup("testdata/fields/b");

  view1 = group_a->getView("vals1");
  view1->setExternalDataPtr(restored_vals1);

  view2 = group_b->getView("vals2");
  view2->setExternalDataPtr(restored_vals2);

  /*
   * Do a piecewise loadExternalData for group_a.  The first argument must
   * match the Group passed in the 'reader.read()' call.
   */
  reader.loadExternalData(testdata, group_a, file_name + ROOT_EXT);

  result = SPIO_TEST_SUCCESS;

  /*
   * Verify that the contents of view1 are restored.
   */
  EXPECT_EQ(view1->getNumElements(), nvals);
  if(view1->getNumElements() != nvals)
  {
    result |= EXT_ARRAY_ERROR;
  }
  else
  {
    for(int i = 0; i < nvals; ++i)
    {
      EXPECT_EQ(orig_vals1[i], restored_vals1[i]);
      if(orig_vals1[i] != restored_vals1[i])
      {
        result |= EXT_ARRAY_ERROR;
        break;
      }
    }
  }
  SLIC_WARNING_IF(result & EXT_ARRAY_ERROR, "External_array was not correctly loaded");

  result = SPIO_TEST_SUCCESS;

  /*
   * Do a piecewise loadExternalData for group_b.  The first argument must
   * match the Group passed in the 'reader.read()' call.
   */
  reader.loadExternalData(testdata, group_b, file_name + ROOT_EXT);

  /*
   * Verify that the contents of view2 are restored.
   */

  EXPECT_EQ(view2->getNumElements(), nvals);
  if(view2->getNumElements() != nvals)
  {
    result |= EXT_ARRAY_ERROR;
  }
  else
  {
    for(int i = 0; i < nvals; ++i)
    {
      EXPECT_EQ(orig_vals2[i], restored_vals2[i]);
      if(orig_vals2[i] != restored_vals2[i])
      {
        result |= EXT_ARRAY_ERROR;
        break;
      }
    }
  }
  SLIC_WARNING_IF(result & EXT_ARRAY_ERROR, "External_array was not correctly loaded");

  delete ds1;
  delete ds2;
  delete ds3;
  delete ds4;
}

//----------------------------------------------------------------------
TEST(spio_parallel, irregular_writeread)
{
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  const int num_output = numOutputFiles(num_ranks);

  /*
   * Create a DataStore and give it a small hierarchy of groups and views.
   *
   * The views are filled with repeatable nonsense data that will vary based
   * on rank.
   */
  DataStore* ds1 = new DataStore();

  Group* root1 = ds1->getRoot();

  int num_fields = my_rank + 2;

  for(int f = 0; f < num_fields; ++f)
  {
    std::ostringstream ostream;
    ostream << "fields" << f;
    Group* flds = root1->createGroup(ostream.str());

    int num_subgroups = ((f + my_rank) % 3) + 1;
    for(int g = 0; g < num_subgroups; ++g)
    {
      std::ostringstream gstream;
      gstream << "subgroup" << g;
      Group* sg = flds->createGroup(gstream.str());

      std::ostringstream vstream;
      vstream << "view" << g;
      if(g % 2)
      {
        sg->createView(vstream.str())->allocate(DataType::c_int(10 + my_rank));
        int* vals = sg->getView(vstream.str())->getData();

        for(int i = 0; i < 10 + my_rank; i++)
        {
          vals[i] = (i + 10) * (404 - my_rank - i - g - f);
        }
      }
      else
      {
        sg->createViewScalar<int>(vstream.str(), 101 * my_rank * (f + g + 1));
      }
    }
  }

  /*
   * Contents of the DataStore written to files with IOManager.
   */
  int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD);

  const std::string file_name = "out_spio_irregular_write_read";
  writer.write(root1, num_files, file_name, PROTOCOL);

  /*
   * Create another DataStore that holds nothing but the root group.
   */
  DataStore* ds2 = new DataStore();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), file_name + ROOT_EXT);

  /*
   * Verify that the contents of ds2 match those written from ds.
   */
  EXPECT_TRUE(ds2->getRoot()->isEquivalentTo(root1));

  for(int f = 0; f < num_fields; ++f)
  {
    std::ostringstream ostream;
    ostream << "fields" << f;
    Group* flds1 = ds1->getRoot()->getGroup(ostream.str());
    Group* flds2 = ds2->getRoot()->getGroup(ostream.str());

    int num_subgroups = ((f + my_rank) % 3) + 1;
    for(int g = 0; g < num_subgroups; ++g)
    {
      std::ostringstream gstream;
      gstream << "subgroup" << g;
      Group* sg1 = flds1->getGroup(gstream.str());
      Group* sg2 = flds2->getGroup(gstream.str());

      std::ostringstream vstream;
      vstream << "view" << g;
      if(g % 2)
      {
        View* view_orig = sg1->getView(vstream.str());
        View* view_restored = sg2->getView(vstream.str());

        int num_elems = view_orig->getNumElements();
        EXPECT_EQ(view_restored->getNumElements(), num_elems);
        int* vals_orig = view_orig->getData();
        int* vals_restored = view_restored->getData();

        for(int i = 0; i < num_elems; ++i)
        {
          EXPECT_EQ(vals_orig[i], vals_restored[i]);
        }
      }
      else
      {
        int testvalue1 = sg1->getView(vstream.str())->getData();
        int testvalue2 = sg2->getView(vstream.str())->getData();

        EXPECT_EQ(testvalue1, testvalue2);
      }
    }
  }

  delete ds1;
  delete ds2;
}

//------------------------------------------------------------------------------
TEST(spio_parallel, preserve_writeread)
{
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  const int num_output = numOutputFiles(num_ranks);

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

  /*
   * Contents of the DataStore written to files with IOManager.
   */
  int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD);

  const std::string file_name = "out_spio_preserve_write_read";
  writer.write(root, num_files, file_name, PROTOCOL);

  std::string root_name = file_name + ROOT_EXT;

  /*
   * Extra stuff to exercise preserve_contents option
   */
  MPI_Barrier(MPI_COMM_WORLD);
  DataStore* dsextra = new DataStore();
  Group* extra = dsextra->getRoot()->createGroup("extra");
  extra->createViewScalar<double>("dval", 1.1);
  Group* child = extra->createGroup("child");
  child->createViewScalar<int>("ival", 7);
  child->createViewString("word0", "hello");
  child->createViewString("word1", "world");

  const std::string extra_file_name = "out_spio_extra";
  writer.write(extra, num_files, extra_file_name, PROTOCOL);
  std::string extra_root = extra_file_name + ROOT_EXT;

  MPI_Barrier(MPI_COMM_WORLD);

  /*
   * Create another DataStore that holds nothing but the root group.
   */
  DataStore* ds2 = new DataStore();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), root_name);

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

  /*
   * Read in extra file while preserving contents
   */
  Group* extra_fields = ds2->getRoot()->getGroup("fields");
  reader.read(extra_fields, extra_root, true);

  /*
   * Test one of the pre-existing Views to show it's unchanged.
   */
  int testvalue_extra = extra_fields->getGroup("a")->getView("i0")->getData();
  EXPECT_EQ(testvalue_extra, testvalue2);

  /*
   * Test the data from the extra file.
   */
  double dval = extra_fields->getView("dval")->getData();

  EXPECT_TRUE(dval > 1.0000009 && dval < 1.1000001);

  int ival = extra_fields->getView("child/ival")->getData();
  EXPECT_EQ(ival, 7);

  EXPECT_EQ(std::string(extra_fields->getView("child/word0")->getString()), "hello");

  EXPECT_EQ(std::string(extra_fields->getView("child/word1")->getString()), "world");

  delete ds;
  delete ds2;
  delete dsextra;
}

TEST(spio_parallel, parallel_increase_procs)
{
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  EXPECT_TRUE(num_ranks > my_rank);

#ifdef AXOM_USE_HDF5
  /*
   * This tests the ability to read data when the run that is doing the
   * reading has more processors than the run that created the files
   * being read.  In this test we split the world communicator into
   * smaller communicators to simulate a run that dumps data from
   * only ranks 0 and 1 (only rank 0 if running on <= 2 processors).
   * Then the read occurs with an IOManager constructed with
   * MPI_COMM_WORLD.
   *
   * This functionality only works with the sidre_hdf5 protocol, so this
   * section of code is inside the #ifdef guards.
   */

  int top_output_rank = 1;
  if(num_ranks <= 2)
  {
    top_output_rank = 0;
  }

  // Split the communicator so that ranks up to and including
  // top_output_rank have their own communicator for the output step.
  MPI_Comm split_comm;

  if(my_rank <= top_output_rank)
  {
    MPI_Comm_split(MPI_COMM_WORLD, 0, my_rank, &split_comm);
  }
  else
  {
    MPI_Comm_split(MPI_COMM_WORLD, my_rank, 0, &split_comm);
  }

  const std::string file_name = "out_spio_parallel_increase_procs";

  DataStore* ds = new DataStore();
  if(my_rank <= top_output_rank)
  {
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

    int num_files = 1;
    axom::sidre::IOManager writer(split_comm);

    writer.write(root, num_files, file_name, PROTOCOL);
  }

  /*
   * The reading section of this test will execute a read on all
   * ranks of MPI_COMM_WORLD, even though the write was only on one or two
   * ranks.  The read will load data on the ranks that wrote data and add
   * nothing on higher ranks.
   */
  DataStore* ds2 = new DataStore();

  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), file_name + ROOT_EXT);

  /*
   * Verify that the contents of ds2 on rank 0 match those written from ds.
   */

  if(my_rank <= top_output_rank)
  {
    EXPECT_TRUE(my_rank != 0 || ds2->getRoot()->isEquivalentTo(ds->getRoot()));

    int testvalue = ds->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();
    int testvalue2 = ds2->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();

    EXPECT_EQ(testvalue, testvalue2);

    View* view_i1_orig = ds->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");
    View* view_i1_restored = ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");

    int num_elems = view_i1_orig->getNumElements();
    EXPECT_EQ(view_i1_restored->getNumElements(), num_elems);
    if(view_i1_restored->getNumElements() == num_elems)
    {
      int* i1_orig = view_i1_orig->getData();
      int* i1_restored = view_i1_restored->getData();

      for(int i = 0; i < num_elems; ++i)
      {
        EXPECT_EQ(i1_orig[i], i1_restored[i]);
      }
    }
  }
  else
  {
    EXPECT_FALSE(ds2->getRoot()->hasGroup("fields"));
    EXPECT_FALSE(ds2->getRoot()->hasGroup("fields2"));
    EXPECT_EQ(ds2->getRoot()->getNumGroups(), 0);
    EXPECT_EQ(ds2->getRoot()->getNumViews(), 0);
  }

  delete ds2;
  delete ds;

  MPI_Comm_free(&split_comm);

#endif
}

TEST(spio_parallel, parallel_decrease_procs)
{
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  EXPECT_TRUE(num_ranks > my_rank);

#ifdef AXOM_USE_HDF5
  /*
   * This tests the ability to read data when the run that is doing the
   * reading has fewer processors than the run that created the files
   * being read.  In this test we dump data using all ranks on the world
   * communicator, and then to read the data we split the world communicator
   * into smaller communicators in order to read data using only ranks 0
   * and 1 (only rank 0 if running on <= 2 processors).
   *
   * This functionality only works with the sidre_hdf5 protocol, so this
   * section of code is inside the #ifdef guards.
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

  int num_files = num_ranks;
  axom::sidre::IOManager writer(MPI_COMM_WORLD);

  const std::string file_name = "out_spio_parallel_decrease_procs";

  writer.write(root, num_files, file_name, PROTOCOL);

  int top_input_rank = 1;
  if(num_ranks <= 2)
  {
    top_input_rank = 0;
  }

  // Split the communicator so that ranks up to and including
  // top_input_rank have their own communicator for the input step.
  MPI_Comm split_comm;

  if(my_rank <= top_input_rank)
  {
    MPI_Comm_split(MPI_COMM_WORLD, 0, my_rank, &split_comm);
  }
  else
  {
    MPI_Comm_split(MPI_COMM_WORLD, my_rank, 0, &split_comm);
  }

  /*
   * The reading section of this test will execute a read on ranks 0 and 1,
   * or only 0 if running on <= 2 processors.
   */
  DataStore* ds2 = new DataStore();

  if(my_rank <= top_input_rank)
  {
    IOManager reader(split_comm);

    reader.read(ds2->getRoot(), file_name + ROOT_EXT);

    Group* ds2_root = ds2->getRoot();

    int num_output_ranks = 1;
    if(num_ranks > 1)
    {
      EXPECT_TRUE(ds2_root->hasView("reduced_input_ranks"));
    }

    if(ds2->getRoot()->hasView("reduced_input_ranks"))
    {
      num_output_ranks = ds2_root->getView("reduced_input_ranks")->getData();
    }

    for(int output_rank = my_rank; output_rank < num_output_ranks;
        output_rank += (top_input_rank + 1))
    {
      /*
       * Verify that the contents of ds2 on rank 0 match those written from ds.
       */
      std::string output_name = axom::fmt::sprintf("rank_%07d/sidre_input", output_rank);

      int testvalue = 101 * output_rank;
      int testvalue2 =
        ds2_root->getGroup(output_name)->getGroup("fields")->getGroup("a")->getView("i0")->getData();

      EXPECT_EQ(testvalue, testvalue2);

      View* view_i1_orig = ds->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");
      View* view_i1_restored =
        ds2_root->getGroup(output_name)->getGroup("fields2")->getGroup("b")->getView("i1");

      int num_elems = view_i1_orig->getNumElements();
      EXPECT_EQ(view_i1_restored->getNumElements(), num_elems);
      if(view_i1_restored->getNumElements() == num_elems)
      {
        int* i1_restored = view_i1_restored->getData();

        for(int i = 0; i < num_elems; ++i)
        {
          EXPECT_EQ((i + 10) * (404 - output_rank - i), i1_restored[i]);
        }
      }
    }
  }
  else
  {
    EXPECT_EQ(ds2->getRoot()->getNumGroups(), 0);
    EXPECT_EQ(ds2->getRoot()->getNumViews(), 0);
  }

  delete ds2;
  delete ds;

  MPI_Comm_free(&split_comm);

#endif  // AXOM_USE_HDF5
}

TEST(spio_parallel, sidre_simple_blueprint_example)
{
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  EXPECT_TRUE(num_ranks > my_rank);

#ifdef AXOM_USE_HDF5
  /*
   This tests creates a simple sidre + mesh blueprint example.

   It is run with 4 MPI tasks, each task creates a domain
   with a very small uniform grid. 

  */

  // create a mesh using conduit's examples
  axom::sidre::Node n_mesh;
  // create a 2x2 element uniform mesh for each domain
  conduit::blueprint::mesh::examples::basic("uniform", 3, 3, 1, n_mesh);

  // shift example mesh in x using rank to create non-overlapping domains
  n_mesh["coordsets/coords/origin/x"] = my_rank * 20.0;

  // add an element-assoced field that contains the rank
  n_mesh["fields/rank/association"] = "element";
  n_mesh["fields/rank/topology"] = "mesh";
  n_mesh["fields/rank/values"].set(conduit::DataType::int64(4));

  // fill rank field values
  std::int64_t* rank_vals_ptr = n_mesh["fields/rank/values"].value();

  for(int i = 0; i < 4; i++)
  {
    rank_vals_ptr[i] = my_rank;
  }

  // setup the data store
  DataStore ds;
  Group* root = ds.getRoot();
  Group* domain_root = root->createGroup("mesh");
  // setup sidre group using the conduit tree
  domain_root->importConduitTree(n_mesh);

  // use SPIO to write to a set of hdf5 files
  // NOTE:
  // I (cyrush) tried using json protocols but hit issues in VisIt.
  // VisIt is tested against general json bp files, but not sidre
  // flavored json, something we should look into.

  axom::sidre::IOManager writer(MPI_COMM_WORLD);
  const std::string file_name = "out_spio_blueprint_example";
  writer.write(root, 4, file_name, "sidre_hdf5");

  // make sure all tasks are done writing before we try
  // to add the blueprint index to the root file
  MPI_Barrier(MPI_COMM_WORLD);

  // Add the bp index to the root file
  writer.writeBlueprintIndexToRootFile(&ds, "mesh", file_name + ROOT_EXT, "mesh");

#endif  // AXOM_USE_HDF5
}
