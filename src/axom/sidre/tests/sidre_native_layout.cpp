// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "conduit_blueprint.hpp"

#include <vector>
#include <string>
#include <iostream>

using axom::sidre::Buffer;
using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::indexIsValid;
using axom::sidre::IndexType;
using axom::sidre::InvalidIndex;
using axom::sidre::InvalidName;
using axom::sidre::nameIsValid;
using axom::sidre::View;

//------------------------------------------------------------------------------

namespace
{
const axom::sidre::TypeID DOUBLE_ID = axom::sidre::DOUBLE_ID;
const axom::sidre::TypeID INT32_ID = axom::sidre::INT32_ID;

/**
 * \brief Simple utility function to initialize an array
 * \param data A pointer to the beginning of the array
 * \param size The number of elements in the array
 * \param initVal The starting value of the element
 * \param intDiv An integer divisor for the data's index
 * \param scaleFac A scaling factor for the data's index
 *
 * The value of element idx will be initVal + ( idx / intDiv) * scaleFac
 */
template <typename T>
void setData(T* data, int size, T initVal = T(0), int intDiv = 1, T scaleFac = T(1))
{
  for(int i = 0; i < size; ++i)
  {
    data[i] = initVal + ((i + 1) / intDiv) * scaleFac;
  }
}

/**
 * \brief Utility function to compare the corresponding node and view
 * \param path The path from the datastore root to the view
 * \param rootNode The root node of the layout
 * \param rootGroup The root group of the datastore
 */
template <typename T>
void checkPointersAndData(const std::string& path,
                          axom::sidre::Node& rootNode,
                          axom::sidre::Group* rootGroup)
{
  axom::sidre::Node& node = rootNode[path];
  T* nD = static_cast<T*>(node.element_ptr(0));

  axom::sidre::View* view = rootGroup->getView(path);
  EXPECT_TRUE(nullptr != view);

  T* vD = view->getData<T*>();
  EXPECT_TRUE(nullptr != vD);

  EXPECT_EQ(nD, vD) << "Error when comparing pointer address between "
                    << "conduit native layout and Sidre View "
                    << "for path " << path << ".\n\t Conduit address: " << nD
                    << "\n\t Datastore view address: " << vD;

  EXPECT_EQ(*nD, *vD) << "Error when comparing values between "
                      << "conduit native layout and Sidre View "
                      << "for path " << path << ".\n\t Conduit value: " << *nD
                      << "\n\t View value: " << *vD;
}

/**
 * \brief Template specialization of checkPointersAndData for scalar views
 *        containing strings
 */
template <>
void checkPointersAndData<std::string>(const std::string& path,
                                       axom::sidre::Node& rootNode,
                                       axom::sidre::Group* rootGroup)
{
  axom::sidre::Node& node = rootNode[path];
  std::string nD = node.as_string();

  axom::sidre::View* view = rootGroup->getView(path);
  EXPECT_TRUE(nullptr != view);
  EXPECT_TRUE(view->isString());

  std::string vD(view->getString());

  EXPECT_EQ(nD, vD) << "Error when comparing strings between "
                    << "conduit native layout and Sidre View "
                    << "for path " << path << ".\n\t Conduit string: " << nD
                    << "\n\t Datastore view string: " << vD;
}
}  // namespace

TEST(sidre_native_layout, empty_layout)
{
  DataStore* ds = new DataStore();

  SLIC_INFO("Testing sidre_layout function on empty datastore");

  SLIC_INFO("****** Visual diagnostics ******");
  SLIC_INFO("Printing datastore");
  ds->print();

  SLIC_INFO("Printing datastore native layout:");
  axom::sidre::Node node;
  ds->getRoot()->createNativeLayout(node);
  node.to_json_stream(std::cout);

  SLIC_INFO("****** done ******");

  if(ds != nullptr)
  {
    delete ds;
    ds = nullptr;
  }
}

TEST(sidre_native_layout, generate_native_layout)
{
  /// Allocate and initialize two external arrays
  // The REAL_BUF has 10 copies of each integer
  const int EXT_REAL_BUF_SIZE = 20;
  double* extRealPtr = new double[EXT_REAL_BUF_SIZE];
  setData(extRealPtr, EXT_REAL_BUF_SIZE, 0., 1, 1 / 10.);

  // The INT_BUF is a decreasing sequence
  const int EXT_INT_BUF_SIZE = 5;
  int* extIntPtr = new int[EXT_INT_BUF_SIZE];
  setData(extIntPtr, EXT_INT_BUF_SIZE, EXT_INT_BUF_SIZE, 1, -1);

  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Setup a buffer that will have two views
  const int REAL_BUF_SIZE = 100;
  Buffer* realBuf = ds->createBuffer(DOUBLE_ID, REAL_BUF_SIZE)->allocate();
  setData<double>(realBuf->getData(), REAL_BUF_SIZE, 0., 10, 1.);

  // create the views using the path syntax
  root->createView("Ga/Va")->attachBuffer(realBuf)->apply(50, 0);
  root->createView("Ga/Vb")->attachBuffer(realBuf)->apply(20, 50);

  root->createView("Ga/Gc/Ve", DOUBLE_ID, 10)->allocate();
  root->createView("Ga/Gc/Vf")->setExternalDataPtr(DOUBLE_ID, EXT_REAL_BUF_SIZE, extRealPtr);

  root->createView("Gb/Vc", INT32_ID, 10)->allocate();
  root->createView("Gb/Vd")->setExternalDataPtr(INT32_ID, EXT_INT_BUF_SIZE, extIntPtr);

  // Set some data in the attached buffers so we can distinguish in our tests
  setData<double>(root->getView("Ga/Gc/Ve")->getData(), 10, 0., 1, 1.1);
  setData<int>(root->getView("Gb/Vc")->getData(), 10, 0, 1, 2);

  SLIC_INFO("****** Visual diagnostics ******");
  SLIC_INFO("Printing datastore");
  ds->print();

  SLIC_INFO("Printing datastore native layout:");
  axom::sidre::Node node;
  ds->getRoot()->createNativeLayout(node);
  node.to_json_stream(std::cout);
  //node.save("nativeLayoutTest.conduit");

  SLIC_INFO("done.");

  SLIC_INFO("*** Checking that addresses and first elements match.");
  checkPointersAndData<double>("Ga/Va", node, root);
  checkPointersAndData<double>("Ga/Va", node, root);
  checkPointersAndData<double>("Ga/Gc/Ve", node, root);
  checkPointersAndData<double>("Ga/Gc/Vf", node, root);
  checkPointersAndData<int>("Gb/Vc", node, root);
  checkPointersAndData<int>("Gb/Vd", node, root);
  SLIC_INFO("done.");

  /// Clean up memory
  if(extRealPtr != nullptr)
  {
    delete[] extRealPtr;
    extRealPtr = nullptr;
  }
  if(extIntPtr != nullptr)
  {
    delete[] extIntPtr;
    extIntPtr = nullptr;
  }

  if(ds != nullptr)
  {
    delete ds;
    ds = nullptr;
  }
}

TEST(sidre_native_layout, native_layout_with_scalars)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  root->createView("Garray/Vdbl20", DOUBLE_ID, 20)->allocate();
  root->createView("Garray/Vint10", INT32_ID, 10)->allocate();

  // Set some data in the attached buffers so we can distinguish in our tests
  setData<double>(root->getView("Garray/Vdbl20")->getData(), 20, 0., 1, 1.1);
  setData<int>(root->getView("Garray/Vint10")->getData(), 10, 0, 1, 2);

  root->createView("Gscalar/Vint")->setScalar(1);
  root->createView("Gscalar/Vdbl")->setScalar(2.);
  root->createView("Gscalar/Vstr")->setString("Hello");

  SLIC_INFO("****** Visual diagnostics ******");
  SLIC_INFO("Printing datastore::info():");
  ds->print();

  SLIC_INFO("Printing datastore::layout():");
  axom::sidre::Node node;
  ds->getRoot()->createNativeLayout(node);
  node.to_json_stream(std::cout);
  //node.save("nativeLayoutTest.conduit");

  SLIC_INFO("done.");
  SLIC_CHECK_MSG(false, "Some message");

  SLIC_INFO("*** Checking that addresses and first elements match.");
  checkPointersAndData<double>("Garray/Vdbl20", node, root);
  checkPointersAndData<int>("Garray/Vint10", node, root);
  checkPointersAndData<int>("Gscalar/Vint", node, root);
  checkPointersAndData<double>("Gscalar/Vdbl", node, root);
  checkPointersAndData<std::string>("Gscalar/Vstr", node, root);

  if(ds != nullptr)
  {
    delete ds;
    ds = nullptr;
  }
}

//------------------------------------------------------------------------------
TEST(sidre_native_layout, export_import_conduit)
{
  // Checks that we can import an exported conduit tree into a separate
  // datastore using Group::createNativeLayout() and Group::importConduitTree()

  // Tests with non-zero- and with zero-sized arrays
  for(int SZ : {10, 0})
  {
    SLIC_INFO("Checking export/import from conduit with array of size " << SZ);

    axom::sidre::DataStore ds1, ds2;
    conduit::Node node1, node2;

    int* conduit_data = nullptr;

    std::string grpName = "some_group";
    std::string viewName = "array_view";

    // Simple lambda to return the expected value at an index
    auto expValue = AXOM_HOST_LAMBDA(int i) { return SZ - i; };

    // Create datastore with an array view and export to conduit node
    {
      // Create views and groups
      auto* root = ds1.getRoot();
      auto* grp = root->createGroup(grpName);
      auto* view = grp->createViewAndAllocate(viewName, INT32_ID, SZ);

      // initialize the data
      int* sidre_data = view->getData();
      EXPECT_NE(nullptr, sidre_data);
      for(int i = 0; i < SZ; ++i)
      {
        sidre_data[i] = expValue(i);
      }

      // Dump to console
      std::cout << "Datastore1 after initialization: " << std::endl;
      root->print();
      std::cout << std::endl;

      // Convert datastore to conduit layout and compare data
      bool success = root->createNativeLayout(node1);
      EXPECT_TRUE(success);

      conduit_data = node1[grpName][viewName].as_int_ptr();
      EXPECT_NE(nullptr, conduit_data);
      for(int i = 0; i < SZ; ++i)
      {
        EXPECT_EQ(conduit_data[i], sidre_data[i]);
      }

      // Dump to console
      std::cout << "Node1 after createNativeLayout: " << std::endl;
      node1.print();
      std::cout << std::endl;
    }

    // import conduit node into new datastore and check validity
    {
      auto* root = ds2.getRoot();
      bool success = root->importConduitTree(node1);
      EXPECT_TRUE(success);

      // Dump to console
      std::cout << "Datastore2 after importing conduit tree: " << std::endl;
      root->print();
      std::cout << std::endl;

      // Check for expected views and groups
      EXPECT_TRUE(root->hasGroup(grpName));
      auto* grp = root->getGroup(grpName);
      EXPECT_NE(nullptr, grp);

      EXPECT_TRUE(grp->hasView(viewName));
      auto* view = grp->getView(viewName);
      EXPECT_NE(nullptr, view);

      EXPECT_TRUE(view->isAllocated());
      EXPECT_EQ(SZ, view->getNumElements());

      // compare data
      int* sidre_data = view->getData();
      EXPECT_NE(nullptr, sidre_data);
      for(int i = 0; i < SZ; ++i)
      {
        EXPECT_EQ(conduit_data[i], sidre_data[i]);
      }

      // confirm that this is a distinct copy of the first data store's data
      auto* otherView = ds1.getRoot()->getGroup(grpName)->getView(viewName);
      EXPECT_NE(sidre_data, otherView->getData<int*>());

      // Convert to conduit layout and compare data
      success = root->createNativeLayout(node2);
      EXPECT_TRUE(success);

      conduit_data = node2[grpName][viewName].as_int_ptr();
      EXPECT_NE(nullptr, conduit_data);
      for(int i = 0; i < SZ; ++i)
      {
        EXPECT_EQ(conduit_data[i], sidre_data[i]);
      }

      // Dump to console
      std::cout << "Node2 after createNativeLayout: " << std::endl;
      node2.print();
      std::cout << std::endl;
    }

    // compare the two conduit nodes
    // first, check compatiblity
    bool compat = node1.compatible(node2);
    EXPECT_TRUE(compat);

    // Finally, check the actual array data
    {
      int* arr1 = node1[grpName][viewName].as_int_ptr();
      EXPECT_NE(nullptr, arr1);
      int* arr2 = node2[grpName][viewName].as_int_ptr();
      EXPECT_NE(nullptr, arr2);

      for(int i = 0; i < SZ; ++i)
      {
        EXPECT_EQ(arr1[i], arr2[i]);
        EXPECT_EQ(expValue(i), arr2[i]);
      }
    }
  }
}

TEST(sidre_native_layout, import_conduit_and_verify_protocol)
{
  // This tests is similar to the 'export_import_conduit' test above,
  // but uses conduit's blueprint_verify on the result
  // which is expected to be a blueprint multicomponent array (mcarray)

  std::string grpName = "multicomponent_array";
  const int NCOMP = 5;  // number of components in the array

  // Tests with non-zero- and with zero-sized arrays
  for(int SZ : {10, 0})
  {
    SLIC_INFO("Checking export of blueprint protocol for 'mcarray' with "
              << NCOMP << " components, each of size " << SZ);

    axom::sidre::DataStore ds1, ds2;
    conduit::Node node1, node2;

    // Simple lambda to get/check expected value at index (i,j)
    auto expValue = AXOM_HOST_LAMBDA(int i, int j) { return j * SZ + (SZ - i); };
    // Simple lambda to get the name of the i^th component group
    auto viewName = AXOM_HOST_LAMBDA(int j)
    {
      return std::string(1, 'a' + j);  // i.e. 'a', 'b', 'c', ...
    };

    // Initialize datastore, export to conduit and verify 'mcarray
    {
      auto* root = ds1.getRoot();
      auto* grp = root->createGroup(grpName);
      for(int j = 0; j < NCOMP; ++j)
      {
        auto* view = grp->createViewAndAllocate(viewName(j), INT32_ID, SZ);

        int* sidre_data = view->getData();
        EXPECT_NE(nullptr, sidre_data);
        for(int i = 0; i < SZ; ++i)
        {
          sidre_data[i] = expValue(i, j);
        }
      }

      // Convert to conduit
      bool success = root->createNativeLayout(node1);
      EXPECT_TRUE(success);

      std::cout << "Conduit native layout of datastore after exporting: \n";
      node1.print();
      std::cout << std::endl;

      // Call blueprint::verify to check that it complies with 'mcarray' schema
      conduit::Node info;
      bool verified = conduit::blueprint::verify("mcarray", node1[grpName], info);
      EXPECT_TRUE(verified);

      std::cout << "Blueprint verify info for 'mcarray' protocol on node 1: \n";
      info.print();
      std::cout << std::endl;
    }

    // Import into sidre, export again to conduit, then verify blueprint
    {
      // Import into sidre
      auto* root = ds1.getRoot();
      bool success = root->importConduitTree(node1);
      EXPECT_TRUE(success);

      // export back to conduit
      success = root->createNativeLayout(node2);
      EXPECT_TRUE(success);

      std::cout << "Conduit layout after export-import-export from Sidre: \n";
      node2.print();
      std::cout << std::endl;

      // verify the mcarray protocol
      conduit::Node info;
      bool verified = conduit::blueprint::verify("mcarray", node2[grpName], info);
      EXPECT_TRUE(verified);

      std::cout << "Blueprint verify info for 'mcarray' protocol on node2: \n";
      info.print();
      std::cout << std::endl;
    }

    // Check the results
    {
      // Verify the two blueprints are compatible
      bool compat = node1.compatible(node2);
      EXPECT_TRUE(compat);

      // Check the data in the final conduit node
      conduit::Node& mcarr = node2[grpName];
      for(int j = 0; j < NCOMP; ++j)
      {
        std::string name = viewName(j);

        EXPECT_TRUE(mcarr.has_child(name));
        EXPECT_EQ(SZ, mcarr[name].dtype().number_of_elements());

        int* arr = mcarr[name].as_int_ptr();
        EXPECT_NE(nullptr, arr);

        for(int i = 0; i < SZ; ++i)
        {
          EXPECT_EQ(expValue(i, j), arr[i]);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

TEST(sidre_native_layout, basic_demo_compare)
{
  //
  // create demo a datastore that includes:
  //  groups
  //  views:
  //    scalar, string, array, external array, array of sliced buffer

  DataStore ds;
  Group* root = ds.getRoot();

  // Note:
  // we use int64 and float64 b/c those types persist even with
  // json or yaml output

  std::int64_t sidre_vals_1[5] = {0, 1, 2, 3, 4};
  axom::float64 sidre_vals_2[6] = {
    1.0,
    2.0,
    1.0,
    2.0,
    1.0,
    2.0,
  };

  Group* group1 = root->createGroup("my_scalars");
  // two scalars
  group1->createViewScalar<conduit::int64>("i64", 1);
  group1->createViewScalar<conduit::float64>("f64", 10.0);

  Group* group2 = root->createGroup("my_strings");
  // two strings
  group2->createViewString("s0", "s0 string");
  group2->createViewString("s1", "s1 string");

  // one basic array
  Group* group3 = root->createGroup("my_arrays");
  View* a5_view = group3->createViewAndAllocate("a5_i64", axom::sidre::DataType::int64(5));

  std::int64_t* a5_view_ptr = a5_view->getData();

  for(int i = 0; i < 5; i++)
  {
    a5_view_ptr[i] = sidre_vals_1[i];
  }

  // one array of size 0
  group3->createViewAndAllocate("a0_i64", axom::sidre::DataType::int64(0));

  // one external array
  group3->createView("a5_i64_ext", conduit::DataType::int64(5))->setExternalDataPtr(sidre_vals_1);

  // change val to make sure this is reflected as external
  sidre_vals_1[4] = -5;

  // create a buffer that we will access via non std offset + stride combo
  Buffer* buff = ds.createBuffer(axom::sidre::FLOAT64_ID, 6);
  buff->allocate();

  conduit::float64* buff_ptr = buff->getData();

  for(int i = 0; i < 6; i++)
  {
    buff_ptr[i] = sidre_vals_2[i];
  }

  View* v0 = group3->createView("b_v0");
  View* v1 = group3->createView("b_v1");
  View* v2 = group3->createView("b_v2");

  // v0 is a view of size zero into this buffer
  v0->attachBuffer(buff);
  v0->apply(conduit::DataType::float64(0));

  // with these settings, bv1 should have 1.0 as all vals
  v1->attachBuffer(buff);
  v1->apply(conduit::DataType::float64(3, 0, 2 * sizeof(conduit::float64)));

  // with these settings, bv2 should have 2.0 as all vals
  v2->attachBuffer(buff);
  v2->apply(conduit::DataType::float64(3, sizeof(conduit::float64), 2 * sizeof(conduit::float64)));
  //
  // save the data store using several protocols
  //
  std::vector<std::string> protocols;

#ifdef AXOM_USE_HDF5
  protocols.push_back("sidre_hdf5");
#endif
  protocols.push_back("sidre_json");
  protocols.push_back("sidre_conduit_json");

  for(size_t i = 0; i < protocols.size(); ++i)
  {
    SLIC_INFO("Testing protocol: " << protocols[i]);
    const std::string file_path = "texample_sidre_basic_ds_demo." + protocols[i];
    // save using current protocol
    ds.getRoot()->save(file_path, protocols[i]);
  }

  axom::sidre::Node n_sidre;
  ds.getRoot()->createNativeLayout(n_sidre);
  SLIC_INFO("Sidre Conduit Native Layout Tree:");
  n_sidre.print();

  //
  // create an equiv conduit tree for testing
  //

  std::int64_t conduit_vals_1[5] = {0, 1, 2, 3, 4};
  axom::float64 conduit_vals_2[6] = {
    1.0,
    2.0,
    1.0,
    2.0,
    1.0,
    2.0,
  };
  // Note: use a std::vector for the empty array as workaround for
  // MSVC error C2466 about allocating arrays of constant size 0
  std::vector<conduit::int64> conduit_vals_0;

  axom::sidre::Node n;
  n["my_scalars/i64"].set_int64(1);
  n["my_scalars/f64"].set_float64(10.0);
  n["my_strings/s0"] = "s0 string";
  n["my_strings/s1"] = "s1 string";
  n["my_arrays/a0_i64"].set(conduit_vals_0.data(), 0);
  n["my_arrays/a5_i64"].set(conduit_vals_1, 5);
  n["my_arrays/a5_i64_ext"].set_external(conduit_vals_1, 5);
  n["my_arrays/b_v0"].set(conduit_vals_2, 0);
  n["my_arrays/b_v1"].set(conduit_vals_2, 3, 0, 2 * sizeof(conduit::float64));
  n["my_arrays/b_v2"].set(conduit_vals_2, 3, sizeof(conduit::float64), 2 * sizeof(conduit::float64));

  // change val to make sure this is reflected as external
  conduit_vals_1[4] = -5;
  SLIC_INFO("Conduit Test Tree:");
  n.print();

  axom::sidre::Node n_info;
  EXPECT_FALSE(n.diff(n_sidre, n_info));
  n_info.print();

  //
  // TODO: When using newer conduit that has relay support for sidre i/o
  //       test round trip with relay here.
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Debug);

  result = RUN_ALL_TESTS();

  return result;
}
