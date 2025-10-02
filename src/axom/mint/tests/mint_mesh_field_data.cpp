// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/slic.hpp"
#include "axom/mint/mesh/FieldData.hpp"         // for mint::FieldData
#include "axom/mint/mesh/FieldAssociation.hpp"  // for FieldAssociation enum

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <algorithm>  // for std::fill()

// namespace aliases
namespace mint = axom::mint;
#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
  #include "axom/mint/deprecated/SidreMCArray.hpp"
namespace sidre = axom::sidre;
#endif

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
template <int ASSOCIATION>
void check_empty_field_data(const mint::FieldData& fd)
{
  EXPECT_TRUE(fd.empty());
  EXPECT_TRUE(fd.getNumFields() == 0);
  EXPECT_TRUE(fd.getAssociation() == ASSOCIATION);
}

//------------------------------------------------------------------------------
void check_resize(mint::FieldData& field_data, int NEW_NUM_TUPLES)
{
  EXPECT_TRUE(field_data.getNumFields() > 1);

  const int numFields = field_data.getNumFields();

  // ensure current number of tuples is smaller the NEW_NUM_TUPLES
  for(int i = 0; i < numFields; ++i)
  {
    mint::Field* f = field_data.getField(i);
    EXPECT_TRUE(f != nullptr);
    EXPECT_TRUE(f->getNumTuples() < NEW_NUM_TUPLES);
  }

  // resize
  field_data.resize(NEW_NUM_TUPLES);

  // check new size
  for(int i = 0; i < numFields; ++i)
  {
    mint::Field* f = field_data.getField(i);
    EXPECT_TRUE(f != nullptr);
    EXPECT_EQ(f->getNumTuples(), NEW_NUM_TUPLES);
  }
}

//------------------------------------------------------------------------------
void check_reserve(mint::FieldData& field_data, int NEW_CAPACITY)
{
  EXPECT_TRUE(field_data.getNumFields() > 1);

  const int numFields = field_data.getNumFields();

  // ensure current number of tuples is smaller the NEW_NUM_TUPLES
  for(int i = 0; i < numFields; ++i)
  {
    mint::Field* f = field_data.getField(i);
    EXPECT_TRUE(f != nullptr);
    EXPECT_TRUE(f->getCapacity() < NEW_CAPACITY);
  }

  // reserve more space
  field_data.reserve(NEW_CAPACITY);

  // check new size
  for(int i = 0; i < numFields; ++i)
  {
    mint::Field* f = field_data.getField(i);
    EXPECT_TRUE(f != nullptr);
    EXPECT_EQ(f->getCapacity(), NEW_CAPACITY);
  }
}

//------------------------------------------------------------------------------
void check_shrink(mint::FieldData& field_data)
{
  EXPECT_TRUE(field_data.getNumFields() > 1);

  const int numFields = field_data.getNumFields();

  // ensure the current capacity is larger than num_tuples
  axom::IndexType num_tuples = field_data.getField(0)->getNumTuples();
  field_data.reserve(num_tuples * 10);
  for(int i = 0; i < numFields; ++i)
  {
    mint::Field* f = field_data.getField(i);
    EXPECT_TRUE(f != nullptr);
    EXPECT_TRUE(f->getNumTuples() < f->getCapacity());
    EXPECT_EQ(f->getCapacity(), num_tuples * 10);
  }

  // shrink the space
  field_data.shrink();

  // check that max capacity is num_tuples
  for(int i = 0; i < numFields; ++i)
  {
    mint::Field* f = field_data.getField(i);
    EXPECT_TRUE(f != nullptr);
    EXPECT_EQ(f->getCapacity(), f->getNumTuples());
  }
}
//------------------------------------------------------------------------------
void check_create_and_access_data(mint::FieldData& field_data,
                                  int NUM_TUPLES,
                                  int NUM_COMPONENTS,
                                  int MAGIC_INT,
                                  double MAGIC_DOUBLE)
{
  // Ensure the data is empty
  EXPECT_TRUE(field_data.empty());

  // create a scalar int field to test with and fill it with MAGIC_INT
  int* f1 = field_data.createField<int>("f1", NUM_TUPLES);
  std::fill(f1, f1 + NUM_TUPLES, MAGIC_INT);

  // create a real-valued vector field, f2, and fill it with MAGIC_DOUBLE
  double* f2 = field_data.createField<double>("f2", NUM_TUPLES, NUM_COMPONENTS);
  std::fill(f2, f2 + (NUM_TUPLES * NUM_COMPONENTS), MAGIC_DOUBLE);

  // check f1 pointer access and parameters
  axom::IndexType N = 0;
  axom::IndexType M = 0;
  int* f1ptr = field_data.getFieldPtr<int>("f1", N, M);
  EXPECT_TRUE(f1ptr != nullptr);
  EXPECT_EQ(N, NUM_TUPLES);
  EXPECT_EQ(M, 1);
  EXPECT_EQ(f1ptr, f1);

  // check f1 field object access by name
  mint::Field* field_by_name = field_data.getField("f1");
  EXPECT_TRUE(field_by_name != nullptr);
  EXPECT_EQ(field_by_name->getName(), "f1");
  EXPECT_EQ(field_by_name->getType(), mint::INT32_FIELD_TYPE);
  EXPECT_EQ(mint::Field::getDataPtr<int>(field_by_name), f1);
  EXPECT_EQ(field_by_name->getNumTuples(), NUM_TUPLES);
  EXPECT_EQ(field_by_name->getNumComponents(), 1);

  // check f1 field object access by index
  mint::Field* field_by_index = field_data.getField(0);
  EXPECT_TRUE(field_by_index != nullptr);
  EXPECT_EQ(field_by_index, field_by_name);
  EXPECT_EQ(field_by_index->getName(), "f1");

  // check f2 pointer access and parameters
  double* f2ptr = field_data.getFieldPtr<double>("f2", N, M);
  EXPECT_TRUE(f2ptr != nullptr);
  EXPECT_EQ(N, NUM_TUPLES);
  EXPECT_EQ(M, NUM_COMPONENTS);
  EXPECT_EQ(f2ptr, f2);

  // check f2 field object access by name
  mint::Field* field2_by_name = field_data.getField("f2");
  EXPECT_TRUE(field2_by_name != nullptr);
  EXPECT_EQ(field2_by_name->getName(), "f2");
  EXPECT_EQ(field2_by_name->getType(), mint::DOUBLE_FIELD_TYPE);
  EXPECT_EQ(mint::Field::getDataPtr<double>(field2_by_name), f2);
  EXPECT_EQ(field2_by_name->getNumTuples(), NUM_TUPLES);
  EXPECT_EQ(field2_by_name->getNumComponents(), NUM_COMPONENTS);

  // check f2 field object access by index
  mint::Field* field2_by_index = field_data.getField(1);
  EXPECT_TRUE(field2_by_index != nullptr);
  EXPECT_EQ(field2_by_index, field2_by_name);
  EXPECT_EQ(field2_by_index->getName(), "f2");

  // check the data
  for(int i = 0; i < NUM_TUPLES; ++i)
  {
    EXPECT_EQ(f1ptr[i], MAGIC_INT);

    const int offset = i * NUM_COMPONENTS;
    EXPECT_EQ(f2ptr[offset], MAGIC_DOUBLE);
    EXPECT_EQ(f2ptr[offset + 1], MAGIC_DOUBLE);
    EXPECT_EQ(f2ptr[offset + 2], MAGIC_DOUBLE);
  }  // END for
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
template <typename T>
void add_field_to_group(sidre::Group* gp,
                        const std::string& name,
                        int numTuples,
                        int numComponents,
                        int centering,
                        T fill_value = T())
{
  SLIC_ASSERT(centering == mint::NODE_CENTERED || centering == mint::CELL_CENTERED);
  SLIC_ASSERT(!gp->hasGroup(name));

  sidre::Group* fg = gp->createGroup(name);

  const std::string association = (centering == mint::NODE_CENTERED) ? "vertex" : "element";

  fg->createView("volume_dependent")->setString("true");
  fg->createView("association")->setString(association);
  fg->createView("topology")->setString("topo");

  sidre::View* fv = fg->createView("values");
  SLIC_ASSERT(fv != nullptr);

  sidre::deprecated::MCArray<T> data(fv, numTuples, numComponents);
  data.fill(fill_value);

  EXPECT_TRUE(gp->hasGroup(name));
}

//------------------------------------------------------------------------------
void check_blueprint(sidre::Group* field_group)
{
  EXPECT_FALSE(field_group->getName().empty());

  EXPECT_TRUE(field_group->hasChildView("association"));
  const char* assoc = field_group->getView("association")->getString();
  const bool isVertex = (strcmp(assoc, "vertex") == 0);
  const bool isElement = (isVertex) ? false : (strcmp(assoc, "element") == 0);
  EXPECT_TRUE(isVertex || isElement);

  EXPECT_TRUE(field_group->hasChildView("volume_dependent"));
  EXPECT_TRUE(field_group->hasChildView("values"));
}

#endif

}  // namespace

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

TEST(mint_mesh_field_data_DeathTest, invalid_construction)
{
  const char* IGNORE_OUTPUT = ".*";
  EXPECT_DEATH_IF_SUPPORTED(mint::FieldData(-1), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(mint::FieldData(42), IGNORE_OUTPUT);

#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* gp = ds.getRoot();

  EXPECT_DEATH_IF_SUPPORTED(mint::FieldData(42, gp, "topo"), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(mint::FieldData(mint::NODE_CENTERED, nullptr, "topo"), IGNORE_OUTPUT);

  // create sidre data that does not conform to the blue-print
  sidre::Group* f1 = gp->createGroup("f1");

  // should fail-- data does not conform to the blueprint
  EXPECT_DEATH_IF_SUPPORTED(mint::FieldData(mint::NODE_CENTERED, gp, "topo"), IGNORE_OUTPUT);

  // should still fail -- doesn't have volume_dependent view
  f1->createView("association")->setString("vertex");
  EXPECT_DEATH_IF_SUPPORTED(mint::FieldData(mint::NODE_CENTERED, gp, "topo"), IGNORE_OUTPUT);

  // should still fail -- doesn't have values view
  f1->createView("volume_dependent")->setString("true");
  EXPECT_DEATH_IF_SUPPORTED(mint::FieldData(mint::NODE_CENTERED, gp, "topo"), IGNORE_OUTPUT);

  // should still fail -- association is foo/bar
  sidre::View* fv = f1->createView("values");
  sidre::deprecated::MCArray<int> data(fv, 4, 1);
  data.fill(42);

  f1->getView("association")->setString("foobar");
  EXPECT_DEATH_IF_SUPPORTED(mint::FieldData(mint::NODE_CENTERED, gp, "topo"), IGNORE_OUTPUT);
#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_data_DeathTest, invalid_operations)
{
  const char* IGNORE_OUTPUT = ".*";

  int data[4] = {1, 2, 3, 4};

  // Create a test field_data instance
  mint::FieldData field_data(mint::NODE_CENTERED);
  field_data.createField<int>("test", data, 4);

  // creating a duplicate field (field has the same name) should fail
  EXPECT_DEATH_IF_SUPPORTED(field_data.createField<int>("test", data, 4), IGNORE_OUTPUT);

  // reserve() and resize() operations are not allowed on fields that point
  // to external buffers.
  EXPECT_DEATH_IF_SUPPORTED(field_data.reserve(10), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(field_data.resize(10), IGNORE_OUTPUT);

  // creating an external field with a null pointer should fail
  EXPECT_DEATH_IF_SUPPORTED(field_data.createField<int>("foo", nullptr, 4), IGNORE_OUTPUT);

  // add a field with an empty name should fail
  EXPECT_DEATH_IF_SUPPORTED(field_data.createField<int>("", data, 4), IGNORE_OUTPUT);

  // remove field that does not exist should fail
  EXPECT_DEATH_IF_SUPPORTED(field_data.removeField("foo"), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_data, empty_constructor)
{
  // check empty native constructor
  mint::FieldData fd(mint::NODE_CENTERED);
  check_empty_field_data<mint::NODE_CENTERED>(fd);
  EXPECT_FALSE(fd.hasSidreGroup());

#ifdef AXOM_MINT_USE_SIDRE
  // check empty constructor from a Sidre group
  sidre::DataStore ds;
  sidre::Group* gp = ds.getRoot();

  mint::FieldData fd2(mint::CELL_CENTERED, gp, "topo");
  check_empty_field_data<mint::CELL_CENTERED>(fd2);
  EXPECT_TRUE(fd2.hasSidreGroup());
#endif
}

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE
TEST(mint_mesh_field_data, sidre_constructor)
{
  constexpr int NUM_TUPLES = 10;
  constexpr int NUM_COMPONENTS = 3;
  constexpr int NTOTAL = NUM_TUPLES * NUM_COMPONENTS;

  constexpr int NUM_NODE_FIELDS = 3;
  constexpr int NUM_CELL_FIELDS = 1;
  constexpr double DOUBLE_MAGIC_NUM = 42.0;
  constexpr int INT_MAGIC_NUM = 42;

  // create data store to test with
  sidre::DataStore ds;
  sidre::Group* gp = ds.getRoot();

  // add node-centered fields
  add_field_to_group<double>(gp, "n1", NUM_TUPLES, NUM_COMPONENTS, mint::NODE_CENTERED, DOUBLE_MAGIC_NUM);
  add_field_to_group<double>(gp, "n2", NUM_TUPLES, NUM_COMPONENTS, mint::NODE_CENTERED, DOUBLE_MAGIC_NUM);
  add_field_to_group<int>(gp, "n3", NUM_TUPLES, NUM_COMPONENTS, mint::NODE_CENTERED, INT_MAGIC_NUM);

  // add cell-centered fields
  add_field_to_group<double>(gp, "c1", NUM_TUPLES, NUM_COMPONENTS, mint::CELL_CENTERED, DOUBLE_MAGIC_NUM);

  // construct node-centered and cell-centered fields
  mint::FieldData node_data(mint::NODE_CENTERED, gp, "topo");
  mint::FieldData cell_data(mint::CELL_CENTERED, gp, "topo");
  EXPECT_EQ(node_data.getNumFields(), NUM_NODE_FIELDS);
  EXPECT_EQ(cell_data.getNumFields(), NUM_CELL_FIELDS);

  // check nodal fields
  for(int i = 0; i < NUM_NODE_FIELDS; ++i)
  {
    std::ostringstream oss;
    oss << "n" << (i + 1);

    const std::string name = oss.str();
    EXPECT_TRUE(node_data.hasField(name));
    EXPECT_FALSE(cell_data.hasField(name));

    mint::Field* f = node_data.getField(i);
    EXPECT_TRUE(f != nullptr);
    EXPECT_EQ(f->getName(), name);
    EXPECT_EQ(f->getNumTuples(), NUM_TUPLES);
    EXPECT_EQ(f->getNumComponents(), NUM_COMPONENTS);

    if(f->getName() == "n3")
    {
      EXPECT_EQ(f->getType(), mint::INT32_FIELD_TYPE);

      const int* field_data = mint::Field::getDataPtr<int>(f);
      EXPECT_TRUE(field_data != nullptr);

      for(int i = 0; i < NTOTAL; ++i)
      {
        EXPECT_EQ(field_data[i], INT_MAGIC_NUM);
      }

    }  // END if
    else
    {
      EXPECT_EQ(f->getType(), mint::DOUBLE_FIELD_TYPE);

      const double* field_data = mint::Field::getDataPtr<double>(f);
      EXPECT_TRUE(field_data != nullptr);

      for(int i = 0; i < NTOTAL; ++i)
      {
        EXPECT_EQ(field_data[i], DOUBLE_MAGIC_NUM);
      }

    }  // END else

  }  // END for

  // check cell-centered fields
  for(int i = 0; i < NUM_CELL_FIELDS; ++i)
  {
    std::ostringstream oss;
    oss << "c" << (i + 1);

    const std::string name = oss.str();
    EXPECT_TRUE(cell_data.hasField(name));
    EXPECT_FALSE(node_data.hasField(name));
  }  // END for
}
#endif

//------------------------------------------------------------------------------
TEST(mint_mesh_field_data, create_and_access_fields)
{
  constexpr int NUM_TUPLES = 4;
  constexpr int NUM_COMPONENTS = 3;
  constexpr int MAGIC_INT = 42;
  constexpr double MAGIC_DOUBLE = 3.14;

  mint::FieldData native_data(mint::NODE_CENTERED);
  check_create_and_access_data(native_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(native_data.getNumFields(), 2);
  EXPECT_FALSE(native_data.hasSidreGroup());

#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* fields_group = ds.getRoot();

  int foobar[] = {7, 7, 7, 7};

  // BEGIN SCOPE
  {
    mint::FieldData sidre_data(mint::NODE_CENTERED, fields_group, "topo");
    check_create_and_access_data(sidre_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);

    EXPECT_EQ(sidre_data.getNumFields(), 2);
    EXPECT_TRUE(sidre_data.hasSidreGroup());

    // create a temp field and don't push to Sidre
    const bool putInSidre = false;
    sidre_data.createField<int>("foo", NUM_TUPLES, 1, NUM_TUPLES, putInSidre);
    EXPECT_TRUE(sidre_data.hasField("foo"));
    EXPECT_FALSE(fields_group->hasChildGroup("foo"));

    // create an external field, doesn't get pushed to Sidre
    sidre_data.createField<int>("foobar", foobar, NUM_TUPLES);
    EXPECT_TRUE(sidre_data.hasField("foobar"));
    EXPECT_FALSE(fields_group->hasChildGroup("foobar"));
  }
  // END SCOPE

  // ensures data persists in sidre after the FieldData instance is destroyed
  EXPECT_EQ(fields_group->getNumGroups(), 2);
  EXPECT_TRUE(fields_group->hasChildGroup("f1"));
  EXPECT_TRUE(fields_group->hasChildGroup("f2"));

  // ensure Sidre hierarchy conforms to the blueprint
  sidre::Group* f1 = fields_group->getGroup("f1");
  EXPECT_TRUE(f1 != nullptr);
  check_blueprint(f1);
  EXPECT_EQ(f1->getView("values")->getTypeID(), sidre::INT32_ID);

  sidre::deprecated::MCArray<int> f1_array(f1->getView("values"));
  EXPECT_EQ(f1_array.size(), NUM_TUPLES);
  EXPECT_EQ(f1_array.numComponents(), 1);

  sidre::Group* f2 = fields_group->getGroup("f2");
  EXPECT_TRUE(f2 != nullptr);
  check_blueprint(f2);
  EXPECT_EQ(f2->getView("values")->getTypeID(), sidre::FLOAT64_ID);

  sidre::deprecated::MCArray<double> f2_array(f2->getView("values"));
  EXPECT_EQ(f2_array.size(), NUM_TUPLES);
  EXPECT_EQ(f2_array.numComponents(), NUM_COMPONENTS);

  // check data
  for(int i = 0; i < NUM_TUPLES; ++i)
  {
    EXPECT_EQ(f1_array(i), MAGIC_INT);

    EXPECT_DOUBLE_EQ(f2_array(i, 0), MAGIC_DOUBLE);
    EXPECT_DOUBLE_EQ(f2_array(i, 1), MAGIC_DOUBLE);
    EXPECT_DOUBLE_EQ(f2_array(i, 2), MAGIC_DOUBLE);
  }

#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_data, remove_field)
{
  constexpr int NUM_TUPLES = 4;
  constexpr int NUM_COMPONENTS = 3;
  constexpr int MAGIC_INT = 42;
  constexpr double MAGIC_DOUBLE = 3.14;

  int foo[] = {7, 7, 7, 7};

  mint::FieldData native_data(mint::NODE_CENTERED);
  check_create_and_access_data(native_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(native_data.getNumFields(), 2);
  EXPECT_FALSE(native_data.hasSidreGroup());

  // remove f2
  native_data.removeField("f2");
  EXPECT_FALSE(native_data.hasField("f2"));
  EXPECT_EQ(native_data.getNumFields(), 1);

  // add an external field
  native_data.createField<int>("foo", foo, NUM_TUPLES);
  EXPECT_TRUE(native_data.hasField("foo"));
  EXPECT_EQ(native_data.getNumFields(), 2);
  EXPECT_EQ(native_data.getField(0)->getName(), "f1");
  EXPECT_EQ(native_data.getField(1)->getName(), "foo");

  // remove foo
  native_data.removeField("foo");
  EXPECT_FALSE(native_data.hasField("foo"));

  // ensure foo is still there
  for(int i = 0; i < NUM_TUPLES; ++i)
  {
    EXPECT_EQ(foo[i], 7);
  }

#ifdef AXOM_MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* fields_group = ds.getRoot();

  mint::FieldData sidre_data(mint::NODE_CENTERED, fields_group, "topo");
  check_create_and_access_data(sidre_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(sidre_data.getNumFields(), 2);
  EXPECT_TRUE(sidre_data.hasSidreGroup());

  // remove f2, ensure it is also remove in Sidre
  sidre_data.removeField("f2");
  EXPECT_FALSE(sidre_data.hasField("f2"));
  EXPECT_FALSE(fields_group->hasChildGroup("f2"));

#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_data, resize)
{
  constexpr int NUM_TUPLES = 4;
  constexpr int NEW_NUM_TUPLES = 10;
  constexpr int NUM_COMPONENTS = 3;
  constexpr int MAGIC_INT = 42;
  constexpr double MAGIC_DOUBLE = 3.14;

  mint::FieldData native_data(mint::NODE_CENTERED);
  check_create_and_access_data(native_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(native_data.getNumFields(), 2);
  check_resize(native_data, NEW_NUM_TUPLES);

#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* fields_group = ds.getRoot();

  mint::FieldData sidre_data(mint::NODE_CENTERED, fields_group, "topo");
  check_create_and_access_data(sidre_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(sidre_data.getNumFields(), 2);
  EXPECT_TRUE(sidre_data.hasSidreGroup());
  check_resize(sidre_data, NEW_NUM_TUPLES);

  // check the raw sidre data
  sidre::deprecated::MCArray<int> f1(fields_group->getView("f1/values"));
  sidre::deprecated::MCArray<double> f2(fields_group->getView("f2/values"));
  EXPECT_EQ(f1.size(), NEW_NUM_TUPLES);
  EXPECT_EQ(f1.numComponents(), 1);
  EXPECT_EQ(f2.size(), NEW_NUM_TUPLES);
  EXPECT_EQ(f2.numComponents(), NUM_COMPONENTS);

#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_data, reserve)
{
  constexpr int NUM_TUPLES = 4;
  constexpr int NEW_CAPACITY = 256;
  constexpr int NUM_COMPONENTS = 3;
  constexpr int MAGIC_INT = 42;
  constexpr double MAGIC_DOUBLE = 3.14;

  mint::FieldData native_data(mint::NODE_CENTERED);
  check_create_and_access_data(native_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(native_data.getNumFields(), 2);
  check_reserve(native_data, NEW_CAPACITY);

#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* fields_group = ds.getRoot();

  mint::FieldData sidre_data(mint::NODE_CENTERED, fields_group, "topo");
  check_create_and_access_data(sidre_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(sidre_data.getNumFields(), 2);
  EXPECT_TRUE(sidre_data.hasSidreGroup());
  check_reserve(sidre_data, NEW_CAPACITY);

  // check the raw sidre data
  sidre::deprecated::MCArray<int> f1(fields_group->getView("f1/values"));
  sidre::deprecated::MCArray<double> f2(fields_group->getView("f2/values"));
  EXPECT_EQ(f1.capacity(), NEW_CAPACITY);
  EXPECT_EQ(f2.capacity(), NEW_CAPACITY);
#endif
}

//------------------------------------------------------------------------------
TEST(mint_mesh_field_data, shrink)
{
  constexpr int NUM_TUPLES = 4;
  constexpr int NUM_COMPONENTS = 3;
  constexpr int MAGIC_INT = 42;
  constexpr double MAGIC_DOUBLE = 3.14;

  mint::FieldData native_data(mint::NODE_CENTERED);
  check_create_and_access_data(native_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(native_data.getNumFields(), 2);

  check_shrink(native_data);

#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* fields_group = ds.getRoot();

  mint::FieldData sidre_data(mint::NODE_CENTERED, fields_group, "topo");
  check_create_and_access_data(sidre_data, NUM_TUPLES, NUM_COMPONENTS, MAGIC_INT, MAGIC_DOUBLE);
  EXPECT_EQ(sidre_data.getNumFields(), 2);
  EXPECT_TRUE(sidre_data.hasSidreGroup());

  check_shrink(sidre_data);

  // check the raw sidre data
  sidre::deprecated::MCArray<int> f1(fields_group->getView("f1/values"));
  sidre::deprecated::MCArray<double> f2(fields_group->getView("f2/values"));
  EXPECT_EQ(f1.capacity(), f1.size());
  EXPECT_EQ(f2.capacity(), f2.size());
#endif
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
