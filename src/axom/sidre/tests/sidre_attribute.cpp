// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/sidre.hpp"

#include "axom/fmt.hpp"

#include "gtest/gtest.h"

using axom::sidre::Attribute;
using axom::sidre::CHAR8_STR_ID;
using axom::sidre::DataStore;
using axom::sidre::DOUBLE_ID;
using axom::sidre::Group;
using axom::sidre::IndexType;
using axom::sidre::INT32_ID;
using axom::sidre::INT64_ID;
using axom::sidre::INT_ID;
using axom::sidre::InvalidIndex;
using axom::sidre::Node;
using axom::sidre::View;

// Create some global attribute values, used by multiple tests
const std::string g_name_color("color");
const std::string g_color_none("white");
const std::string g_color_red("red");
const std::string g_color_blue("blue");

const std::string g_name_animal("animal");
const std::string g_animal_none("human");
const std::string g_animal_cat("cat");
const std::string g_animal_dog("dog");

const std::string g_namea("a");
const std::string g_nameb("b");

const std::string g_name_dump("dump");
const int g_dump_no = 0;
const int g_dump_yes = 1;

const std::string g_name_size("size");
const double g_size_small = 1.2;
const double g_size_medium = 2.3;
const double g_size_large = 3.4;

// intel has a problem overloading 'Attribute *' and 'IndexType'.
const Attribute* g_attr_null = nullptr;

// Test protocols
#ifdef AXOM_USE_HDF5
const int g_nprotocols = 3;
const std::string g_protocols[] = {"sidre_json", "sidre_hdf5", "json"};
#else
const int g_nprotocols = 2;
const std::string g_protocols[] = {"sidre_json", "json"};
#endif

const std::set<std::string> g_protocol_saves_attributes {"sidre_hdf5",
                                                         "sidre_conduit_json",
                                                         "sidre_json",
                                                         "sidre_layout_json"};

//------------------------------------------------------------------------------
// Create attribute in a Datastore
//
TEST(sidre_attribute, create_attr)
{
  SLIC_INFO("Some warnings are expected in the 'create_attr' test");

  bool ok;

  DataStore* ds = new DataStore();

  int nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);

  bool has_index = ds->hasAttribute(0);
  EXPECT_FALSE(has_index);
  bool has_name = ds->hasAttribute(g_name_color);
  EXPECT_FALSE(has_name);

  // Create string attribute
  Attribute* color = ds->createAttributeString(g_name_color, g_color_none);
  EXPECT_TRUE(color != nullptr);
  EXPECT_EQ(CHAR8_STR_ID, color->getTypeID());

  IndexType attr_index = color->getIndex();
  EXPECT_EQ(0, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(1, nattrs);

  has_name = ds->hasAttribute(g_name_color);
  EXPECT_TRUE(has_name);
  has_index = ds->hasAttribute(0);
  EXPECT_TRUE(has_index);

  // Try to change default to a different type.
  // Check template of setDefaultScalar.
  ok = color->setDefaultScalar(1);
  EXPECT_FALSE(ok);
  ok = color->setDefaultScalar(3.14);
  EXPECT_FALSE(ok);

  // Change to legal values.
  ok = color->setDefaultString("unknown");
  EXPECT_TRUE(ok);
#if 0
  // XXX - unsupported overload
  ok = color->setDefaultString('u');
  EXPECT_TRUE(ok);
#endif
  ok = color->setDefaultString(std::string("string"));  // non-const string
  EXPECT_TRUE(ok);

  Attribute* attr = ds->getAttribute(g_name_color);
  EXPECT_EQ(attr, color);
  const Attribute* attrc = ds->getAttribute(g_name_color);
  EXPECT_EQ(attrc, color);

  attr = ds->getAttribute(0);
  EXPECT_EQ(attr, color);
  attrc = ds->getAttribute(0);
  EXPECT_EQ(attrc, color);

  ds->destroyAttribute(color);
  nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);
  has_name = ds->hasAttribute(g_name_color);
  EXPECT_FALSE(has_name);
  // At this point color points to deallocated memory

  // Create additional attributes
  Attribute* dump = ds->createAttributeScalar(g_name_dump, g_dump_no);
  EXPECT_TRUE(dump != nullptr);

  attr_index = dump->getIndex();
  EXPECT_EQ(0, attr_index);

  Attribute* size = ds->createAttributeScalar(g_name_size, g_size_small);
  EXPECT_TRUE(size != nullptr);

  attr_index = size->getIndex();
  EXPECT_EQ(1, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(2, nattrs);

  ok = dump->setDefaultScalar(1);
  EXPECT_TRUE(ok);
  // try to change default to a different type
  ok = dump->setDefaultString(g_name_dump);
  EXPECT_FALSE(ok);
  ok = dump->setDefaultString("yes");
  EXPECT_FALSE(ok);
  ok = dump->setDefaultScalar(3.1415);
  EXPECT_FALSE(ok);

  ds->destroyAllAttributes();
  nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);

  delete ds;
}

//------------------------------------------------------------------------------
// Set attributes on a view

TEST(sidre_attribute, view_attr)
{
  SLIC_INFO("Some warnings are expected in the 'view_attr' test");

  // Note: This test relies on re-wiring conduit error handlers
  DataStore::setConduitSLICMessageHandlers();

  bool ok;

  DataStore* ds = new DataStore();

  // Create all attributes for DataStore
  Attribute* attr_color = ds->createAttributeString(g_name_color, g_color_none);
  EXPECT_TRUE(attr_color != nullptr);

  Attribute* attr_animal = ds->createAttributeString(g_name_animal, g_animal_none);
  EXPECT_TRUE(attr_animal != nullptr);

  Group* root = ds->getRoot();

  //----------------------------------------
  // Set the first attribute in a Group
  Group* grp1 = root->createGroup("grp1");
  View* view1a = grp1->createView(g_namea);
  EXPECT_TRUE(view1a != nullptr);

  EXPECT_FALSE(view1a->hasAttributeValue(g_attr_null));
  EXPECT_FALSE(view1a->hasAttributeValue(attr_color));

  // Check values of unset attributes
  const std::string out1x = view1a->getAttributeString(attr_color);
  EXPECT_EQ(g_color_none, out1x);

  const std::string out1y = view1a->getAttributeString(attr_animal);
  EXPECT_EQ(g_animal_none, out1y);

  ok = view1a->setAttributeString(attr_color, g_color_red);
  EXPECT_TRUE(ok);

  EXPECT_TRUE(view1a->hasAttributeValue(attr_color));

  const std::string out = view1a->getAttributeString(attr_color);
  EXPECT_EQ(g_color_red, out);

  // reset attribute value
  ok = view1a->setAttributeString(attr_color, g_color_blue);
  EXPECT_TRUE(ok);

  const std::string out1b = view1a->getAttributeString(attr_color);
  EXPECT_EQ(g_color_blue, out1b);

  // Check second, unset attribute. Should be default value
  EXPECT_FALSE(view1a->hasAttributeValue(attr_animal));
  const std::string out1d = view1a->getAttributeString(attr_animal);
  EXPECT_EQ(g_animal_none, out1d);

  // Now set second attribute
  ok = view1a->setAttributeString(attr_animal, g_animal_dog);
  EXPECT_TRUE(ok);

  const std::string out1c = view1a->getAttributeString(attr_animal);
  EXPECT_EQ(g_animal_dog, out1c);

  //----------------------------------------
  // Set the second attribute in a Group
  Group* grp2 = root->createGroup("grp2");

  View* view2a = grp2->createView(g_namea);
  EXPECT_TRUE(view2a != nullptr);

  EXPECT_FALSE(view2a->hasAttributeValue(attr_color));
  EXPECT_FALSE(view2a->hasAttributeValue(attr_animal));

  ok = view2a->setAttributeString(attr_animal, g_animal_dog);
  EXPECT_TRUE(ok);

  EXPECT_FALSE(view2a->hasAttributeValue(attr_color));
  EXPECT_TRUE(view2a->hasAttributeValue(attr_animal));

  const std::string out2a = view2a->getAttributeString(attr_animal);
  EXPECT_EQ(g_animal_dog, out2a);

  // Get the first, unset, attribute
  const std::string out2b = view2a->getAttributeString(attr_color);
  EXPECT_EQ(g_color_none, out2b);

  // Now set first attribute
  ok = view2a->setAttributeString(attr_color, g_color_red);
  EXPECT_TRUE(ok);

  EXPECT_TRUE(view2a->hasAttributeValue(attr_color));
  EXPECT_TRUE(view2a->hasAttributeValue(attr_animal));

  const std::string out2c = view2a->getAttributeString(attr_color);
  EXPECT_EQ(g_color_red, out2c);

  // Try to get a scalar from string
  int novalue = view2a->getAttributeScalar(attr_color);
  EXPECT_EQ(0, novalue);

  //----------------------------------------
  // Set attribute on second View in a Group
  Group* grp3 = root->createGroup("grp3");
  View* view3a = grp3->createView(g_namea);
  EXPECT_TRUE(view3a != nullptr);
  View* view3b = grp3->createView(g_nameb);
  EXPECT_TRUE(view3b != nullptr);

  ok = view3b->setAttributeString(attr_animal, g_animal_dog);
  EXPECT_TRUE(ok);

  EXPECT_FALSE(view3b->hasAttributeValue(attr_color));
  EXPECT_TRUE(view3b->hasAttributeValue(attr_animal));

  const std::string& out3a = view3b->getAttributeString(attr_animal);
  EXPECT_EQ(g_animal_dog, out3a);

  //----------------------------------------
  // Moving a view should preserve attributes
  Group* grp4 = root->createGroup("grp4");

  grp4->moveView(view3b);

  const std::string& out4a = view3b->getAttributeString(attr_animal);
  EXPECT_EQ(g_animal_dog, out4a);

  // Create an attribute which will be destroyed
  view3a->setAttributeString(attr_animal, g_animal_dog);

  grp3->destroyView(g_namea);
  grp4->destroyView(g_nameb);

  delete ds;

  // restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
// Use different type of attributes

TEST(sidre_attribute, view_int_and_double)
{
  SLIC_INFO("Some warnings are expected in the 'view_int_and_double' test");

  // Note: This test relies on re-wiring conduit error handlers
  DataStore::setConduitSLICMessageHandlers();

  bool ok;

  DataStore* ds = new DataStore();

  // Create all attributes for DataStore
  Attribute* attr_dump = ds->createAttributeScalar(g_name_dump, g_dump_no);
  EXPECT_TRUE(attr_dump != nullptr);
  EXPECT_EQ(INT_ID, attr_dump->getTypeID());

  Attribute* attr_size = ds->createAttributeScalar(g_name_size, g_size_small);
  EXPECT_TRUE(attr_size != nullptr);
  EXPECT_EQ(DOUBLE_ID, attr_size->getTypeID());

  Group* root = ds->getRoot();

  //----------------------------------------
  // Create a View
  Group* grp1 = root->createGroup("grp1");
  View* view1a = grp1->createView(g_namea);
  EXPECT_TRUE(view1a != nullptr);

  // Get default values
  int dump = view1a->getAttributeScalar(attr_dump);
  EXPECT_EQ(g_dump_no, dump);

  double size = view1a->getAttributeScalar(attr_size);
  EXPECT_EQ(g_size_small, size);

  // Set values
  ok = view1a->setAttributeScalar(attr_dump, g_dump_yes);
  EXPECT_TRUE(ok);
  dump = -1;  // clear value
  dump = view1a->getAttributeScalar(attr_dump);
  EXPECT_EQ(g_dump_yes, dump);

  ok = view1a->setAttributeScalar(attr_size, g_size_medium);
  EXPECT_TRUE(ok);
  size = 0.0;  // clear value
  size = view1a->getAttributeScalar(attr_size);
  EXPECT_EQ(g_size_medium, size);

  // Set values with incorrect types
  ok = view1a->setAttributeScalar(attr_dump, g_size_small);
  EXPECT_FALSE(ok);
  ok = view1a->setAttributeString(attr_dump, g_namea);
  EXPECT_FALSE(ok);
#if 0
  ok = view1a->setAttributeString(attr_dump, 'a');
  EXPECT_FALSE( ok );
#endif
  ok = view1a->setAttributeString(attr_dump, "g_namea");
  EXPECT_FALSE(ok);

  // Try to get a string from a scalar
  const char* nostr = view1a->getAttributeString(attr_dump);
  EXPECT_EQ(nullptr, nostr);

  int i = -1;
  i = view1a->getAttributeScalar(g_attr_null);
  EXPECT_EQ(0, i);

  delete ds;

  // restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
// Reset attribute to default

TEST(sidre_attribute, set_default)
{
  SLIC_INFO("Some warnings are expected in the 'set_default' test");

  bool ok;

  DataStore* ds = new DataStore();

  // Create all attributes for DataStore
  Attribute* attr_dump = ds->createAttributeScalar(g_name_dump, g_dump_no);
  EXPECT_TRUE(attr_dump != nullptr);
  EXPECT_EQ(INT_ID, attr_dump->getTypeID());

  Attribute* attr_size = ds->createAttributeScalar(g_name_size, g_size_small);
  EXPECT_TRUE(attr_size != nullptr);
  EXPECT_EQ(DOUBLE_ID, attr_size->getTypeID());

  Group* root = ds->getRoot();

  //----------------------------------------
  // Create a View
  Group* grp1 = root->createGroup("grp1");
  View* view1a = grp1->createView(g_namea);
  EXPECT_TRUE(view1a != nullptr);

  // reset unset attribute 1
  EXPECT_FALSE(view1a->hasAttributeValue(attr_dump));

  ok = view1a->setAttributeToDefault(attr_dump);
  EXPECT_TRUE(ok);

  EXPECT_FALSE(view1a->hasAttributeValue(attr_dump));

  // Set value
  ok = view1a->setAttributeScalar(attr_dump, g_dump_yes);
  EXPECT_TRUE(ok);
  EXPECT_TRUE(view1a->hasAttributeValue(attr_dump));

  // reset set attribute 1
  ok = view1a->setAttributeToDefault(attr_dump);
  EXPECT_TRUE(ok);
  EXPECT_FALSE(view1a->hasAttributeValue(attr_dump));

  // reset unset attribute 2
  EXPECT_FALSE(view1a->hasAttributeValue(attr_size));

  ok = view1a->setAttributeToDefault(attr_size);
  EXPECT_TRUE(ok);

  EXPECT_FALSE(view1a->hasAttributeValue(attr_size));

  // Check errors
  ok = view1a->setAttributeToDefault(g_attr_null);
  EXPECT_FALSE(ok);

  delete ds;
}

//------------------------------------------------------------------------------
// get attribute as Conduit::Node

TEST(sidre_attribute, as_node)
{
  SLIC_INFO("Some warnings are expected in the 'as_node' test");

  bool ok;

  DataStore* ds = new DataStore();

  // Create attributes for DataStore
  Attribute* attr_color = ds->createAttributeString(g_name_color, g_color_none);
  EXPECT_TRUE(attr_color != nullptr);

  Attribute* attr_dump = ds->createAttributeScalar(g_name_dump, g_dump_no);
  EXPECT_TRUE(attr_dump != nullptr);

  Group* root = ds->getRoot();

  //----------------------------------------
  // Set the first attribute in a Group
  Group* grp1 = root->createGroup("grp1");
  View* view1a = grp1->createView(g_namea);
  EXPECT_TRUE(view1a != nullptr);

  ok = view1a->setAttributeString(attr_color, g_color_red);
  EXPECT_TRUE(ok);

  const Node& node1 = view1a->getAttributeNodeRef(attr_color);
  EXPECT_EQ(g_color_red, node1.as_string());

  const Node& node2 = view1a->getAttributeNodeRef(attr_dump);
  EXPECT_EQ(g_dump_no, node2.as_int());

  const Node& node3 = view1a->getAttributeNodeRef(g_attr_null);
  EXPECT_TRUE(node3.schema().dtype().is_empty());

  delete ds;
}

//------------------------------------------------------------------------------
// Access attributes by name or index

TEST(sidre_attribute, overloads)
{
  SLIC_INFO("Some warnings are expected in the 'overloads' test");

  // Note: This test relies on re-wiring conduit error handlers
  DataStore::setConduitSLICMessageHandlers();

  bool ok;
  DataStore* ds = new DataStore();

  // Create string and scalar attributes
  Attribute* attr_color = ds->createAttributeString(g_name_color, g_color_none);
  EXPECT_TRUE(attr_color != nullptr);
  IndexType icolor = attr_color->getIndex();
  EXPECT_EQ(0, icolor);

  Attribute* attr_dump = ds->createAttributeScalar(g_name_dump, g_dump_no);
  EXPECT_TRUE(attr_dump != nullptr);
  IndexType idump = attr_dump->getIndex();
  EXPECT_EQ(1, idump);

  EXPECT_EQ(attr_color, ds->getAttribute(g_name_color));
  EXPECT_EQ(attr_color, ds->getAttribute(icolor));

  //----------------------------------------
  Group* root = ds->getRoot();
  View* view = root->createView("view1");

  // string
  ok = view->setAttributeString(attr_color, g_color_red);
  EXPECT_TRUE(ok);
  ok = view->setAttributeString(icolor, g_color_red);
  EXPECT_TRUE(ok);
  ok = view->setAttributeString(g_name_color, g_color_red);
  EXPECT_TRUE(ok);

  const char* attr1a = view->getAttributeString(attr_color);
  EXPECT_EQ(g_color_red, attr1a);
  const char* attr2a = view->getAttributeString(icolor);
  EXPECT_EQ(g_color_red, attr2a);
  const char* attr3a = view->getAttributeString(g_name_color);
  EXPECT_EQ(g_color_red, attr3a);

  // scalar
  ok = view->setAttributeScalar(attr_dump, g_dump_yes);
  EXPECT_TRUE(ok);
  ok = view->setAttributeScalar(idump, g_dump_yes);
  EXPECT_TRUE(ok);
  ok = view->setAttributeScalar(g_name_dump, g_dump_yes);
  EXPECT_TRUE(ok);

  int attr1b = view->getAttributeScalar(attr_dump);
  EXPECT_EQ(g_dump_yes, attr1b);
  int attr2b = view->getAttributeScalar(idump);
  EXPECT_EQ(g_dump_yes, attr2b);
  int attr3b = view->getAttributeScalar(g_name_dump);
  EXPECT_EQ(g_dump_yes, attr3b);

  EXPECT_EQ(g_dump_yes, view->getAttributeScalar<int>(attr_dump));
  EXPECT_EQ(g_dump_yes, view->getAttributeScalar<int>(idump));
  EXPECT_EQ(g_dump_yes, view->getAttributeScalar<int>(g_name_dump));

  const Node& node1 = view->getAttributeNodeRef(attr_dump);
  EXPECT_EQ(g_dump_yes, node1.as_int());
  const Node& node2 = view->getAttributeNodeRef(idump);
  EXPECT_EQ(g_dump_yes, node2.as_int());
  const Node& node3 = view->getAttributeNodeRef(g_name_dump);
  EXPECT_EQ(g_dump_yes, node3.as_int());

  EXPECT_TRUE(view->hasAttributeValue(attr_dump));
  EXPECT_TRUE(view->hasAttributeValue(idump));
  EXPECT_TRUE(view->hasAttributeValue(g_name_dump));

  ok = view->setAttributeToDefault(attr_dump);
  EXPECT_TRUE(ok);
  ok = view->setAttributeToDefault(idump);
  EXPECT_TRUE(ok);
  ok = view->setAttributeToDefault(g_name_dump);
  EXPECT_TRUE(ok);

  // Attribute no longer set
  EXPECT_FALSE(view->hasAttributeValue(attr_dump));
  EXPECT_FALSE(view->hasAttributeValue(idump));
  EXPECT_FALSE(view->hasAttributeValue(g_name_dump));

  // Check some errors
  EXPECT_EQ(0, view->getAttributeScalar<int>(g_attr_null));
  EXPECT_EQ(0, view->getAttributeScalar<int>(InvalidIndex));
  EXPECT_EQ(0, view->getAttributeScalar<int>("noname"));

  delete ds;

  // restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
// Test looping over Attributes and Attribute Values.

TEST(sidre_attribute, loop_attributes)
{
  DataStore* ds = new DataStore();

  // Create attributes for DataStore
  Attribute* color = ds->createAttributeString(g_name_color, g_color_none);
  IndexType icolor = color->getIndex();
  EXPECT_EQ(0, icolor);

  Attribute* dump = ds->createAttributeScalar(g_name_dump, g_dump_no);
  IndexType idump = dump->getIndex();
  EXPECT_EQ(1, idump);

  Attribute* size = ds->createAttributeScalar(g_name_size, g_size_small);
  IndexType isize = size->getIndex();
  EXPECT_EQ(2, isize);

  {
    IndexType idx1 = ds->getFirstValidAttributeIndex();
    EXPECT_EQ(0, idx1);
    IndexType idx2 = ds->getNextValidAttributeIndex(idx1);
    EXPECT_EQ(1, idx2);
    IndexType idx3 = ds->getNextValidAttributeIndex(idx2);
    EXPECT_EQ(2, idx3);
    IndexType idx4 = ds->getNextValidAttributeIndex(idx3);
    EXPECT_EQ(InvalidIndex, idx4);
    IndexType idx5 = ds->getNextValidAttributeIndex(idx4);
    EXPECT_EQ(InvalidIndex, idx5);
  }

  //----------------------------------------
  Group* root = ds->getRoot();

  // set all attributes
  View* view1 = root->createView("view1");
  view1->setAttributeString(color, g_color_red);
  view1->setAttributeScalar(dump, g_dump_yes);
  view1->setAttributeScalar(size, g_size_large);

  {
    IndexType idx1 = view1->getFirstValidAttrValueIndex();
    EXPECT_EQ(0, idx1);
    IndexType idx2 = view1->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(1, idx2);
    IndexType idx3 = view1->getNextValidAttrValueIndex(idx2);
    EXPECT_EQ(2, idx3);
    IndexType idx4 = view1->getNextValidAttrValueIndex(idx3);
    EXPECT_EQ(InvalidIndex, idx4);
  }

  // set first attribute
  View* view2 = root->createView("view2");
  view2->setAttributeString(color, g_color_red);

  {
    IndexType idx1 = view2->getFirstValidAttrValueIndex();
    EXPECT_EQ(0, idx1);
    IndexType idx2 = view2->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(InvalidIndex, idx2);
  }

  // set last attribute
  View* view3 = root->createView("view3");
  view3->setAttributeScalar(size, g_size_large);

  {
    IndexType idx1 = view3->getFirstValidAttrValueIndex();
    EXPECT_EQ(2, idx1);
    IndexType idx2 = view3->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(InvalidIndex, idx2);
  }

  // set first and last attributes
  View* view4 = root->createView("view4");
  view4->setAttributeString(color, g_color_red);
  view4->setAttributeScalar(size, g_size_large);

  {
    IndexType idx1 = view4->getFirstValidAttrValueIndex();
    EXPECT_EQ(0, idx1);
    IndexType idx2 = view4->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(2, idx2);
    IndexType idx3 = view4->getNextValidAttrValueIndex(idx2);
    EXPECT_EQ(InvalidIndex, idx3);
  }

  // no attributes
  View* view5 = root->createView("view5");

  {
    IndexType idx1 = view5->getFirstValidAttrValueIndex();
    EXPECT_EQ(InvalidIndex, idx1);
    IndexType idx2 = view5->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(InvalidIndex, idx2);
  }

  // XXX - now delete attribute and check again

  delete ds;
}

TEST(sidre_attribute, iterate_attributes)
{
  DataStore* ds = new DataStore();

  // Create attributes for DataStore
  Attribute* color = ds->createAttributeString(g_name_color, g_color_none);
  IndexType icolor = color->getIndex();
  EXPECT_EQ(0, icolor);

  Attribute* dump = ds->createAttributeScalar(g_name_dump, g_dump_no);
  IndexType idump = dump->getIndex();
  EXPECT_EQ(1, idump);

  Attribute* size = ds->createAttributeScalar(g_name_size, g_size_small);
  IndexType isize = size->getIndex();
  EXPECT_EQ(2, isize);

  for(Attribute& attr : ds->attributes())
  {
    switch(attr.getIndex())
    {
    case 0:
      EXPECT_EQ(g_name_color, attr.getName());
      break;
    case 1:
      EXPECT_EQ(g_name_dump, attr.getName());
      break;
    case 2:
      EXPECT_EQ(g_name_size, attr.getName());
      break;
    default:
      FAIL() << axom::fmt::format("Unexpected attribute: {{id:{}, name:'{}'}}",
                                  attr.getIndex(),
                                  attr.getName());
      break;
    }
  }

  for(const Attribute& attr : ds->attributes())
  {
    switch(attr.getIndex())
    {
    case 0:
      EXPECT_EQ(g_name_color, attr.getName());
      break;
    case 1:
      EXPECT_EQ(g_name_dump, attr.getName());
      break;
    case 2:
      EXPECT_EQ(g_name_size, attr.getName());
      break;
    default:
      FAIL() << axom::fmt::format("Unexpected attribute: {{id:{}, name:'{}'}}",
                                  attr.getIndex(),
                                  attr.getName());
      break;
    }
  }

  // Destroy all attributes and check for empty
  ds->destroyAllAttributes();
  EXPECT_EQ(0, ds->getNumAttributes());

  for(Attribute& attr : ds->attributes())
  {
    switch(attr.getIndex())
    {
    default:
      FAIL() << axom::fmt::format("Expected no attributes, but found: {{id:{}, name:'{}'}}",
                                  attr.getIndex(),
                                  attr.getName());
      break;
    }
  }

  delete ds;
}

//------------------------------------------------------------------------------
// save and load attributes from a file

TEST(sidre_attribute, save_attributes)
{
  //  bool ok;
  int idata[5], *bdata;

  const std::string file_path_base("sidre_attribute_datastore_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  // Create attributes for DataStore
  Attribute* color = ds1->createAttributeString(g_name_color, g_color_none);
  EXPECT_TRUE(color != nullptr);

  Attribute* dump = ds1->createAttributeScalar(g_name_dump, g_dump_no);
  EXPECT_TRUE(dump != nullptr);

  Attribute* size = ds1->createAttributeScalar(g_name_size, g_size_small);
  EXPECT_TRUE(size != nullptr);

  EXPECT_EQ(3, ds1->getNumAttributes());

  // empty
  View* view1a = root1->createView("empty");
  view1a->setAttributeString(color, "color-empty");
  view1a->setAttributeScalar(dump, g_dump_yes);
  view1a->setAttributeScalar(size, g_size_small);

  // buffer
  View* view1b = root1->createViewAndAllocate("buffer", INT_ID, 5);
  bdata = view1b->getData();
  view1b->setAttributeString(color, "color-buffer");
  view1b->setAttributeScalar(size, g_size_medium);

  // external
  View* view1c = root1->createView("external", INT_ID, 5, idata);
  view1c->setAttributeScalar(size, g_size_large);

  // scalar
  View* view1d = root1->createViewScalar("scalar", 1);
  view1d->setAttributeString(color, "color-scalar");

  // string
  View* view1e = root1->createViewString("string", "value");
  view1e->setAttributeString(color, "color-string");

  // empty without attributes
  root1->createView("empty-no-attributes");

  for(int i = 0; i < 5; i++)
  {
    idata[i] = i;
    bdata[i] = i;
  }

  //----------------------------------------

  for(int i = 0; i < g_nprotocols; ++i)
  {
    const std::string file_path = file_path_base + g_protocols[i];
    root1->save(file_path, g_protocols[i]);
  }

  delete ds1;

  //----------------------------------------
  for(int i = 0; i < g_nprotocols; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(g_protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + g_protocols[i];

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    root2->load(file_path, g_protocols[i]);
    EXPECT_EQ(3, ds2->getNumAttributes());

    // Check available attributes

    Attribute* attr_color = ds2->getAttribute(g_name_color);
    EXPECT_EQ(g_color_none, attr_color->getDefaultNodeRef().as_string());

    Attribute* attr_dump = ds2->getAttribute(g_name_dump);
    EXPECT_EQ(g_dump_no, attr_dump->getDefaultNodeRef().as_int());

    Attribute* attr_size = ds2->getAttribute(g_name_size);
    EXPECT_EQ(g_size_small, attr_size->getDefaultNodeRef().as_double());

    // Check attributes assigned to Views

    View* view2a = root2->getView("empty");
    EXPECT_TRUE(view2a->hasAttributeValue(g_name_color));
    EXPECT_TRUE(view2a->hasAttributeValue(g_name_dump));
    EXPECT_TRUE(view2a->hasAttributeValue(g_name_size));
    EXPECT_TRUE(strcmp("color-empty", view2a->getAttributeString(attr_color)) == 0);
    EXPECT_EQ(g_dump_yes, view2a->getAttributeScalar<int>(attr_dump));
    EXPECT_EQ(g_size_small, view2a->getAttributeScalar<double>(attr_size));

    View* view2b = root2->getView("buffer");
    EXPECT_TRUE(view2b->hasAttributeValue(g_name_color));
    EXPECT_FALSE(view2b->hasAttributeValue(g_name_dump));
    EXPECT_TRUE(view2b->hasAttributeValue(g_name_size));
    EXPECT_TRUE(strcmp("color-buffer", view2b->getAttributeString(attr_color)) == 0);
    EXPECT_EQ(g_size_medium, view2b->getAttributeScalar<double>(attr_size));

    View* view2c = root2->getView("external");
    EXPECT_FALSE(view2c->hasAttributeValue(g_name_color));
    EXPECT_FALSE(view2c->hasAttributeValue(g_name_dump));
    EXPECT_TRUE(view2c->hasAttributeValue(g_name_size));
    EXPECT_EQ(g_size_large, view2c->getAttributeScalar<double>(attr_size));

    View* view2d = root2->getView("scalar");
    EXPECT_TRUE(view2d->hasAttributeValue(g_name_color));
    EXPECT_FALSE(view2d->hasAttributeValue(g_name_dump));
    EXPECT_FALSE(view2d->hasAttributeValue(g_name_size));
    EXPECT_TRUE(strcmp("color-scalar", view2d->getAttributeString(attr_color)) == 0);

    View* view2e = root2->getView("string");
    EXPECT_TRUE(view2e->hasAttributeValue(g_name_color));
    EXPECT_FALSE(view2e->hasAttributeValue(g_name_dump));
    EXPECT_FALSE(view2e->hasAttributeValue(g_name_size));
    EXPECT_TRUE(strcmp("color-string", view2e->getAttributeString(attr_color)) == 0);

    View* view2f = root2->getView("empty-no-attributes");
    EXPECT_FALSE(view2f->hasAttributeValue(g_name_color));
    EXPECT_FALSE(view2f->hasAttributeValue(g_name_dump));
    EXPECT_FALSE(view2f->hasAttributeValue(g_name_size));

    delete ds2;
  }
}

//------------------------------------------------------------------------------
// save views with a specific attribute

TEST(sidre_attribute, save_by_attribute)
{
  int idata[5], jdata[5];

  const std::string file_path_base("sidre_attribute_by_attribute_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  // Create attributes for DataStore
  Attribute* dump = ds1->createAttributeScalar(g_name_dump, g_dump_no);
  EXPECT_TRUE(dump != nullptr);

  // scalar
  root1->createViewScalar("view1", 1)->setAttributeScalar(dump, g_dump_yes);

  root1->createViewScalar("view2", 2);

  // Create a deep path with and without attribute
  root1->createViewScalar("grp1a/grp1b/view3", 3);

  root1->createViewScalar("grp2a/view4", 4);  // make sure empty "views" not saved
  root1->createViewScalar("grp2a/grp2b/view5", 5)->setAttributeScalar(dump, g_dump_yes);

  root1->createView("view6", INT_ID, 5, idata)->setAttributeScalar(dump, g_dump_yes);

  // nested external view without dump, do not create intermediate Groups.
  root1->createView("grp3a/grp3b/view7", INT_ID, 5, jdata);

  for(int i = 0; i < 5; i++)
  {
    idata[i] = i;
    jdata[i] = i + 10;
  }

  //----------------------------------------

  for(int i = 0; i < g_nprotocols; ++i)
  {
    const std::string file_path = file_path_base + g_protocols[i];
    root1->save(file_path, g_protocols[i], dump);
  }

  delete ds1;

  //----------------------------------------
  for(int i = 0; i < g_nprotocols; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(g_protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + g_protocols[i];

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    root2->load(file_path, g_protocols[i]);

    // Only views with the dump attribute should exist.

    EXPECT_TRUE(root2->hasView("view1"));
    EXPECT_FALSE(root2->hasView("view2"));
    EXPECT_FALSE(root2->hasView("grp1a/grp1b/view3"));
    EXPECT_FALSE(root2->hasView("grp2a/view4"));
    EXPECT_TRUE(root2->hasView("grp2a/grp2b/view5"));
    EXPECT_TRUE(root2->hasView("view6"));
    EXPECT_FALSE(root2->hasView("grp3a/grp3b/view7"));

    delete ds2;
  }
}

TEST(sidre_attribute, save_load_group_with_attributes_new_ds)
{
  using axom::utilities::string::endsWith;

  for(const auto& protocol : g_protocols)
  {
    axom::sidre::DataStore ds1, ds2;

    const std::string ext = endsWith(protocol, "hdf5") ? "hdf5" : "json";
    const std::string filename = axom::fmt::format("saveFile_{}.{}", protocol, ext);

    SLIC_INFO(axom::fmt::format("Checking attribute save/load w/ protocol '{}' using file '{}'",
                                protocol,
                                filename));

    // set up first datastore and save to disk
    {
      // create the attributes
      ds1.createAttributeScalar("attr", 10);
      ds1.createAttributeString(g_name_color, g_color_none);

      // create groups/views and attach attributes
      auto* gr = ds1.getRoot()->createGroup("gr");
      gr->createViewScalar("scalar1", 1);
      gr->createViewScalar("scalar2", 2)->setAttributeString(g_name_color, g_color_red);
      gr->createViewScalar("scalar3", 3)->setAttributeString(g_name_color, g_color_blue);

      EXPECT_EQ(2, ds1.getNumAttributes());
      EXPECT_TRUE(INT32_ID == ds1.getAttribute("attr")->getTypeID() ||
                  INT64_ID == ds1.getAttribute("attr")->getTypeID());
      EXPECT_EQ(CHAR8_STR_ID, ds1.getAttribute(g_name_color)->getTypeID());

      EXPECT_FALSE(gr->getView("scalar1")->hasAttributeValue(g_name_color));

      EXPECT_TRUE(gr->getView("scalar2")->hasAttributeValue(g_name_color));
      EXPECT_EQ(g_color_red, gr->getView("scalar2")->getAttributeString(g_name_color));

      EXPECT_TRUE(gr->getView("scalar3")->hasAttributeValue(g_name_color));
      EXPECT_EQ(g_color_blue, gr->getView("scalar3")->getAttributeString(g_name_color));

      gr->save(filename, protocol);
    }

    if(g_protocol_saves_attributes.find(protocol) == g_protocol_saves_attributes.end())
    {
      SLIC_INFO(
        axom::fmt::format("Skipping attribute load tests for protocol '{}' -- "
                          "it doesn't support saving attributes",
                          protocol));
      continue;
    }

    // load second datastore from saved data
    {
      auto* gr = ds2.getRoot()->createGroup("gr");
      gr->load(filename, protocol);

      EXPECT_EQ(2, ds2.getNumAttributes());
      EXPECT_TRUE(INT32_ID == ds2.getAttribute("attr")->getTypeID() ||
                  INT64_ID == ds2.getAttribute("attr")->getTypeID());

      EXPECT_EQ(CHAR8_STR_ID, ds2.getAttribute(g_name_color)->getTypeID());

      EXPECT_TRUE(gr->hasView("scalar1"));
      EXPECT_FALSE(gr->getView("scalar1")->hasAttributeValue(g_name_color));

      EXPECT_TRUE(gr->hasView("scalar2"));
      EXPECT_TRUE(gr->getView("scalar2")->hasAttributeValue(g_name_color));
      EXPECT_EQ(g_color_red, gr->getView("scalar2")->getAttributeString(g_name_color));

      EXPECT_TRUE(gr->hasView("scalar3"));
      EXPECT_TRUE(gr->getView("scalar3")->hasAttributeValue(g_name_color));
      EXPECT_EQ(g_color_blue, gr->getView("scalar3")->getAttributeString(g_name_color));
    }

    EXPECT_EQ(ds1.getNumAttributes(), ds2.getNumAttributes());
    EXPECT_EQ(ds1.getAttribute(g_name_color)->getName(), ds2.getAttribute(g_name_color)->getName());
    EXPECT_EQ(ds1.getAttribute(g_name_color)->getTypeID(),
              ds2.getAttribute(g_name_color)->getTypeID());
    EXPECT_EQ(ds1.getAttribute(g_name_color)->getDefaultNodeRef().to_string(),
              ds2.getAttribute(g_name_color)->getDefaultNodeRef().to_string());
  }
}

TEST(sidre_attribute, save_load_group_with_attributes_same_ds)
{
  using axom::utilities::string::endsWith;

  for(const auto& protocol : g_protocols)
  {
    axom::sidre::DataStore ds;

    const std::string ext = endsWith(protocol, "hdf5") ? "hdf5" : "json";
    const std::string filename = axom::fmt::format("saveFile_{}.{}", protocol, ext);

    SLIC_INFO(axom::fmt::format("Checking attribute save/load w/ protocol '{}' using file '{}'",
                                protocol,
                                filename));

    // create the attributes
    ds.createAttributeScalar("attr", 10);
    ds.createAttributeString(g_name_color, g_color_none);

    // attach some attributes to views
    auto* gr = ds.getRoot()->createGroup("gr");
    gr->createViewScalar("scalar1", 1);
    gr->createViewScalar("scalar2", 2)->setAttributeString(g_name_color, g_color_red);
    gr->createViewScalar("scalar3", 3)->setAttributeString(g_name_color, g_color_blue);

    EXPECT_EQ(2, ds.getNumAttributes());

    EXPECT_FALSE(gr->getView("scalar1")->hasAttributeValue(g_name_color));

    EXPECT_TRUE(gr->getView("scalar2")->hasAttributeValue(g_name_color));
    EXPECT_EQ(g_color_red, gr->getView("scalar2")->getAttributeString(g_name_color));

    EXPECT_TRUE(gr->getView("scalar3")->hasAttributeValue(g_name_color));
    EXPECT_EQ(g_color_blue, gr->getView("scalar3")->getAttributeString(g_name_color));

    gr->save(filename, protocol);

    if(g_protocol_saves_attributes.find(protocol) == g_protocol_saves_attributes.end())
    {
      SLIC_INFO(
        axom::fmt::format("Skipping attribute load tests for protocol '{}' -- "
                          "it doesn't support saving attributes",
                          protocol));
      continue;
    }

    {
      gr->load(filename, protocol);

      // Check that things are still as expected after loading
      EXPECT_EQ(2, ds.getNumAttributes());

      EXPECT_FALSE(gr->getView("scalar1")->hasAttributeValue(g_name_color));

      EXPECT_TRUE(gr->getView("scalar2")->hasAttributeValue(g_name_color));
      EXPECT_EQ(g_color_red, gr->getView("scalar2")->getAttributeString(g_name_color));

      EXPECT_TRUE(gr->getView("scalar3")->hasAttributeValue(g_name_color));
      EXPECT_EQ(g_color_blue, gr->getView("scalar3")->getAttributeString(g_name_color));
    }

    // check that changes to attributes get overwritten when loading
    {
      // modify/remove some attributes before loading
      ds.destroyAttribute("attr");
      gr->getView("scalar1")->setAttributeString(g_name_color, g_color_red);
      gr->getView("scalar2")->setAttributeToDefault(g_name_color);
      gr->getView("scalar3")->setAttributeString(g_name_color, g_color_red);

      EXPECT_EQ(1, ds.getNumAttributes());

      // reload group from file; this should revert changes
      gr->load(filename, protocol);

      // Check that things are reverted after loading
      EXPECT_EQ(2, ds.getNumAttributes());

      EXPECT_FALSE(gr->getView("scalar1")->hasAttributeValue(g_name_color));

      EXPECT_TRUE(gr->getView("scalar2")->hasAttributeValue(g_name_color));
      EXPECT_EQ(g_color_red, gr->getView("scalar2")->getAttributeString(g_name_color));

      EXPECT_TRUE(gr->getView("scalar3")->hasAttributeValue(g_name_color));
      EXPECT_EQ(g_color_blue, gr->getView("scalar3")->getAttributeString(g_name_color));
    }
  }
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
