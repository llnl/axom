// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/sidre.hpp"
#include "axom/sidre/core/SidreTypes.hpp"

#include "gtest/gtest.h"

#include <cstring>
#include <vector>

using axom::sidre::Buffer;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::DOUBLE_ID;
using axom::sidre::FLOAT64_ID;
using axom::sidre::Group;
using axom::sidre::indexIsValid;
using axom::sidre::IndexType;
using axom::sidre::INT_ID;
using axom::sidre::InvalidIndex;
using axom::sidre::nameIsValid;
using axom::sidre::TypeID;
using axom::sidre::View;

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(sidre_group, get_name)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("test");

  EXPECT_TRUE(root->isRoot());
  EXPECT_FALSE(group->isRoot());

  EXPECT_TRUE(group->getName() == std::string("test"));

  Group* group2 = root->getGroup("foo");
  EXPECT_TRUE(group2 == nullptr);

  delete ds;
}

//------------------------------------------------------------------------------
// getPath(), getPathName()
//------------------------------------------------------------------------------
TEST(sidre_group, get_path_name)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  EXPECT_EQ(root->getParent(), root);
  const Group* croot = ds->getRoot();
  EXPECT_EQ(croot, root);
  EXPECT_EQ(root->getName(), "");
  Group* group = root->createGroup("test/a/b/c");
  Group* grp2 = root->getGroup("test/a");
  Group* grp3 = root->getGroup("test");

  EXPECT_EQ(root->getName(), std::string(""));
  EXPECT_EQ(root->getPath(), std::string(""));
  EXPECT_EQ(root->getPathName(), std::string(""));

  EXPECT_EQ(grp2->getName(), std::string("a"));
  EXPECT_EQ(grp2->getPath(), std::string("test"));
  EXPECT_EQ(grp2->getPathName(), std::string("test/a"));

  EXPECT_EQ(grp3->getName(), std::string("test"));
  EXPECT_EQ(grp3->getPath(), std::string(""));
  EXPECT_EQ(grp3->getPathName(), std::string("test"));

  EXPECT_EQ(group->getName(), std::string("c"));
  EXPECT_EQ(group->getPath(), std::string("test/a/b"));
  EXPECT_EQ(group->getPathName(), std::string("test/a/b/c"));

  delete ds;
}

//------------------------------------------------------------------------------
// createGroup(), getGroup(), hasGroup()  with path strings
//------------------------------------------------------------------------------
TEST(sidre_group, group_with_path)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Test full path access when building incrementally
  Group* group = root->createGroup("test1")->createGroup("test2")->createGroup("test3");
  Group* group2 = root->getGroup("test1/test2/test3");

  EXPECT_TRUE(nullptr != group2);
  EXPECT_EQ(group, group2);

  // Test incremental access when building full path
  Group* groupP = root->createGroup("testA/testB/testC");
  Group* groupP2 = root->getGroup("testA")->getGroup("testB")->getGroup("testC");

  EXPECT_TRUE(nullptr != groupP2);
  EXPECT_EQ(groupP, groupP2);
  // test non-const getGroup() with path
  Group* groupPParent = root->getGroup("testA/testB");
  EXPECT_EQ(groupP->getParent(), groupPParent);
  EXPECT_EQ(groupP->getParent()->getName(), "testB");

  // Now verify that code will not create missing groups.

  root->createGroup("testa")->createGroup("testb")->createGroup("testc");
  Group* group_bada = root->getGroup("BAD/testb/testc");
  Group* group_badb = root->getGroup("testa/BAD/testc");
  Group* group_badc = root->getGroup("testa/testb/BAD");

  EXPECT_EQ(group_bada, static_cast<void*>(nullptr));
  EXPECT_EQ(group_badb, static_cast<void*>(nullptr));
  EXPECT_EQ(group_badc, static_cast<void*>(nullptr));

  // Test hasGroup with paths.

  EXPECT_FALSE(root->hasGroup("BAD/testb/testc"));
  EXPECT_FALSE(root->hasGroup("testa/BAD/testc"));
  EXPECT_FALSE(root->hasGroup("testa/testb/BAD"));

  EXPECT_TRUE(root->hasGroup("test1"));
  EXPECT_TRUE(root->hasGroup("test1/test2"));
  EXPECT_TRUE(root->hasGroup("test1/test2/test3"));
  Group* group_testa = root->getGroup("testa");
  EXPECT_TRUE(group_testa->hasGroup("testb"));
  EXPECT_TRUE(group_testa->hasGroup("testb/testc"));
  EXPECT_FALSE(group_testa->hasGroup("testb/BAD"));
  EXPECT_FALSE(group_testa->hasGroup("testb/testc/BAD"));

  EXPECT_EQ(root->getNumGroups(), 3u);
  EXPECT_TRUE(root->hasGroup(0));
  EXPECT_TRUE(root->hasGroup(1));
  EXPECT_TRUE(root->hasGroup(2));
  EXPECT_FALSE(root->hasGroup(3));
  EXPECT_FALSE(root->hasGroup(InvalidIndex));

  unsigned int testbnumgroups = group_testa->getGroup("testb")->getNumGroups();
  Group* group_cdup = group_testa->createGroup("testb/testc");

  EXPECT_EQ(group_cdup, static_cast<void*>(nullptr));
  EXPECT_EQ(group_testa->getGroup("testb")->getNumGroups(), testbnumgroups);

  delete ds;
}

//------------------------------------------------------------------------------
// createGroup(), destroyGroup()  with path strings
//------------------------------------------------------------------------------
TEST(sidre_group, destroy_group_with_path)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Test full path access when building incrementally
  Group* group = root->createGroup("test1/test2/test3");
  (void)group;

  const std::size_t exp_no_groups = 0;
  const std::size_t exp_one_group = 1;

  EXPECT_EQ(exp_one_group, root->getNumGroups());
  EXPECT_EQ(exp_one_group, root->getGroup("test1")->getNumGroups());
  EXPECT_EQ(exp_one_group, root->getGroup("test1/test2")->getNumGroups());
  EXPECT_EQ(exp_no_groups, root->getGroup("test1/test2/test3")->getNumGroups());

  root->destroyGroup("test1/test2");

  EXPECT_EQ(exp_one_group, root->getNumGroups());
  EXPECT_EQ(exp_no_groups, root->getGroup("test1")->getNumGroups());
  EXPECT_FALSE(root->hasGroup("test1/test2/test3"));
  EXPECT_FALSE(root->hasGroup("test1/test2"));

  root->destroyGroup("test1/BAD");

  EXPECT_EQ(exp_one_group, root->getNumGroups());
  EXPECT_EQ(exp_no_groups, root->getGroup("test1")->getNumGroups());

  delete ds;
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(sidre_group, get_parent)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* parent = root->createGroup("parent");
  Group* child = parent->createGroup("child");

  EXPECT_TRUE(child->getParent() == parent);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(sidre_group, get_datastore)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("parent");

  EXPECT_TRUE(group->getDataStore() == ds);

  DataStore const* const_ds = group->getDataStore();
  EXPECT_TRUE(const_ds == ds);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getGroup()
//------------------------------------------------------------------------------
TEST(sidre_group, get_group)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  Group* child = parent->createGroup("child");
  EXPECT_TRUE(child->getParent() == parent);

  EXPECT_TRUE(parent->getGroup("child") == child);
  EXPECT_TRUE(parent->getGroup(0) == child);

  // check error condition
  EXPECT_TRUE(parent->getGroup("non-existant group") == nullptr);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getView()
//------------------------------------------------------------------------------
TEST(sidre_group, get_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  View* view = parent->createView("view");

  EXPECT_TRUE(parent->getView("view") == view);
  EXPECT_TRUE(parent->getView(0) == view);

  // check error condition
  EXPECT_TRUE(parent->getView("non-existant view") == nullptr);

  delete ds;
}

//------------------------------------------------------------------------------
// createView, hasView(), getView(), destroyView() with path strings
//------------------------------------------------------------------------------
TEST(sidre_group, view_with_path)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Test with full path access when building incrementally
  View* view = root->createGroup("group1")->createGroup("group2")->createView("view1");
  View* view2 = root->getView("group1/group2/view1");

  EXPECT_TRUE(nullptr != view2);
  EXPECT_EQ(view, view2);

  // Test incremental access when building with full path
  View* viewP = root->createView("groupA/groupB/viewA");
  View* viewP2 = root->getGroup("groupA")->getGroup("groupB")->getView("viewA");

  EXPECT_TRUE(nullptr != viewP2);
  EXPECT_EQ(viewP, viewP2);

  // Now verify that bad paths just return null, and don't create missing groups
  View* v_bad1 = root->getView("BAD/groupB/viewA");
  View* v_bad2 = root->getView("groupA/BAD/viewA");
  View* v_bad3 = root->getView("groupA/groupB/BAD");

  EXPECT_EQ(v_bad1, static_cast<void*>(nullptr));
  EXPECT_EQ(v_bad2, static_cast<void*>(nullptr));
  EXPECT_EQ(v_bad3, static_cast<void*>(nullptr));

  const std::size_t exp_no_groups = 0;
  const std::size_t exp_one_group = 1;
  const std::size_t exp_two_group = 2;

  EXPECT_EQ(exp_two_group, root->getNumGroups());
  EXPECT_TRUE(root->hasGroup("group1"));
  EXPECT_TRUE(root->hasGroup("groupA"));
  EXPECT_EQ(exp_one_group, root->getGroup("group1")->getNumGroups());
  EXPECT_TRUE(root->hasGroup("group1/group2"));
  EXPECT_EQ(exp_no_groups, root->getGroup("group1/group2")->getNumGroups());
  EXPECT_EQ(exp_one_group, root->getGroup("group1/group2")->getNumViews());
  EXPECT_EQ(root->getGroup("group1/group2")->getView("view1"), view);
  EXPECT_EQ(root->getGroup("group1")->getView("group2/view1"), view);

  EXPECT_EQ(exp_one_group, root->getGroup("groupA")->getNumGroups());
  EXPECT_TRUE(root->hasGroup("groupA/groupB"));
  EXPECT_EQ(exp_no_groups, root->getGroup("groupA/groupB")->getNumGroups());
  EXPECT_EQ(exp_one_group, root->getGroup("groupA/groupB")->getNumViews());
  EXPECT_EQ(root->getGroup("groupA/groupB")->getView("viewA"), viewP);
  EXPECT_EQ(root->getGroup("groupA")->getView("groupB/viewA"), viewP);

  root->destroyView("group1/group2/view1");

  EXPECT_EQ(exp_no_groups, root->getGroup("group1/group2")->getNumViews());
  EXPECT_FALSE(root->getGroup("group1/group2")->hasView("view1"));
  EXPECT_EQ(root->getGroup("group1/group2")->getView("view1"), static_cast<void*>(nullptr));
  EXPECT_FALSE(root->hasView("group1/group2/view1"));
  EXPECT_EQ(root->getView("group1/group2/view1"), static_cast<void*>(nullptr));

  Group* groupA = root->getGroup("groupA");
  EXPECT_TRUE(groupA->hasView("groupB/viewA"));
  EXPECT_EQ(groupA->getView("groupB/viewA"), viewP);
  EXPECT_TRUE(root->hasView("groupA/groupB/viewA"));
  EXPECT_EQ(root->getView("groupA/groupB/viewA"), viewP);

  groupA->destroyView("groupB/viewA");

  EXPECT_EQ(exp_no_groups, groupA->getGroup("groupB")->getNumViews());
  EXPECT_FALSE(groupA->getGroup("groupB")->hasView("viewA"));
  EXPECT_EQ(groupA->getGroup("groupB")->getView("viewA"), static_cast<void*>(nullptr));
  EXPECT_FALSE(groupA->hasView("groupB/viewA"));
  EXPECT_EQ(groupA->getView("groupB/viewA"), static_cast<void*>(nullptr));
  EXPECT_EQ(root->getView("groupA/groupB/viewA"), static_cast<void*>(nullptr));

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(sidre_group, get_view_names_and_indicies)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  View* view1 = parent->createView("view1");
  View* view2 = parent->createView("view2");

  EXPECT_EQ(parent->getNumViews(), 2u);

  IndexType idx1 = parent->getViewIndex("view1");
  IndexType idx2 = parent->getViewIndex("view2");
  EXPECT_EQ(idx1, view1->getIndex());
  EXPECT_EQ(idx2, view2->getIndex());

  const std::string& name1 = parent->getViewName(idx1);
  const std::string& name2 = parent->getViewName(idx2);

  EXPECT_EQ(name1, std::string("view1"));
  EXPECT_EQ(view1->getName(), name1);

  EXPECT_EQ(name2, std::string("view2"));
  EXPECT_EQ(view2->getName(), name2);

  // check error conditions
  IndexType idx3 = parent->getViewIndex("view3");
  EXPECT_FALSE(indexIsValid(idx3));

  const std::string& name3 = parent->getViewName(idx3);
  EXPECT_TRUE(name3.empty());
  EXPECT_FALSE(nameIsValid(name3));

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getFirstValidGroupIndex, getNextValidGroupIndex
//------------------------------------------------------------------------------
TEST(sidre_group, get_first_and_next_group_index)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  Group* group1 = parent->createGroup("group1");
  Group* group2 = parent->createGroup("group2");
  EXPECT_EQ(parent->getNumGroups(), 2u);

  IndexType idx1 = parent->getFirstValidGroupIndex();
  IndexType idx2 = parent->getNextValidGroupIndex(idx1);
  IndexType idx3 = parent->getNextValidGroupIndex(idx2);
  EXPECT_EQ(0, idx1);
  EXPECT_EQ(1, idx2);
  EXPECT_EQ(InvalidIndex, idx3);

  Group* group1out = parent->getGroup(idx1);
  Group* group2out = parent->getGroup(idx2);
  EXPECT_EQ(group1, group1out);
  EXPECT_EQ(group2, group2out);

  // check error conditions
  Group* emptyGroup = root->createGroup("emptyGroup");
  IndexType badidx1 = emptyGroup->getFirstValidGroupIndex();
  IndexType badidx2 = emptyGroup->getNextValidGroupIndex(badidx1);

  EXPECT_FALSE(indexIsValid(badidx1));
  EXPECT_FALSE(indexIsValid(badidx2));

  delete ds;
}

//------------------------------------------------------------------------------
// Verify Groups holding items in the list format
//------------------------------------------------------------------------------
TEST(sidre_group, child_lists)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // parent is a Group in list format.
  Group* parent = root->createGroup("parent", true);

  // Create 10 unnamed Groups as children of parent.
  for(IndexType i = 0; i < 10; ++i)
  {
    Group* unnamed_group = parent->createUnnamedGroup();
    unnamed_group->createViewScalar("val", i);
  }

  // Create 15 unnamed Views as children of parent.
  for(IndexType i = 0; i < 15; ++i)
  {
    View* unnamed_view;
    if(i % 3 == 0)
    {
      unnamed_view = parent->createView("");
    }
    else if(i % 3 == 1)
    {
      unnamed_view = parent->createViewScalar("", i * i);
    }
    else
    {
      unnamed_view = parent->createViewString("", "foo");
    }
    if(!unnamed_view->isApplied())
    {
      unnamed_view->apply(INT_ID, i);
      unnamed_view->allocate(INT_ID, i);
      int* vdata = unnamed_view->getData();
      for(IndexType j = 0; j < i; ++j)
      {
        vdata[j] = j + 3;
      }
    }
  }

  // Create Group not in list format, show that it can't create
  // unnamed children.
  Group* not_list = root->createGroup("not_list", false);
  Group* dummy_group = not_list->createUnnamedGroup();
  View* dummy_view = not_list->createView("");
  EXPECT_TRUE(not_list->isUsingMap());
  EXPECT_EQ(not_list->getNumGroups(), 0);
  EXPECT_EQ(not_list->getNumViews(), 0);
  EXPECT_TRUE(dummy_group == nullptr);
  EXPECT_TRUE(dummy_view == nullptr);

  // Access data from unnamed Groups held by parent.
  std::set<int> scalars;
  for(IndexType idx = parent->getFirstValidGroupIndex(); indexIsValid(idx);
      idx = parent->getNextValidGroupIndex(idx))
  {
    Group* unnamed_group = parent->getGroup(idx);
    View* val_view = unnamed_group->getView("val");
    IndexType val = val_view->getScalar();
    EXPECT_TRUE(val >= 0 && val < 10);

    scalars.insert(val);
  }

  EXPECT_TRUE(parent->isUsingList());
  EXPECT_TRUE(scalars.size() == 10);
  EXPECT_TRUE(parent->hasGroup(6));
  EXPECT_FALSE(parent->hasGroup(20));

  // Destroy five of the unnamed Groups held by parent.
  for(IndexType idx = parent->getFirstValidGroupIndex(); indexIsValid(idx);
      idx = parent->getNextValidGroupIndex(idx))
  {
    if(idx % 2 == 1)
    {
      parent->destroyGroup(idx);
    }
  }

  // Add one more unnamed Group, so there should be six child Groups.
  (void)parent->createUnnamedGroup();
  EXPECT_EQ(parent->getNumGroups(), 6);

  // Access data from the unnamed Views.
  for(IndexType idx = parent->getFirstValidViewIndex(); indexIsValid(idx);
      idx = parent->getNextValidViewIndex(idx))
  {
    View* unnamed_view = parent->getView(idx);
    if(idx % 3 == 0)
    {
      EXPECT_EQ(unnamed_view->getTypeID(), INT_ID);
      IndexType num_elems = unnamed_view->getNumElements();
      EXPECT_EQ(num_elems, idx);
      int* vdata = unnamed_view->getData();
      for(IndexType j = 0; j < num_elems; ++j)
      {
        EXPECT_EQ(vdata[j], j + 3);
      }
    }
    else if(idx % 3 == 1)
    {
      EXPECT_TRUE(unnamed_view->isScalar());
      IndexType val = unnamed_view->getScalar();
      EXPECT_EQ(val, idx * idx);
    }
    else
    {
      EXPECT_TRUE(unnamed_view->isString());
      std::string vstr = unnamed_view->getString();
      EXPECT_EQ(vstr, std::string("foo"));
    }
  }

  root->destroyGroup("parent");

  delete ds;
}

//------------------------------------------------------------------------------
// Verify results with various path arguments for items in list
//------------------------------------------------------------------------------
TEST(sidre_group, list_item_names)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Create a group that uses the list format.
  Group* list_test = root->createGroup("list_test", true);

  // It is recommended that all items held in a Group that uses the list
  // format be unnamed, as the names are not useful to access those items
  // from the parent Group.  Nonetheless it is allowed, and this test
  // verifies that the names are created as expected.

  // Test a group created with createUnnamedGroup, the recommended way to
  // create groups for a list.
  Group* unnamed_group = list_test->createUnnamedGroup();

  ASSERT_TRUE(unnamed_group->getName().empty());

  // Test groups created with empty string, a simple name, and a path.

  // The API says that groups cannot be created with an empty string argument,
  // so a nullptr is returned for this one.
  Group* blank_group = list_test->createGroup("");

  // A simple name with no path syntax will be assigned to the group.
  Group* named_group = list_test->createGroup("named");

  // With a path, the leading names are ignored and the final name is used.
  Group* path_group_a = list_test->createGroup("testing/path");
  Group* path_group_b = list_test->createGroup("testing/longer/path");
  Group* found_group_a = list_test->getGroup("testing/path");
  Group* found_group_b = list_test->getGroup("testing/longer/path");

  ASSERT_EQ(blank_group, nullptr);
  ASSERT_EQ(named_group->getName(), "named");
  ASSERT_EQ(path_group_a, nullptr);
  ASSERT_EQ(path_group_b, nullptr);
  ASSERT_EQ(found_group_a, nullptr);
  ASSERT_EQ(found_group_b, nullptr);

  // Similar tests for views

  // An empty string is the recommended way to create an unnamed view for
  // a list.
  View* blank_view = list_test->createView("");

  // A simple name with no path syntax will be assigned to the view.
  View* named_view = list_test->createView("named");

  // With a path, the leading names are ignored and the final name is used.
  View* path_view_a = list_test->createView("testing/path");
  View* path_view_b = list_test->createView("testing/longer/path");
  View* found_view_a = list_test->getView("testing/path");
  View* found_view_b = list_test->getView("testing/longer/path");

  ASSERT_TRUE(blank_view->getName().empty());
  ASSERT_EQ(named_view->getName(), "named");
  ASSERT_EQ(path_view_a, nullptr);
  ASSERT_EQ(path_view_b, nullptr);
  ASSERT_EQ(found_view_a, nullptr);
  ASSERT_EQ(found_view_b, nullptr);

  root->destroyGroup("list_test");

  delete ds;
}
//------------------------------------------------------------------------------
// Test Group's list format for holding contents of a vector of strings.
//------------------------------------------------------------------------------
TEST(sidre_group, string_list)
{
  // Round-trip test from std::vector<std::string> to Group and back.
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  std::vector<std::string> str_vec;
  str_vec.push_back("This");
  str_vec.push_back("is");
  str_vec.push_back("a");
  str_vec.push_back("vector");
  str_vec.push_back("to");
  str_vec.push_back("test");
  str_vec.push_back("strings");
  str_vec.push_back("in");
  str_vec.push_back("sidre::Group's");
  str_vec.push_back("list");
  str_vec.push_back("format");

  // my_strings is a Group in list format.
  const bool use_list_collection = true;
  Group* my_strings = root->createGroup("my_strings", use_list_collection);

  // Put strings into the Group.
  for(auto itr = str_vec.begin(); itr != str_vec.end(); ++itr)
  {
    // The first parameter will be ignored when creating a View in a Group
    // that uses list collections, so we use the empty string
    View* str_view = my_strings->createViewString("", *itr);
    EXPECT_FALSE(str_view == nullptr);
    EXPECT_TRUE(str_view->isString());
  }

  std::vector<std::string> test_vec;

  // Get strings from the Group.
  for(IndexType idx = my_strings->getFirstValidViewIndex(); indexIsValid(idx);
      idx = my_strings->getNextValidViewIndex(idx))
  {
    View* str_view = my_strings->getView(idx);
    EXPECT_FALSE(str_view == nullptr);
    EXPECT_TRUE(str_view->isString());
    std::string vstr = str_view->getString();
    test_vec.push_back(vstr);
  }

  EXPECT_TRUE(str_vec == test_vec);

  delete ds;
}

//------------------------------------------------------------------------------
// Iterate Views with getFirstValidViewIndex, getNextValidViewIndex
//------------------------------------------------------------------------------
TEST(sidre_group, iterate_groups)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  (void)parent->createView("view1");
  (void)parent->createView("view2");
  (void)parent->createView("view3");
  (void)parent->createGroup("g1");
  (void)parent->createView("view4");
  (void)parent->createView("view5");
  (void)parent->createView("view6");
  (void)parent->createGroup("g2");
  (void)parent->createGroup("g3");
  (void)parent->createView("view7");
  (void)parent->createView("view8");
  (void)parent->createView("view9");
  (void)parent->createGroup("g4");

  int groupcount = 0;
  for(IndexType idx = parent->getFirstValidGroupIndex(); indexIsValid(idx);
      idx = parent->getNextValidGroupIndex(idx))
  {
    groupcount += 1;
  }
  EXPECT_EQ(4, groupcount);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, iterate_groups_with_iterator)
{
  using axom::utilities::string::endsWith;

  axom::sidre::DataStore ds;
  auto* foo_group = ds.getRoot()->createGroup("foo");
  foo_group->createGroup("bar_group");
  foo_group->createGroup("bar_group/child_1");
  foo_group->createGroup("bar_group/child_2");

  foo_group->createGroup("baz_group");
  foo_group->createGroup("baz_group/child_1");
  foo_group->createGroup("baz_group/child_2");
  foo_group->createGroup("baz_group/child_3");

  foo_group->createGroup("qux_group");
  foo_group->createGroup("qux_group/child_1");

  foo_group->createView("bar_view");
  foo_group->createView("baz_view");
  foo_group->createView("qux_view");
  foo_group->createView("quux_view");

  constexpr int numExpGroups = 3;
  constexpr int numExpViews = 4;

  // iterate through groups and views of 'foo' using range-for syntax
  {
    int nGroups = 0;
    for(auto& group : foo_group->groups())
    {
      EXPECT_EQ(foo_group, group.getParent());
      EXPECT_TRUE(endsWith(group.getName(), "_group"));
      ++nGroups;
    }
    EXPECT_EQ(numExpGroups, nGroups);

    int nViews = 0;
    for(auto& view : foo_group->views())
    {
      EXPECT_EQ(foo_group, view.getOwningGroup());
      EXPECT_TRUE(endsWith(view.getName(), "_view"));
      ++nViews;
    }
    EXPECT_EQ(numExpViews, nViews);
  }

  // const iterate through groups and views of 'foo' using range-for syntax
  // starting from a non-const group
  {
    int nGroups = 0;
    for(const auto& group : foo_group->groups())
    {
      EXPECT_EQ(foo_group, group.getParent());
      EXPECT_TRUE(endsWith(group.getName(), "_group"));
      ++nGroups;
    }
    EXPECT_EQ(numExpGroups, nGroups);

    int nViews = 0;
    for(const auto& view : foo_group->views())
    {
      EXPECT_EQ(foo_group, view.getOwningGroup());
      EXPECT_TRUE(endsWith(view.getName(), "_view"));
      ++nViews;
    }
    EXPECT_EQ(numExpViews, nViews);
  }

  // iterate through groups and views of 'foo' using range-for syntax
  // starting from a const group
  {
    const auto* cfoo_group = ds.getRoot()->getGroup("foo");

    int nGroups = 0;
    for(const auto& group : cfoo_group->groups())
    {
      EXPECT_EQ(foo_group, group.getParent());
      EXPECT_TRUE(endsWith(group.getName(), "_group"));
      ++nGroups;
    }
    EXPECT_EQ(numExpGroups, nGroups);

    int nViews = 0;
    for(const auto& view : foo_group->views())
    {
      EXPECT_EQ(foo_group, view.getOwningGroup());
      EXPECT_TRUE(endsWith(view.getName(), "_view"));
      ++nViews;
    }
    EXPECT_EQ(numExpViews, nViews);
  }

  // iterate though groups and views of 'foo' using iterator syntax
  {
    int nGroups = 0;
    for(auto it = foo_group->groups().begin(), itEnd = foo_group->groups().end(); it != itEnd; ++it)
    {
      auto& grp = *it;
      EXPECT_EQ(grp.getPath(), it->getPath());

      EXPECT_EQ(foo_group, it->getParent());
      EXPECT_TRUE(endsWith(it->getName(), "_group"));
      ++nGroups;
    }
    EXPECT_EQ(numExpGroups, nGroups);

    int nViews = 0;
    for(auto it = foo_group->views().begin(), itEnd = foo_group->views().end(); it != itEnd; ++it)
    {
      auto& view = *it;
      EXPECT_EQ(view.getPath(), it->getPath());

      EXPECT_EQ(foo_group, it->getOwningGroup());
      EXPECT_TRUE(endsWith(it->getName(), "_view"));
      ++nViews;
    }
    EXPECT_EQ(numExpViews, nViews);
  }

  // iterate though groups and views of 'foo' using const iterator syntax
  {
    int nGroups = 0;
    for(auto it = foo_group->groups().cbegin(), itEnd = foo_group->groups().cend(); it != itEnd; ++it)
    {
      EXPECT_EQ(foo_group, it->getParent());
      EXPECT_TRUE(endsWith(it->getName(), "_group"));
      ++nGroups;
    }
    EXPECT_EQ(numExpGroups, nGroups);

    int nViews = 0;
    for(auto it = foo_group->views().cbegin(), itEnd = foo_group->views().cend(); it != itEnd; ++it)
    {
      EXPECT_EQ(foo_group, it->getOwningGroup());
      EXPECT_TRUE(endsWith(it->getName(), "_view"));
      ++nViews;
    }
    EXPECT_EQ(numExpViews, nViews);
  }
}

//------------------------------------------------------------------------------
// Verify getFirstValidViewIndex, getNextValidViewIndex
//------------------------------------------------------------------------------
TEST(sidre_group, get_first_and_next_view_index)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  View* view1 = parent->createView("view1");
  View* view2 = parent->createView("view2");
  EXPECT_EQ(parent->getNumViews(), 2u);

  IndexType idx1 = parent->getFirstValidViewIndex();
  IndexType idx2 = parent->getNextValidViewIndex(idx1);
  IndexType idx3 = parent->getNextValidViewIndex(idx2);
  EXPECT_EQ(0, idx1);
  EXPECT_EQ(1, idx2);
  EXPECT_EQ(InvalidIndex, idx3);

  View* view1out = parent->getView(idx1);
  View* view2out = parent->getView(idx2);
  EXPECT_EQ(view1, view1out);
  EXPECT_EQ(view2, view2out);

  // check error conditions
  Group* emptyGroup = root->createGroup("emptyGroup");
  IndexType badidx1 = emptyGroup->getFirstValidViewIndex();
  IndexType badidx2 = emptyGroup->getNextValidViewIndex(badidx1);

  EXPECT_FALSE(indexIsValid(badidx1));
  EXPECT_FALSE(indexIsValid(badidx2));

  delete ds;
}

//------------------------------------------------------------------------------
// Iterate Views with getFirstValidViewIndex, getNextValidViewIndex
//------------------------------------------------------------------------------
TEST(sidre_group, iterate_views)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  (void)parent->createView("view1");
  (void)parent->createView("view2");
  (void)parent->createView("view3");
  (void)parent->createGroup("g1");
  (void)parent->createView("view4");
  (void)parent->createView("view5");
  (void)parent->createView("view6");
  (void)parent->createGroup("g2");
  (void)parent->createGroup("g3");
  (void)parent->createView("view7");
  (void)parent->createView("view8");
  (void)parent->createView("view9");
  (void)parent->createGroup("g4");

  int viewcount = 0;
  IndexType idx = parent->getFirstValidViewIndex();
  while(indexIsValid(idx))
  {
    viewcount += 1;

    idx = parent->getNextValidViewIndex(idx);
  }
  EXPECT_EQ(9, viewcount);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getGroupName(), getGroupIndex()
//------------------------------------------------------------------------------
TEST(sidre_group, get_group_name_index)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  Group* group1 = parent->createGroup("group1");
  Group* group2 = parent->createGroup("group2");

  EXPECT_EQ(parent->getNumGroups(), 2u);

  IndexType idx1 = parent->getGroupIndex("group1");
  IndexType idx2 = parent->getGroupIndex("group2");
  EXPECT_EQ(idx1, group1->getIndex());
  EXPECT_EQ(idx2, group2->getIndex());

  const std::string& name1 = parent->getGroupName(idx1);
  const std::string& name2 = parent->getGroupName(idx2);

  EXPECT_EQ(name1, std::string("group1"));
  EXPECT_EQ(group1->getName(), name1);

  EXPECT_EQ(name2, std::string("group2"));
  EXPECT_EQ(group2->getName(), name2);

  // check error conditions
  IndexType idx3 = parent->getGroupIndex("group3");
  EXPECT_FALSE(indexIsValid(idx3));

  const std::string& name3 = parent->getGroupName(idx3);
  EXPECT_TRUE(name3.empty());
  EXPECT_FALSE(nameIsValid(name3));

  delete ds;
}

//------------------------------------------------------------------------------
// createView()
// createViewAndAllocate()
// destroyView()
// destroyViewAndData()
// hasView()
//------------------------------------------------------------------------------
TEST(sidre_group, create_destroy_has_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("parent");

  View* view = group->createView("view");
  EXPECT_TRUE(group->getParent() == root);
  EXPECT_FALSE(view->hasBuffer());
  EXPECT_TRUE(group->hasView("view"));
  IndexType iview = group->getViewIndex("view");
  EXPECT_EQ(0, iview);
  iview = view->getIndex();
  EXPECT_EQ(0, iview);

  // try creating view again, should be a no-op.
  EXPECT_TRUE(group->createView("view") == nullptr);

  // Create another view to make sure destroyView only destroys one view
  group->createView("viewfiller");
  EXPECT_EQ(2u, group->getNumViews());
  IndexType iviewfiller = group->getViewIndex("viewfiller");
  EXPECT_EQ(1, iviewfiller);

  group->destroyView("view");
  EXPECT_EQ(1u, group->getNumViews());
  // Check if index changed
  EXPECT_EQ(iviewfiller, group->getViewIndex("viewfiller"));

  // destroy already destroyed group.  Should be a no-op, not a failure
  group->destroyView("view");
  EXPECT_EQ(1u, group->getNumViews());
  EXPECT_FALSE(group->hasView("view"));

  // try api call that specifies specific type and length
  group->createViewAndAllocate("viewWithLength1", INT_ID, 50);
  IndexType iview2 = group->getViewIndex("viewWithLength1");
  EXPECT_EQ(iview, iview2);  // reuse slot

  // error condition check - try again with duplicate name, should be a no-op
  EXPECT_TRUE(group->createViewAndAllocate("viewWithLength1", FLOAT64_ID, 50) == nullptr);
  group->destroyViewAndData("viewWithLength1");
  EXPECT_FALSE(group->hasView("viewWithLength1"));

  EXPECT_TRUE(group->createViewAndAllocate("viewWithLengthBadLen", FLOAT64_ID, -1) == nullptr);

  // try api call that specifies data type in another way
  group->createViewAndAllocate("viewWithLength2", DataType::float64(50));
  EXPECT_TRUE(group->createViewAndAllocate("viewWithLength2", DataType::float64(50)) == nullptr);
  // destroy view and its buffer using index
  IndexType indx = group->getFirstValidViewIndex();
  IndexType bindx = group->getView(indx)->getBuffer()->getIndex();
  group->destroyViewAndData(indx);
  EXPECT_TRUE(ds->getBuffer(bindx) == nullptr);

  // Destroy view but not the buffer
  view = group->createViewAndAllocate("viewWithLength2", INT_ID, 50);
  Buffer* buff = view->getBuffer();
  group->destroyView("viewWithLength2");
  EXPECT_TRUE(buff->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------
// createViewAndAllocate() with zero-sized array
//------------------------------------------------------------------------------
TEST(sidre_group, create_zero_sized_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  View* zeroSizedView = root->createViewAndAllocate("foo", INT_ID, 0);
  EXPECT_TRUE(zeroSizedView->isDescribed());
  EXPECT_TRUE(zeroSizedView->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------
// createGroup()
// destroyGroup()
// hasGroup()
//------------------------------------------------------------------------------
TEST(sidre_group, create_destroy_has_group)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("group");
  EXPECT_TRUE(group->getParent() == root);

  EXPECT_TRUE(root->hasGroup("group"));

  root->destroyGroup("group");
  EXPECT_FALSE(root->hasGroup("group"));

  // should be a no-op, not a failure
  root->destroyGroup("group");

  Group* group2 = root->createGroup("group2");
  // shut up compiler about unused variable
  (void)group2;
  root->destroyGroup(root->getFirstValidGroupIndex());

  delete ds;
}

TEST(sidre_group, destroy_group_and_data)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group0 = root->createGroup("group0");
  Group* group1 = root->createGroup("group1");
  Group* child0 = group0->createGroup("child0");
  Group* child1 = group0->createGroup("child1");
  Group* child2 = group1->createGroup("child2");
  Group* child3 = group1->createGroup("child3");
  Group* child4 = group1->createGroup("child4");

  child0->createViewAndAllocate("intview", INT_ID, 15);
  Group* foo0 = child0->createGroup("foo");
  child0->createGroup("empty");
  child0->createViewScalar("sclview", 3.14159);
  child0->createViewString("strview", "Hello world.");
  foo0->createViewAndAllocate("fooview", FLOAT64_ID, 12);

  int* int0_vals = child0->getView("intview")->getData();
  for(int i = 0; i < 15; ++i)
  {
    int0_vals[i] = i;
  }

  conduit::float64* flt0_vals = foo0->getView("fooview")->getData();
  for(int i = 0; i < 12; ++i)
  {
    flt0_vals[i] = (conduit::float64)(-i);
  }

  Buffer* intbuf = child0->getView("intview")->getBuffer();
  Buffer* fltbuf = foo0->getView("fooview")->getBuffer();

  //Store each Buffer's index for later testing.
  IndexType int_idx = intbuf->getIndex();
  IndexType flt_idx = fltbuf->getIndex();

  child1->createView("intview", INT_ID, 15, intbuf);
  Group* foo1 = child1->createGroup("foo");
  child1->createGroup("empty");
  child1->createViewScalar("sclview", 3.14159);
  child1->createViewString("strview", "Hello world.");
  foo1->createView("fooview", FLOAT64_ID, 12, fltbuf);

  child2->createView("intview", INT_ID, 15, intbuf);
  Group* foo2 = child2->createGroup("foo");
  child2->createGroup("empty");
  child2->createViewScalar("sclview", 3.14159);
  child2->createViewString("strview", "Hello world.");
  foo2->createView("fooview", FLOAT64_ID, 12, fltbuf);

  child3->createView("intview", INT_ID, 15, intbuf);
  Group* foo3 = child3->createGroup("foo");
  child3->createGroup("empty");
  child3->createViewScalar("sclview", 3.14159);
  child3->createViewString("strview", "Hello world.");
  foo3->createView("fooview", FLOAT64_ID, 12, fltbuf);

  child4->createView("intview", INT_ID, 15, intbuf);
  Group* foo4 = child4->createGroup("foo");
  child4->createGroup("empty");
  child4->createViewScalar("sclview", 3.14159);
  child4->createViewString("strview", "Hello world.");
  foo4->createView("fooview", FLOAT64_ID, 12, fltbuf);

  //Beginning state:  there are 2 Buffers, each attached to 5 Views.
  EXPECT_EQ(ds->getNumBuffers(), 2);
  EXPECT_EQ(intbuf->getNumViews(), 5);
  EXPECT_EQ(fltbuf->getNumViews(), 5);

  //Destroy "child0/foo" by path. This destroys the View that created fltbuf.
  EXPECT_EQ(child0->getNumGroups(), 2);
  group0->destroyGroupAndData("child0/foo");

  //Verify that fltbuf now is attached to 4 Views and child0 has only one
  //group "empty".
  EXPECT_EQ(fltbuf->getNumViews(), 4);
  EXPECT_EQ(child0->getNumGroups(), 1);
  EXPECT_TRUE(child0->hasGroup("empty"));

  //Destroy child3 using index argument
  IndexType idx3 = group1->getGroupIndex("child3");
  group1->destroyGroupAndData(idx3);

  //intbuf and fltbuf both lose one attached View.
  EXPECT_EQ(intbuf->getNumViews(), 4);
  EXPECT_EQ(fltbuf->getNumViews(), 3);

  //Verify that the Buffers' data can be accessed by other Views.
  int* int4_vals = child4->getView("intview")->getData();
  for(int i = 0; i < 15; ++i)
  {
    EXPECT_EQ(int4_vals[i], i);
  }

  conduit::float64* flt4_vals = foo4->getView("fooview")->getData();
  for(int i = 0; i < 12; ++i)
  {
    EXPECT_NEAR(flt4_vals[i], (conduit::float64)(-i), 1.0e-12);
  }

  //Destroy Groups held by child1. This removes "foo" and "empty" but leaves
  //the Views held by child1 in place.
  EXPECT_EQ(child1->getNumGroups(), 2);
  EXPECT_EQ(child1->getNumViews(), 3);
  child1->destroyGroupsAndData();
  EXPECT_EQ(child1->getNumGroups(), 0);
  EXPECT_EQ(child1->getNumViews(), 3);

  //That removed one more View from fltbuf but left intbuf unchanged.
  EXPECT_EQ(intbuf->getNumViews(), 4);
  EXPECT_EQ(fltbuf->getNumViews(), 2);

  //Destroy the entire subtree of child4
  EXPECT_EQ(child4->getNumGroups(), 2);
  EXPECT_EQ(child4->getNumViews(), 3);
  child4->destroyGroupSubtreeAndData();
  EXPECT_EQ(child4->getNumGroups(), 0);
  EXPECT_EQ(child4->getNumViews(), 0);

  //Both buffers lost one more View.
  EXPECT_EQ(intbuf->getNumViews(), 3);
  EXPECT_EQ(fltbuf->getNumViews(), 1);

  //THe View at "group1/child2/foo/fooview" is the only View
  //still attached to fltbuf.
  EXPECT_EQ(ds->getNumBuffers(), 2);
  EXPECT_TRUE(group1->hasView("child2/foo/fooview"));
  EXPECT_TRUE(group1->getView("child2/foo/fooview")->hasBuffer());
  EXPECT_EQ(group1->getView("child2/foo/fooview")->getBuffer(), fltbuf);

  //Destroy entire subtree of group1. This will detach the last View
  //from fltbuf and cause fltbuf to be destroyed.
  group1->destroyGroupSubtreeAndData();
  EXPECT_EQ(ds->getNumBuffers(), 1);
  EXPECT_TRUE(ds->hasBuffer(int_idx));
  EXPECT_FALSE(ds->hasBuffer(flt_idx));
  EXPECT_EQ(group1->getNumViews(), 0);
  EXPECT_EQ(group1->getNumGroups(), 0);

  //intbuf still is attached to the "intview" Views in child0 and child1.
  EXPECT_EQ(intbuf->getNumViews(), 2);
  EXPECT_EQ(group0->getView("child0/intview")->getBuffer(), intbuf);
  EXPECT_EQ(group0->getView("child1/intview")->getBuffer(), intbuf);

  //Destroy everything below root, and the remaining buffer will be destroyed.
  root->destroyGroupSubtreeAndData();
  EXPECT_FALSE(ds->hasBuffer(int_idx));
  EXPECT_EQ(ds->getNumBuffers(), 0);
  EXPECT_EQ(root->getNumViews(), 0);
  EXPECT_EQ(root->getNumGroups(), 0);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, group_name_collisions)
{
  DataStore* ds = new DataStore();
  Group* flds = ds->getRoot()->createGroup("fields");
  flds->createView("a");

  EXPECT_TRUE(flds->hasChildView("a"));

  // attempt to create duplicate group name
  EXPECT_TRUE(ds->getRoot()->createGroup("fields") == nullptr);

  // attempt to create duplicate view name.
  EXPECT_TRUE(flds->createView("a") == nullptr);

  // attempt to create a group named the same as an existing view
  EXPECT_TRUE(flds->createGroup("a") == nullptr);

  // attempt to create a view named the same as an existing group
  EXPECT_TRUE(ds->getRoot()->createView("fields") == nullptr);

  ds->getRoot()->createGroup("here//is/path");
  ds->getRoot()->createGroup("éch≈o/Ωd");
  ds->getRoot()->createGroup("../group/..");

  IndexType idx = ds->getRoot()->getFirstValidGroupIndex();
  while(indexIsValid(idx))
  {
    std::cout << ds->getRoot()->getGroup(idx)->getName() << std::endl;
    idx = ds->getRoot()->getNextValidGroupIndex(idx);
  }

  delete ds;
}
//------------------------------------------------------------------------------

TEST(sidre_group, view_copy_move)
{
  DataStore* ds = new DataStore();
  Group* flds = ds->getRoot()->createGroup("fields");
  int* buffdata;
  int extdata[10];

  View* views[6];
  std::string names[6];

  // Create view in different states
  views[0] = flds->createView("empty0");
  views[1] = flds->createView("empty1", INT_ID, 10);
  views[2] = flds->createViewAndAllocate("buffer", INT_ID, 10);
  views[3] = flds->createView("external", INT_ID, 10)->setExternalDataPtr(extdata);
  views[4] = flds->createViewScalar("scalar", 25);
  views[5] = flds->createViewString("string", "I am string");

  buffdata = flds->getView("buffer")->getData();
  for(int i = 0; i < 10; ++i)
  {
    extdata[i] = i;
    buffdata[i] = i + 100;
  }

  for(int i = 0; i < 6; ++i)
  {
    names[i] = views[i]->getName();
    EXPECT_TRUE(flds->hasView(names[i]));
  }

  // test moving a view from flds to sub1
  Group* sub1 = flds->createGroup("sub1");

  // flds->print();

  for(int i = 0; i < 6; ++i)
  {
    sub1->moveView(views[i]);
    EXPECT_FALSE(flds->hasView(names[i]));
    EXPECT_TRUE(sub1->hasView(names[i]));

    // moving to same group is a no-op
    sub1->moveView(views[i]);
    EXPECT_TRUE(sub1->hasView(names[i]));
  }

  // flds->print();

  Group* sub2 = flds->createGroup("sub2");

  for(int i = 0; i < 6; ++i)
  {
    sub2->copyView(views[i]);
    EXPECT_TRUE(sub1->hasView(names[i]));
    EXPECT_TRUE(sub2->hasView(names[i]));
  }

  // Check copies
  View* view1 = sub1->getView("empty0");
  View* view2 = sub2->getView("empty0");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isEmpty());
  EXPECT_FALSE(view2->isDescribed());
  EXPECT_FALSE(view2->isAllocated());
  EXPECT_FALSE(view2->isApplied());

  view1 = sub1->getView("empty1");
  view2 = sub2->getView("empty1");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isEmpty());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_FALSE(view2->isAllocated());
  EXPECT_FALSE(view2->isApplied());

  view1 = sub1->getView("buffer");
  view2 = sub2->getView("buffer");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->hasBuffer());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_TRUE(view2->isAllocated());
  EXPECT_TRUE(view2->isApplied());
  EXPECT_TRUE(view1->getBuffer() == view2->getBuffer());
  EXPECT_EQ(2, view1->getBuffer()->getNumViews());

  view1 = sub1->getView("external");
  view2 = sub2->getView("external");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isExternal());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_TRUE(view2->isAllocated());
  EXPECT_TRUE(view2->isApplied());
  EXPECT_EQ(view1->getVoidPtr(), view2->getVoidPtr());

  view1 = sub1->getView("scalar");
  view2 = sub2->getView("scalar");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isScalar());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_TRUE(view2->isAllocated());
  EXPECT_TRUE(view2->isApplied());
  EXPECT_EQ(view1->getData<int>(), view2->getData<int>());

  view1 = sub1->getView("string");
  view2 = sub2->getView("string");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isString());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_TRUE(view2->isAllocated());
  EXPECT_TRUE(view2->isApplied());
  const char* svalue = view1->getString();
  EXPECT_TRUE(strcmp("I am string", svalue) == 0);

  // flds->print();

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, groups_move_copy)
{
  DataStore* ds = new DataStore();
  Group* flds = ds->getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  Group* gc = flds->createGroup("c");

  const double f0value = 100.0;
  const double val = 101.0;

  Group* bschild = gb->createGroup("childOfB");

  ga->createView("i0")->setScalar(1);
  gb->createView("f0")->setScalar(f0value);
  gc->createView("d0")->setScalar(3000.0);
  bschild->createView("val")->setScalar(val);

  int buffercount = ds->getNumBuffers();

  // check that all sub groups exist
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("b"));
  EXPECT_TRUE(flds->hasGroup("c"));

  // move "b" to a child of "sub"
  EXPECT_EQ(1, gb->getIndex());
  EXPECT_EQ(flds, gb->getParent());
  Group* gsub = flds->createGroup("sub");
  Group* gb0 = gsub->moveGroup(gb);

  // gb0 is an alias to gb
  EXPECT_EQ(gb, gb0);
  EXPECT_EQ(0, gb->getIndex());
  EXPECT_EQ(gsub, gb->getParent());
  EXPECT_EQ(gb0->getNumGroups(), 1);
  EXPECT_EQ(gb0->getGroup("childOfB"), bschild);
  EXPECT_EQ(bschild->getNumGroups(), 0);
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  EXPECT_EQ(gb0->getNumViews(), 1);
  EXPECT_TRUE(gb0->hasChildView("f0"));
  if(gb0->hasChildView("f0"))
  {
    EXPECT_EQ((double)gb0->getView("f0")->getScalar(), f0value);
  }

  EXPECT_EQ(bschild->getNumViews(), 1);
  EXPECT_TRUE(bschild->hasChildView("val"));
  if(gb0->hasChildView("val"))
  {
    EXPECT_EQ((double)gb0->getView("val")->getScalar(), val);
  }

  EXPECT_EQ(flds->getNumGroups(), 3);
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("sub"));
  EXPECT_TRUE(flds->hasGroup("c"));

  EXPECT_EQ(flds->getGroup("sub")->getGroup("b"), gb);

  // verify that we can copy a group into an empty group
  Group* containCopy = ds->getRoot()->createGroup("containCopy");
  Group* theCopy = containCopy->copyGroup(flds);
  EXPECT_TRUE(theCopy->isEquivalentTo(flds));
  EXPECT_EQ(containCopy->getNumGroups(), 1);
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  // verify that we can copy a group, when there is no name clash
  Group* anotherCopy = ds->getRoot()->createGroup("anotherCopy");
  anotherCopy->createGroup("futureSiblingGroup");
  Group* theOtherCopy = anotherCopy->copyGroup(flds);
  EXPECT_EQ(anotherCopy->getNumGroups(), 2);
  EXPECT_TRUE(theOtherCopy->isEquivalentTo(flds));
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  // verify that we cannot copy a group when there is a name clash
  Group* otherB = containCopy->createGroup("b");
  otherB->createView("f1")->setScalar(42.0);
  otherB->createGroup("Q");
  Group* triedCopy = gsub->copyGroup(otherB);
  EXPECT_EQ(triedCopy, static_cast<void*>(nullptr));
  EXPECT_EQ(gsub->getNumGroups(), 1);
  EXPECT_TRUE(gsub->hasChildGroup("b"));
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  EXPECT_EQ(gb0->getNumGroups(), 1);
  EXPECT_EQ(gb0->getGroup("childOfB"), bschild);
  EXPECT_EQ(gb0->getNumViews(), 1);
  EXPECT_TRUE(gb0->hasChildView("f0"));
  if(gb0->hasChildView("f0"))
  {
    EXPECT_EQ((double)gb0->getView("f0")->getScalar(), f0value);
  }

  delete ds;
}

// Return vector of known allocator ids.
std::vector<int> getKnownAllocIds()
{
  std::vector<int> allocIds(1, axom::MALLOC_ALLOCATOR_ID);
#ifdef AXOM_USE_UMPIRE
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Host>());
  #ifdef AXOM_USE_GPU
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Device>());
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Unified>());
  // Does it make sense to check Pinned and Constant memory spaces?
  #endif
#endif
  return allocIds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, scalar_memory_allocator)
{
  std::vector<int> allocIds = getKnownAllocIds();

  DataStore ds;
  Group* grp = ds.getRoot()->createGroup("grp");

  for(auto allocId : allocIds)
  {
    grp->setDefaultTupleAllocator(allocId);

    std::int32_t scalarInt = 12345;
    grp->createViewScalar("scalarInt", scalarInt);

    std::int32_t* scalarIntPtr = static_cast<std::int32_t*>(grp->getView("scalarInt")->getVoidPtr());
    auto allocIdScalarInt = axom::getAllocatorIDFromPointer(scalarIntPtr);
    EXPECT_EQ(allocIdScalarInt, allocId);

    std::int32_t scalarIntTmp;
    axom::copy(&scalarIntTmp, scalarIntPtr, sizeof(scalarInt));
    EXPECT_EQ(scalarIntTmp, scalarInt);

    grp->destroyViewsAndData();
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, group_deep_copy)
{
  DataStore* ds = new DataStore();
  Group* flds = ds->getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");

  double dval0 = 100.0;
  double dval1 = 301.0;

  ga->createView("i0")->setScalar(1);
  ga->createView("d0")->setScalar(dval0);
  gb->createView("d1")->setScalar(dval1);
  gb->createView("s0")->setString("my string");

  // check that all sub groups exist
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("b"));

  int viewlen = 8;
  View* ownsbuf = ga->createViewAndAllocate("ownsbuf", DataType::c_int(viewlen));

  int* int_vals = ownsbuf->getData();
  for(int i = 0; i < viewlen; ++i)
  {
    int_vals[i] = i + 1;
  }

  IndexType buflen = 24;
  Buffer* dbuff = ds->createBuffer()->allocate(FLOAT64_ID, buflen);

  dbuff->allocate();
  conduit::float64* buf_ptr = dbuff->getData();

  for(int i = 0; i < buflen; ++i)
  {
    buf_ptr[i] = 2.0 * conduit::float64(i);
  }

  const int NUM_VIEWS = 4;
  int size[NUM_VIEWS] = {5, 4, 10, 11};
  int stride[NUM_VIEWS] = {3, 2, 2, 1};
  int offset[NUM_VIEWS] = {2, 9, 0, 10};
  std::string names[NUM_VIEWS] = {"viewa", "viewb", "viewc", "viewd"};

  for(int i = 0; i < NUM_VIEWS; ++i)
  {
    ga->createView(names[i], dbuff)->apply(size[i], offset[i], stride[i]);
  }

  const int extlen = 30;
  conduit::float64 ext_array[extlen];
  for(int i = 0; i < extlen; ++i)
  {
    ext_array[i] = -1.0 * conduit::float64(i);
  }

  for(int i = 0; i < NUM_VIEWS; ++i)
  {
    gb->createView(names[i], ext_array)->apply(FLOAT64_ID, size[i], offset[i], stride[i]);
  }

  Group* deep_copy = ds->getRoot()->createGroup("deep_copy");
  deep_copy->deepCopyGroup(flds);

  EXPECT_TRUE(deep_copy->hasGroup("fields/a"));
  EXPECT_TRUE(deep_copy->hasGroup("fields/b"));

  Group* copy_ga = deep_copy->getGroup("fields/a");
  Group* copy_gb = deep_copy->getGroup("fields/b");

  int io_val = copy_ga->getView("i0")->getScalar();
  EXPECT_EQ(io_val, 1);
  EXPECT_NEAR(copy_ga->getView("d0")->getScalar(), dval0, 1.0e-12);
  EXPECT_NEAR(copy_gb->getView("d1")->getScalar(), dval1, 1.0e-12);
  EXPECT_EQ(copy_gb->getView("s0")->getString(), std::string("my string"));

  for(int i = 0; i < NUM_VIEWS; ++i)
  {
    EXPECT_TRUE(copy_ga->hasView(names[i]));
    View* copy_view = copy_ga->getView(names[i]);
    EXPECT_TRUE(copy_view->hasBuffer());

    // The deep copy creates a compact buffer in the copied View, associated
    // only with that View.
    Buffer* buffer = copy_view->getBuffer();
    EXPECT_EQ(buffer->getNumViews(), 1);
    EXPECT_EQ(buffer->getNumElements(), size[i]);
    EXPECT_EQ(copy_view->getOffset(), 0);
    EXPECT_EQ(copy_view->getStride(), 1);

    conduit::float64* fdata = copy_view->getData();
    for(int j = 0; j < size[i]; ++j)
    {
      EXPECT_NEAR(fdata[j], 2.0 * (offset[i] + j * stride[i]), 1.0e-12);
    }
  }

  for(int i = 0; i < NUM_VIEWS; ++i)
  {
    EXPECT_TRUE(copy_gb->hasView(names[i]));
    View* copy_view = copy_gb->getView(names[i]);
    EXPECT_TRUE(copy_view->hasBuffer());

    // The deep copy creates a compact buffer in the copied View, associated
    // only with that View.  Here external data was copied to internal
    // buffers.
    Buffer* buffer = copy_view->getBuffer();
    EXPECT_EQ(buffer->getNumViews(), 1);
    EXPECT_EQ(buffer->getNumElements(), size[i]);
    EXPECT_EQ(copy_view->getOffset(), 0);
    EXPECT_EQ(copy_view->getStride(), 1);

    conduit::float64* fdata = copy_view->getData();
    for(int j = 0; j < size[i]; ++j)
    {
      EXPECT_NEAR(fdata[j], -1.0 * (offset[i] + j * stride[i]), 1.0e-12);
    }
  }

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, deep_copy_interspace)
{
  // Test deep copies with a change in memory space.

  std::vector<int> allocIds = getKnownAllocIds();

  DataStore ds;

  constexpr int N = 5;
  auto dtypeDouble = conduit::DataType::float64(1);
  auto dtypeIntArray = conduit::DataType::int32(N);

  // Original values on host
  const double doubleValue = 10.13456;
  axom::Array<std::int32_t> intArray(N, N);
  for(int i = 0; i < N; ++i)
  {
    intArray[i] = 1001.5 + i;
  }
  std::string stringValue = "a string";

  // Temporary space for correctness checks
  double tmpDoubleValue = doubleValue;
  axom::Array<std::int32_t> tmpIntArray(N, N);
  axom::Array<char> tmpCharArray;

  Group* srcGrandparent = ds.getRoot()->createGroup("srcGrandparent");
  Group* dstGrandparent = ds.getRoot()->createGroup("dstGrandparent");

  // For each combination of source and destination alloc ids,
  // allocate src and copy to dst.

  for(size_t si = 0; si < allocIds.size(); ++si)
  {
    int srcArrayAllocId = allocIds[si];
    int srcTupleAllocId = allocIds[allocIds.size() - 1 - si];

    Group* srcParent = srcGrandparent->createGroup("testParent");
    srcParent->setDefaultArrayAllocator(srcArrayAllocId);
    srcParent->setDefaultTupleAllocator(srcTupleAllocId);
    Group* src = srcParent->createGroup("testGrp");
    EXPECT_EQ(srcParent->getDefaultAllocatorID(), srcArrayAllocId);
    EXPECT_EQ(src->getDefaultAllocatorID(), srcArrayAllocId);

    //
    // Create, initialize and sanity-check src objects to test.
    //

    View* srcString = src->createViewString("aString", stringValue);
    const char* srcStringPtr = static_cast<char*>(srcString->getNode().data_ptr());
    EXPECT_EQ(axom::getAllocatorIDFromPointer(srcStringPtr), srcTupleAllocId);

    View* srcScalar = src->createViewScalar("aDouble", doubleValue);
    double* srcScalarPtr = static_cast<double*>(srcScalar->getNode().data_ptr());
    EXPECT_EQ(axom::getAllocatorIDFromPointer(srcScalarPtr), srcTupleAllocId);
    axom::copy(srcScalarPtr, &doubleValue, sizeof(double));

    View* srcArray = src->createViewAndAllocate("aIntArray", dtypeIntArray);
    std::int32_t* srcArrayPtr = static_cast<std::int32_t*>(srcArray->getNode().data_ptr());
    EXPECT_EQ(axom::getAllocatorIDFromPointer(srcArrayPtr), srcArrayAllocId);
    axom::copy(srcArrayPtr, intArray.data(), N * sizeof(std::int32_t));

    if(axom::execution_space<axom::SEQ_EXEC>::usesAllocId(srcArrayAllocId) &&
       axom::execution_space<axom::SEQ_EXEC>::usesAllocId(srcTupleAllocId))
    {
      std::cout << "srcGrandparent group:" << std::endl;
      srcGrandparent->print();
    }

    //
    // Copy the source into destinations of different alloc ids
    //

    for(size_t di = 0; di < allocIds.size(); ++di)
    {
      int dstArrayAllocId = allocIds[di];
      int dstTupleAllocId = allocIds[allocIds.size() - 1 - di];

      std::cout << "Testing copying allocator id " << srcArrayAllocId << " and " << srcTupleAllocId
                << " to " << dstArrayAllocId << " and " << dstTupleAllocId << std::endl;

      dstGrandparent->setDefaultArrayAllocator(dstArrayAllocId);
      dstGrandparent->setDefaultTupleAllocator(dstTupleAllocId);
      dstGrandparent->deepCopyGroup(srcParent);

      axom::sidre::Group* dstParent = dstGrandparent->getGroup(srcParent->getName());
      Group* dst = dstParent->getGroup(src->getName());

      //
      // Check pointers.  Copy data to temporary host buffers and check.
      //

      double* dstScalarPtr = static_cast<double*>(dst->getView(srcScalar->getName())->getVoidPtr());
      EXPECT_NE(dstScalarPtr, nullptr);
      EXPECT_NE(dstScalarPtr, srcScalarPtr);
      EXPECT_EQ(axom::getAllocatorIDFromPointer(dstScalarPtr), dstTupleAllocId);
      axom::copy(&tmpDoubleValue, dstScalarPtr, sizeof(double));
      EXPECT_EQ(tmpDoubleValue, doubleValue);

      std::int32_t* dstArrayPtr =
        static_cast<std::int32_t*>(dst->getView(srcArray->getName())->getVoidPtr());
      EXPECT_NE(dstArrayPtr, nullptr);
      EXPECT_NE(dstArrayPtr, srcArrayPtr);
      EXPECT_EQ(axom::getAllocatorIDFromPointer(dstArrayPtr), dstArrayAllocId);
      axom::copy(tmpIntArray.data(), srcArrayPtr, N * sizeof(std::int32_t));
      for(int i = 0; i < N; ++i)
      {
        EXPECT_EQ(tmpIntArray[i], intArray[i]);
      }

      char* dstStringPtr = static_cast<char*>(dst->getView(srcString->getName())->getVoidPtr());
      EXPECT_NE(dstStringPtr, nullptr);
      EXPECT_NE(dstStringPtr, srcStringPtr);
      EXPECT_EQ(axom::getAllocatorIDFromPointer(dstStringPtr), dstTupleAllocId);
      tmpCharArray.resize(srcString->getNumElements());
      axom::copy(tmpCharArray.data(), srcStringPtr, srcString->getNumElements());
      for(int i = 0; i < N; ++i)
      {
        EXPECT_EQ(tmpCharArray[i], stringValue[i]);
      }

      dstGrandparent->destroyGroupAndData(srcParent->getName());
    }

    srcGrandparent->destroyGroupAndData(srcParent->getName());
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, create_destroy_view_and_buffer2)
{
  DataStore* const ds = new DataStore();
  Group* const grp = ds->getRoot()->createGroup("grp");

  std::string viewName1("viewBuffer1");
  std::string viewName2("viewBuffer2");

  View* view1 = grp->createViewAndAllocate(viewName1, INT_ID, 1);
  View* view2 = grp->createViewAndAllocate(viewName2, INT_ID, 1);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ(grp->getView(viewName1), view1);

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ(grp->getView(viewName2), view2);

  IndexType const bufferId1 = view1->getBuffer()->getIndex();

  grp->destroyViewAndData(viewName1);

  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumBuffers(), 1u);

  Buffer const* const buffer1 = ds->getBuffer(bufferId1);
  EXPECT_TRUE(buffer1 == nullptr);

  View const* const view3 = grp->createView("viewBuffer3");
  grp->destroyViewsAndData();
  // should be no-op
  grp->destroyViewsAndData();
  // shut up compiler about unused variable
  (void)view3;

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, create_destroy_alloc_view_and_buffer)
{
  DataStore* const ds = new DataStore();
  Group* const grp = ds->getRoot()->createGroup("grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  // use create + alloc convenience methods
  // this one is the DataType & method
  View* const view1 = grp->createViewAndAllocate(viewName1, DataType::c_int(10));

  EXPECT_TRUE(grp->hasChildView(viewName1));
  EXPECT_EQ(grp->getView(viewName1), view1);

  int* v1_vals = view1->getData();

  for(int i = 0; i < 10; i++)
  {
    v1_vals[i] = i;
  }

  EXPECT_EQ(view1->getNumElements(), 10u);
  EXPECT_EQ(view1->getTotalBytes(), static_cast<axom::sidre::IndexType>(10 * sizeof(int)));

  grp->destroyViewAndData(viewName1);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, create_view_of_buffer_with_schema)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  // use create + alloc convenience methods
  // this one is the DataType & method
  View* base = root->createViewAndAllocate("base", DataType::c_int(10));
  int* base_vals = base->getData();
  for(int i = 0; i < 10; i++)
  {
    if(i < 5)
    {
      base_vals[i] = 10;
    }
    else
    {
      base_vals[i] = 20;
    }
  }

  Buffer* base_buff = base->getBuffer();

  // create two views into this buffer
  //
  // view for the first 5 values
  root->createView("sub_a", base_buff)->apply(DataType::c_int(5));

  int* sub_a_vals = root->getView("sub_a")->getData();

  for(int i = 0; i < 5; i++)
  {
    EXPECT_EQ(sub_a_vals[i], 10);
  }

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, copy_to_conduit_node)
{
  DataStore* ds1 = new DataStore();

  // group_a uses map format, group_b uses list format.
  Group* group_a = ds1->getRoot()->createGroup("group_a", false);
  Group* group_b = ds1->getRoot()->createGroup("group_b", true);

  // add child groups and views to group_a
  Group* foo_a = group_a->createGroup("foo_a");
  Group* bar_a = group_a->createGroup("bar_a");
  group_a->createViewScalar("i0", 1);

  foo_a->createViewScalar<int>("i1", 5);
  bar_a->createViewScalar<double>("d0", 11.0);

  // add child groups and views to group_b
  Group* foo_b = group_b->createUnnamedGroup();
  Group* bar_b = group_b->createUnnamedGroup();
  group_b->createViewScalar("", -1);

  foo_b->createViewScalar<int32_t>("i3", 15);
  bar_b->createViewScalar<double>("d1", 33.0);

  conduit::Node n;
  ds1->getRoot()->copyToConduitNode(n);

  // copyToConduitNode converts the View contents into strings in a
  // JSON layout, so we verify the values as strings

  // We can directly access group_a by name
  EXPECT_EQ(n["groups/group_a/views/i0/value"].as_string(), std::string("1"));
  EXPECT_EQ(n["groups/group_a/groups/foo_a/views/i1/value"].as_string(), std::string("5"));
  EXPECT_EQ(n["groups/group_a/groups/bar_a/views/d0/value"].as_string(), std::string("11.0"));

  // group_b has lists, so we iterate
  conduit::Node& b_views = n["groups/group_b/views"];
  conduit::Node& b_groups = n["groups/group_b/groups"];
  EXPECT_TRUE(b_views.dtype().is_list());
  EXPECT_TRUE(b_groups.dtype().is_list());
  EXPECT_EQ(b_views.number_of_children(), 1);
  EXPECT_EQ(b_groups.number_of_children(), 2);

  // Verify the View held directly in group_b
  conduit::NodeIterator v_itr = b_views.children();
  while(v_itr.has_next())
  {
    conduit::Node& chld = v_itr.next();
    EXPECT_EQ(chld["value"].as_string(), std::string("-1"));
  }

  // Verify the Views held in group_b's child Groups
  conduit::NodeIterator g_itr = b_groups.children();
  while(g_itr.has_next())
  {
    conduit::Node& chld = g_itr.next();
    EXPECT_TRUE(chld.has_child("views"));

    conduit::Node& chld_views = chld["views"];
    EXPECT_TRUE(chld_views.number_of_children() == 1);
    EXPECT_TRUE(chld_views.has_child("i3") || chld_views.has_child("d1"));

    if(chld_views.has_child("i3"))
    {
      EXPECT_EQ(chld_views["i3/value"].as_string(), std::string("15"));
    }
    else
    {
      EXPECT_EQ(chld_views["d1/value"].as_string(), std::string("33.0"));
    }
  }

  delete ds1;
}

//------------------------------------------------------------------------------
TEST(sidre_group, copy_array_to_conduit_node)
{
  DataStore ds1;

  const axom::IndexType count = 5;

  // group_a uses map format, group_b uses list format.
  Group* group_a = ds1.getRoot()->createGroup("group_a", false);
  Group* group_b = ds1.getRoot()->createGroup("group_b", true);

  // add child groups and views to group_a
  View* view_a =
    group_a->createViewAndAllocate("i0", axom::sidre::detail::SidreTT<std::int32_t>::id, count);
  View* view_b =
    group_b->createViewAndAllocate("d0", axom::sidre::detail::SidreTT<axom::float64>::id, count);

  std::int32_t* v_a = view_a->getArray();
  axom::float64* v_b = view_b->getArray();

  for(IndexType i = 0; i < count; ++i)
  {
    v_a[i] = 100 + i;
    v_b[i] = 10 + i;
  }
  EXPECT_NE(v_a, nullptr);
  EXPECT_NE(v_b, nullptr);

  conduit::Node shallowDst;
  ds1.getRoot()->createNativeLayout(shallowDst);
  EXPECT_TRUE(shallowDst.has_path("group_a/i0"));
  EXPECT_TRUE(shallowDst.has_path("group_b"));
  EXPECT_TRUE(shallowDst.fetch_existing("group_b").number_of_children() == 1);

  {
    std::int32_t* n_a = shallowDst.fetch_existing("group_a/i0").value();
    axom::float64* n_b = shallowDst.fetch_existing("group_b").child(0).value();
    EXPECT_NE(n_a, nullptr);
    EXPECT_NE(n_b, nullptr);
    EXPECT_EQ(n_a, v_a);
    EXPECT_EQ(n_b, v_b);
  }

  conduit::Node deepDst;
  ds1.getRoot()->deepCopyToConduit(deepDst);
  EXPECT_TRUE(deepDst.has_path("group_a/i0"));
  EXPECT_TRUE(deepDst.has_path("group_b"));
  EXPECT_TRUE(deepDst.fetch_existing("group_b").number_of_children() == 1);

  {
    std::int32_t* n_a = deepDst.fetch_existing("group_a/i0").value();
    axom::float64* n_b = deepDst.fetch_existing("group_b").child(0).value();
    EXPECT_NE(n_a, nullptr);
    EXPECT_NE(n_b, nullptr);
    EXPECT_NE(n_a, v_a);
    EXPECT_NE(n_b, v_b);
    // Verify correct deep-copy values
    const int hostAllocId = axom::MALLOC_ALLOCATOR_ID;
    axom::Array<std::int32_t> n_a_host {axom::ArrayView<std::int32_t> {n_a, {count}}, hostAllocId};
    axom::Array<axom::float64> n_b_host {axom::ArrayView<axom::float64> {n_b, {count}}, hostAllocId};
    axom::IndexType errCount = 0;
    for(IndexType i = 0; i < count; ++i)
    {
      errCount += n_a_host[i] != v_a[i];
      errCount += n_b_host[i] != v_b[i];
    }
    EXPECT_EQ(errCount, 0);
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_empty_datastore)
{
  const std::string file_path_base("sidre_empty_datastore_");
  DataStore* ds1 = new DataStore();

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(ds1->getRoot()->save(file_path, protocol));
  }

  delete ds1;

  // Only restore default protocol
  {
    const std::string& default_protocol = Group::getDefaultIOProtocol();
    const std::string file_path = file_path_base + default_protocol;

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();
    EXPECT_TRUE(root2->load(file_path, default_protocol));

    EXPECT_TRUE(ds2->getNumBuffers() == 0);
    EXPECT_TRUE(root2->getNumGroups() == 0);
    EXPECT_TRUE(root2->getNumViews() == 0);

    delete ds2;
  }
}

#ifdef AXOM_USE_HDF5
//------------------------------------------------------------------------------
// make sure the hdf5 methods are consistent with the path based methods
//------------------------------------------------------------------------------
TEST(sidre_group, save_load_via_hdf5_ids)
{
  DataStore ds_save;
  // populate the datastore
  Group* root = ds_save.getRoot();
  root->createViewScalar<int>("i0", 1);
  root->createViewAndAllocate("vals", INT_ID, 5);
  // set values for the "vals" array
  int* vals_ptr = root->getView("vals")->getData();
  for(int i = 0; i < 5; ++i)
  {
    vals_ptr[i] = i;
  }

  // save using the sidre_hdf5 protocol
  EXPECT_TRUE(root->save("out_save_load_via_hdf5_ids.sidre_hdf5", "sidre_hdf5"));

  // load via path based
  DataStore ds_load_generic;
  EXPECT_TRUE(
    ds_load_generic.getRoot()->load("out_save_load_via_hdf5_ids.sidre_hdf5", "sidre_hdf5"));

  // load via hdf5 id
  DataStore ds_load_hdf5;

  hid_t h5_id = H5Fopen("out_save_load_via_hdf5_ids.sidre_hdf5", H5F_ACC_RDWR, H5P_DEFAULT);
  EXPECT_TRUE(h5_id >= 0);

  // this implies protocol == "sidre_hdf5"
  EXPECT_TRUE(ds_load_hdf5.getRoot()->load(h5_id));

  // ? Does isEquivalentTo check values?
  // check path based with source
  EXPECT_TRUE(ds_load_generic.getRoot()->isEquivalentTo(ds_save.getRoot()));

  // check hdf5 based with source
  EXPECT_TRUE(ds_load_hdf5.getRoot()->isEquivalentTo(ds_save.getRoot()));

  // check path based vs hdf5 based
  EXPECT_TRUE(ds_load_generic.getRoot()->isEquivalentTo(ds_load_hdf5.getRoot()));

  // close hdf5 handle
  EXPECT_TRUE(H5Fclose(h5_id) >= 0);
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_root_restore_as_child)
{
  // We'll save the DataStore's root Group into a file then restore it
  // into a child group of another DataStore.

  // Create a DataStore; put some groups and views into it
  const std::string file_path_base("sidre_save_root_restore_as_child_");
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* child1 = root->createGroup("g_a");
  Group* child2 = root->createGroup("g_b");
  root->createGroup("g_c");  // We don't put anything into g_c

  child1->createViewScalar<int>("i0", 1);
  child1->createViewString("s0", "I am a string");

  const int num_elems = 4;
  int* pa0 = child2->createViewAndAllocate("a0", INT_ID, num_elems)->getArray();
  double* pa1 = child2->createViewAndAllocate("a1", FLOAT64_ID, num_elems)->getArray();

  const double factor = 2.3;
  const double offset = -0.23;
  for(int i = 0; i < num_elems; ++i)
  {
    pa0[i] = i;
    pa1[i] = offset + i * factor;
  }

  // Save the DataStore's root
  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(root->save(file_path, protocol));
  }

  // Restore the original DataStore into a child group
  for(const auto& protocol : protocols)
  {
    // Only restore sidre_hdf5 protocol
    if(protocol != "sidre_hdf5")
    {
      continue;
    }

    DataStore* dscopy = new DataStore();
    Group* dsroot = dscopy->getRoot();
    const std::string file_path = file_path_base + protocol;

    const std::string group_base("group_");
    const std::string group_name = group_base + protocol;
    Group* cg = dsroot->createGroup(group_name);

    if(axom::utilities::filesystem::pathExists(file_path))
    {
      std::cout << "loading " << file_path << std::endl;
      EXPECT_TRUE(cg->load(file_path, protocol));

      EXPECT_TRUE(cg->isEquivalentTo(root, false));
      EXPECT_TRUE(root->isEquivalentTo(cg, false));
    }
    else
    {
      FAIL() << "file not present: " << file_path;
    }

    delete dscopy;
  }

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_child_restore_as_root)
{
  // We'll save a child Group into a file then restore it as the root of
  // a new DataStore.

  // Create a DataStore; put some groups and views into it
  const std::string file_path_base("sidre_save_child_restore_as_root_");
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  // child1 and all its descendents will get saved into files
  Group* child1 = root->createGroup("g_a");
  Group* child2 = child1->createGroup("g_b");
  child1->createViewScalar<int>("i0", 1);
  child1->createViewString("s0", "I am a string");
  // everything else won't be put into the files
  Group* child1a = root->createGroup("g_notSaved");
  child1a->createViewScalar<int>("i1", 42);
  root->createViewString("s1", "string view off the root, not saved");

  const int num_elems = 4;
  // included in files
  int* pa0 = child2->createViewAndAllocate("a0", INT_ID, num_elems)->getArray();
  double* pa1 = child2->createViewAndAllocate("a1", FLOAT64_ID, num_elems)->getArray();
  // not included
  int* pa2 = child1a->createViewAndAllocate("a2", INT_ID, num_elems)->getArray();
  int* pa3 = root->createViewAndAllocate("a3", INT_ID, num_elems)->getArray();

  const double factor = 2.3;
  const double offset = -0.23;
  for(int i = 0; i < num_elems; ++i)
  {
    pa0[i] = i;
    pa1[i] = offset + i * factor;
    pa2[i] = i + 2;
    pa3[i] = 4 - i;
  }

  // Save the Group in question (child1) into an archive
  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(child1->save(file_path, protocol));
  }

  // Restore the saved child1 into a root group
  for(const auto& protocol : protocols)
  {
    // Only restore sidre_hdf5 protocol
    if(protocol != "sidre_hdf5")
    {
      continue;
    }

    DataStore* dscopy = new DataStore();
    const std::string file_path = file_path_base + protocol;
    if(axom::utilities::filesystem::pathExists(file_path))
    {
      EXPECT_TRUE(dscopy->getRoot()->load(file_path, protocol));

      EXPECT_TRUE(dscopy->getRoot()->isEquivalentTo(child1, false));
      EXPECT_TRUE(child1->isEquivalentTo(dscopy->getRoot(), false));
    }
    else
    {
      FAIL() << "file not present: " << file_path;
    }

    delete dscopy;
  }

  delete ds;
}
#endif  // AXOM_USE_HDF5

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_api)
{
  const std::string file_path_base("sidre_save_subtree_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createViewScalar<int>("i0", 1);

  // These should be produce identical files.

  // No group provided, defaults to root group
  EXPECT_TRUE(root1->save("sidre_save_fulltree_conduit", "json"));

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(root1->save(file_path, protocol));
  }

  //These are commented out because createViewScalar<int> creates a scalar View
  //with a 32-bit int, but when conduit reads a raw integer from json file, it
  //stores it as a 64-bit int, so the isEquivalentTo test fails.

#if 0
  DataStore* ds2 = new DataStore();
  Group* root2 = ds2->getRoot();
  EXPECT_TRUE(root2->load("sidre_save_fulltree_conduit", "json"));
  EXPECT_TRUE( ds2->getRoot()->isEquivalentTo(root1) );
  delete ds2;

  DataStore* ds3 = new DataStore();
  Group* root3 = ds3->getRoot();
  EXPECT_TRUE(root3->load("sidre_save_subtree_sidre_json", "sidre_json"));
  EXPECT_TRUE( ds3->getRoot()->isEquivalentTo(root1) );
  delete ds3;
#endif

#ifdef AXOM_USE_HDF5
  DataStore* ds4 = new DataStore();
  Group* root4 = ds4->getRoot();
  EXPECT_TRUE(root4->load("sidre_save_subtree_sidre_hdf5", "sidre_hdf5"));
  EXPECT_TRUE(ds4->getRoot()->isEquivalentTo(root1));
  delete ds4;
#endif

#if 0
  DataStore* ds5 = new DataStore();
  Group* root5 = ds5->getRoot();
  EXPECT_TRUE(root5->load("sidre_save_subtree_json", "json"));
  EXPECT_TRUE( ds5->getRoot()->isEquivalentTo(root1) );
  delete ds5;
#endif
  delete ds1;

  // Test loading of same subtree to different parts of a Group tree
  DataStore* ds_new = new DataStore();
  Group* tree1 = ds_new->getRoot()->createGroup("api1");
  Group* tree2 = ds_new->getRoot()->createGroup("api2");

  Group* load1 = tree1->createGroup("subtree");
  Group* load2 = tree2->createGroup("subtree");

  EXPECT_TRUE(load1->load("sidre_save_subtree_sidre_json", "sidre_json"));
  EXPECT_TRUE(load2->load("sidre_save_subtree_sidre_json", "sidre_json"));

  EXPECT_TRUE(load1->isEquivalentTo(load2));

  std::string newgroupname = "in case of blank";
  std::string groupname = newgroupname;
  bool loadSuccess = false;
  Group* load3 =
    load2->createGroupAndLoad(groupname, "sidre_save_subtree_sidre_json", "sidre_json", loadSuccess);

  EXPECT_NE(load3, (Group*)nullptr);
  EXPECT_TRUE(loadSuccess);
  EXPECT_EQ(groupname, "");
  if(load3 != nullptr)
  {
    EXPECT_EQ(newgroupname, load3->getName());
    EXPECT_TRUE(load1->isEquivalentTo(load3, false));
  }

  std::string anothergroupname = "another group";
  groupname = anothergroupname;
  loadSuccess = false;
  Group* load4 =
    load2->createGroupAndLoad(groupname, "sidre_save_subtree_sidre_json", "sidre_json", loadSuccess);

  EXPECT_NE(load4, (Group*)nullptr);
  EXPECT_TRUE(loadSuccess);
  if(load4 != nullptr)
  {
    EXPECT_EQ(anothergroupname, load4->getName());
    EXPECT_TRUE(load3->isEquivalentTo(load4, false));
  }

  groupname = anothergroupname;
  loadSuccess = false;
  Group* load5 =
    load2->createGroupAndLoad(groupname, "sidre_save_subtree_sidre_json", "sidre_json", loadSuccess);
  EXPECT_EQ(load5, (Group*)nullptr);
  EXPECT_FALSE(loadSuccess);

  delete ds_new;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_scalars_and_strings)
{
  const std::string file_path_base("sidre_save_scalars_and_strings_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createViewScalar<int>("i0", 1);
  root1->createViewScalar<float>("f0", 1.0);
  root1->createViewScalar<double>("d0", 10.0);
  root1->createViewString("s0", "I am a string");

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    //      if ( protocol == "conduit_hdf5")
    //    continue;   // XXX - Does not work
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(root1->save(file_path, protocol));
  }

  for(const auto& protocol : protocols)
  {
    // Only restore sidre_hdf5 protocol
    if(protocol != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocol;

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    EXPECT_TRUE(root2->load(file_path, protocol));

    EXPECT_TRUE(root1->isEquivalentTo(root2));

    int i0 = root2->getView("i0")->getScalar();
    float f0 = root2->getView("f0")->getScalar();
    double d0 = root2->getView("d0")->getScalar();
    const char* s0 = root2->getView("s0")->getString();

    EXPECT_EQ(1, i0);
    EXPECT_EQ(1.0, f0);
    EXPECT_EQ(10.0, d0);
    EXPECT_EQ(std::string(s0), "I am a string");

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
TEST(sidre_group, rename_group)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* child1 = root->createGroup("g_a");
  Group* child2 = root->createGroup("g_b");
  Group* child3 = root->createGroup("g_c");

  // rename should not change the index
  EXPECT_EQ(0, child1->getIndex());
  bool success = child1->rename("g_r");
  EXPECT_TRUE(success);
  EXPECT_EQ("g_r", child1->getName());
  EXPECT_EQ(0, child1->getIndex());
  EXPECT_TRUE(root->hasGroup("g_r"));
  EXPECT_FALSE(root->hasGroup("g_a"));

  // try to rename to path
  success = child2->rename("fields/g_s");
  EXPECT_FALSE(success);
  EXPECT_EQ("g_b", child2->getName());

  // Try to rename to existing group name
  success = child3->rename("g_b");
  EXPECT_FALSE(success);
  EXPECT_EQ("g_c", child3->getName());

  // Rename root group
  EXPECT_FALSE(indexIsValid(root->getIndex()));
  EXPECT_EQ(root, root->getParent());
  EXPECT_EQ("", root->getName());
  root->rename("newroot");
  EXPECT_FALSE(indexIsValid(root->getIndex()));
  EXPECT_EQ(root, root->getParent());
  EXPECT_EQ("newroot", root->getName());

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_name_change)
{
  const std::string file_path_base("sidre_save_name_change_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();
  Group* child1 = root1->createGroup("child1");

  child1->createViewScalar<int>("i0", 1);
  child1->createViewString("s0", "I am a string");

  bool success = child1->getView("s0")->rename("s0_renamed");

  EXPECT_TRUE(success);
  EXPECT_FALSE(child1->hasView("s0"));
  EXPECT_TRUE(child1->hasView("s0_renamed"));

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    //      if ( protocol == "conduit_hdf5")
    //    continue;   // XXX - Does not work
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(child1->save(file_path, protocol));
  }

  std::string groupname;
  for(const auto& protocol : protocols)
  {
    // Only restore sidre_hdf5 protocol
    if(protocol != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocol;

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();
    Group* child2 = root2->createGroup("child2");

    EXPECT_EQ(child2->getName(), "child2");

    EXPECT_TRUE(child2->load(file_path, protocol, false, groupname));

    EXPECT_EQ(child2->getName(), "child2");

    child2->rename(groupname);

    EXPECT_TRUE(root1->isEquivalentTo(root2));

    int i0 = child2->getView("i0")->getScalar();
    const char* s0 = child2->getView("s0_renamed")->getString();

    EXPECT_EQ(1, i0);
    EXPECT_EQ(std::string(s0), "I am a string");

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_external_data)
{
  const std::string file_path_base("sidre_save_external_");

  const int nfoo = 10;
  int foo1[nfoo], foo2[nfoo], *foo3, foo4[nfoo];
  int int2d1[nfoo * 2], int2d2[nfoo * 2];
  IndexType shape[] = {nfoo, 2};

  for(int i = 0; i < nfoo; ++i)
  {
    foo1[i] = i;
    foo2[i] = 0;
    foo4[i] = i;
  }
  for(int i = 0; i < 2 * nfoo; ++i)
  {
    int2d1[i] = i;
    int2d2[i] = 0;
  }
  foo3 = nullptr;

  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createView("external_array", INT_ID, nfoo, foo1);
  root1->createView("empty_array", INT_ID, nfoo, foo3);
  // XXX this falls into createView(name, type, ndims, shape)
  // root1->createView("empty_array", INT_ID, nfoo, NULL);
  root1->createView("external_undescribed")->setExternalDataPtr(foo4);
  root1->createViewWithShape("int2d", INT_ID, 2, shape, int2d1);

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(root1->save(file_path, protocol));
  }

  delete ds1;

  // Now load back in.
  for(const auto& protocol : protocols)
  {
    // Only restore sidre_hdf5 protocol
    if(protocol != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocol;
    IndexType extents[7];
    int rank;

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    EXPECT_TRUE(root2->load(file_path, protocol));

    // load has set the type and size of the view.
    // Now set the external address before calling loadExternal.
    View* view1 = root2->getView("external_array");
    EXPECT_TRUE(view1->isExternal());
    EXPECT_TRUE(view1->isDescribed());
    EXPECT_EQ(view1->getTypeID(), INT_ID);
    EXPECT_EQ(view1->getNumElements(), nfoo);
    EXPECT_TRUE(view1->getVoidPtr() == nullptr);
    view1->setExternalDataPtr(foo2);

    View* view2 = root2->getView("empty_array");
    EXPECT_TRUE(view2->isEmpty());
    EXPECT_TRUE(view2->isDescribed());
    EXPECT_EQ(view2->getTypeID(), INT_ID);
    EXPECT_TRUE(view2->getVoidPtr() == nullptr);
    view2->setExternalDataPtr(foo3);

    View* view3 = root2->getView("external_undescribed");
    EXPECT_TRUE(view3->isEmpty());
    EXPECT_FALSE(view3->isDescribed());
    EXPECT_TRUE(view3->getVoidPtr() == nullptr);
    // Set "external_array" and "external_undescribed" to the same external
    // array
    // since it was created that way.  However, "external_undescribed" was not
    // written to the dump since it is undescribed.
    view3->setExternalDataPtr(foo2);

    View* view4 = root2->getView("int2d");
    EXPECT_FALSE(view4->isEmpty());
    EXPECT_TRUE(view4->isDescribed());
    EXPECT_TRUE(view4->getVoidPtr() == nullptr);
    EXPECT_EQ(view4->getTypeID(), INT_ID);
    EXPECT_EQ(view4->getNumElements(), nfoo * 2);
    EXPECT_EQ(view4->getNumDimensions(), 2);
    rank = view4->getShape(7, extents);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(extents[0] == nfoo && extents[1] == 2);
    view4->setExternalDataPtr(int2d2);

    // Read external data into views
    EXPECT_TRUE(root2->loadExternalData(file_path));

    // Make sure addresses have not changed
    EXPECT_TRUE(view1->getVoidPtr() == static_cast<void*>(foo2));
    EXPECT_TRUE(view2->getVoidPtr() == static_cast<void*>(foo3));  // nullptr
    EXPECT_TRUE(view3->getVoidPtr() == static_cast<void*>(foo2));
    EXPECT_TRUE(view4->getVoidPtr() == static_cast<void*>(int2d2));

    for(int j = 0; j < nfoo; ++j)
    {
      EXPECT_TRUE(foo1[j] == foo2[j]);
    }
    for(int j = 0; j < 2 * nfoo; ++j)
    {
      EXPECT_TRUE(int2d1[j] == int2d2[j]);
    }

    delete ds2;

    // Now try to load external from an incorrect path.  This will fail,
    // causing loadExternalData() to return false.
    DataStore* ds3 = new DataStore();
    Group* root3 = ds3->getRoot();
    root3->load(file_path, protocol);
    View* view31 = root3->getView("external_array");
    view31->setExternalDataPtr(foo2);
    View* view32 = root3->getView("empty_array");
    view32->setExternalDataPtr(foo3);
    View* view33 = root3->getView("external_undescribed");
    view33->setExternalDataPtr(foo2);
    View* view34 = root3->getView("int2d");
    view34->setExternalDataPtr(int2d2);

    const std::string bad_file_path = "garbage_" + file_path;
    EXPECT_FALSE(root3->loadExternalData(bad_file_path));

    delete ds3;
  }
}

//------------------------------------------------------------------------------

// Check the association between views and buffers to make sure it is what we
// expect.
// This checks more than isEquivalentTo.

static void save_restore_buffer_association(const std::string& msg, DataStore* ds)
{
  const IndexType len = 10;

  SCOPED_TRACE(msg);

  // Make sure all buffers were created
  ASSERT_EQ(ds->getNumBuffers(), 4);

  Group* root = ds->getRoot();

  // Get all views and their buffers
  View* view1 = root->getView("undescribed_attached_buffer");
  ASSERT_TRUE(view1->hasBuffer());
  Buffer* buff1a = view1->getBuffer();
  ASSERT_FALSE(buff1a->isDescribed());
  ASSERT_FALSE(buff1a->isAllocated());

  View* view2 = root->getView("unallocated_attached_buffer");
  ASSERT_TRUE(view2->hasBuffer());
  Buffer* buff2a = view2->getBuffer();
  ASSERT_TRUE(buff2a->isDescribed());
  ASSERT_FALSE(buff2a->isAllocated());

  View* view3 = root->getView("undescribed_view_described_buffer");
  ASSERT_TRUE(view3->hasBuffer());
  Buffer* buff3a = view3->getBuffer();
  ASSERT_TRUE(buff3a->isDescribed());
  ASSERT_TRUE(buff3a->isAllocated());

  View* view4 = root->getView("describe_view_described_buffer");
  ASSERT_TRUE(view4->hasBuffer());
  Buffer* buff3b = view4->getBuffer();

  View* view5 = root->getView("even");
  ASSERT_TRUE(view5->hasBuffer());
  Buffer* buff3c = view5->getBuffer();

  View* view6 = root->getView("odd");
  ASSERT_TRUE(view6->hasBuffer());
  Buffer* buff3d = view6->getBuffer();

  View* view7 = root->getView("empty_described");
  ASSERT_FALSE(view7->hasBuffer());
  ASSERT_TRUE(view7->isEmpty());

  View* view8 = root->getView("allocated");
  ASSERT_TRUE(view8->hasBuffer());
  Buffer* buff4a = view8->getBuffer();
  ASSERT_TRUE(buff4a->isDescribed());
  ASSERT_TRUE(buff4a->isAllocated());

  // Check shared buffer  3a == 3b == 3c == 3d
  ASSERT_EQ(buff3a, buff3b);
  ASSERT_EQ(buff3b, buff3c);
  ASSERT_EQ(buff3c, buff3d);

  // Check unique buffers
  ASSERT_NE(buff1a, buff2a);
  ASSERT_NE(buff1a, buff3a);
  ASSERT_NE(buff1a, buff4a);

  ASSERT_NE(buff2a, buff3a);
  ASSERT_NE(buff2a, buff4a);

  ASSERT_NE(buff3a, buff4a);

  // Check contents of buffers
  ASSERT_EQ(buff3a->getNumElements(), len);
  int* idata = buff3a->getData();
  for(int ii = 0; ii < len; ++ii)
  {
    ASSERT_EQ(idata[ii], ii + 100);
  }

  ASSERT_EQ(buff4a->getNumElements(), len);
  idata = buff4a->getData();
  for(int ii = 0; ii < len; ++ii)
  {
    ASSERT_EQ(idata[ii], ii + 200);
  }
}

TEST(sidre_group, save_restore_buffer)
{
  const std::string file_path_base("sidre_save_buffer_");
  const IndexType len = 10;

  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();
  Buffer* buff1 = ds1->createBuffer();
  Buffer* buff2 = ds1->createBuffer(INT_ID, len);
  Buffer* buff3 = ds1->createBuffer(INT_ID, len)->allocate();

  int* idata = buff3->getData();
  for(IndexType ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii + 100;
  }

  root1->createView("undescribed_attached_buffer", buff1);
  root1->createView("unallocated_attached_buffer", buff2);

  // These views share a buffer
  root1->createView("undescribed_view_described_buffer", buff3);
  root1->createView("describe_view_described_buffer", INT_ID, len, buff3);
  root1->createView("even", buff3)->apply(INT_ID, 5, 0, 2);
  root1->createView("odd", buff3)->apply(INT_ID, 5, 1, 2);

  root1->createView("empty_described", INT_ID, len);
  View* view = root1->createViewAndAllocate("allocated", INT_ID, len);

  idata = view->getData();
  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii + 200;
  }

  save_restore_buffer_association("original datastore", ds1);

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(root1->save(file_path, protocol));
  }

  // Now load back in.
  for(const auto& protocol : protocols)
  {
    // Only restore sidre_hdf5 protocol
    if(protocol != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocol;

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    EXPECT_TRUE(root2->load(file_path, protocol));

    bool isequivalent = root1->isEquivalentTo(root2);
    EXPECT_TRUE(isequivalent);
    if(isequivalent)
    {
      save_restore_buffer_association("loaded datastore", ds2);
    }

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------

TEST(sidre_group, save_restore_other)
{
  const std::string file_path_base("sidre_save_other_");
  const int ndata = 10;
  IndexType shape1[] = {ndata, 2};
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createView("empty_view");
  root1->createView("empty_described", INT_ID, ndata);
  root1->createViewWithShape("empty_shape", INT_ID, 2, shape1);

  auto* view = root1->createViewWithShapeAndAllocate("buffer_shape", INT_ID, 2, shape1);

  // fill the data
  {
    int* data = view->getData();
    for(int i = 0; i < ndata * 2; ++i)
    {
      data[i] = i;
    }
  }

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(root1->save(file_path, protocol));
  }

  delete ds1;

  // Now load back in.
  for(const auto& protocol : protocols)
  {
    // Only restore sidre_hdf5 protocol
    if(protocol != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocol;
    IndexType shape2[7];
    int rank;

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    EXPECT_TRUE(root2->load(file_path, protocol));

    View* view1 = root2->getView("empty_view");
    EXPECT_TRUE(view1->isEmpty());
    EXPECT_FALSE(view1->isDescribed());

    View* view2 = root2->getView("empty_described");
    EXPECT_TRUE(view2->isEmpty());
    EXPECT_TRUE(view2->isDescribed());
    EXPECT_EQ(view2->getTypeID(), INT_ID);
    EXPECT_EQ(view2->getNumElements(), ndata);

    View* view3 = root2->getView("empty_shape");
    EXPECT_TRUE(view3->isEmpty());
    EXPECT_TRUE(view3->isDescribed());
    EXPECT_EQ(view3->getTypeID(), INT_ID);
    EXPECT_EQ(view3->getNumElements(), ndata * 2);
    shape2[0] = 0;
    shape2[1] = 0;
    rank = view3->getShape(7, shape2);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(shape2[0] == ndata && shape2[1] == 2);

    View* view4 = root2->getView("buffer_shape");
    EXPECT_TRUE(view4->hasBuffer());
    EXPECT_TRUE(view4->isDescribed());
    EXPECT_EQ(view4->getTypeID(), INT_ID);
    EXPECT_EQ(view4->getNumElements(), ndata * 2);
    shape2[0] = 0;
    shape2[1] = 0;
    rank = view4->getShape(7, shape2);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(shape2[0] == ndata && shape2[1] == 2);

    delete ds2;
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_complex)
{
  const std::string file_path_base("sidre_mixed_types_");
  DataStore* ds1 = new DataStore();
  Group* flds = ds1->getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  Group* gc = flds->createGroup("c");
  int ndata = 10;

  ga->createViewScalar<int>("i0", 100.0);
  ga->createViewScalar<double>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", INT_ID, ndata);
  int* data_ptr = gc->getView("int10")->getArray();
  for(int i = 0; i < ndata; ++i)
  {
    data_ptr[i] = i;
  }

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    const std::string file_path = file_path_base + protocol;
    EXPECT_TRUE(ds1->getRoot()->save(file_path, protocol));
  }

  for(const auto& protocol : protocols)
  {
    // Only restore sidre_hdf5 protocol
    if(protocol != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocol;

    DataStore* ds2 = new DataStore();

    EXPECT_TRUE(ds2->getRoot()->load(file_path, protocol));

    EXPECT_TRUE(ds1->getRoot()->isEquivalentTo(ds2->getRoot()));

    flds = ds2->getRoot()->getGroup("fields");

    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getData<int>(), 100.0);
    EXPECT_NEAR(flds->getGroup("a")->getView("d0")->getData<double>(), 3000.0, 1e-12);

    int* new_data_ptr = flds->getGroup("c")->getView("int10")->getArray();
    for(int i = 0; i < ndata; ++i)
    {
      EXPECT_TRUE(new_data_ptr[i] == i);
    }

    const char* char_ptr = flds->getView("b/s0")->getString();
    EXPECT_TRUE(std::string(char_ptr) == "foo");

    //ds2->print();

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
// isEquivalentTo()
//------------------------------------------------------------------------------
TEST(sidre_group, is_equivalent_to)
{
  DataStore* ds = new DataStore();

  //These are the parents for two separate subtrees of the root group.
  //Everything below them will be created identically.
  Group* parent1 = ds->getRoot()->createGroup("parent1");
  Group* parent2 = ds->getRoot()->createGroup("parent2");

  //The flds1 and flds2 groups will be compared for equivalence
  Group* flds1 = parent1->createGroup("fields");
  Group* flds2 = parent2->createGroup("fields");

  Group* ga1 = flds1->createGroup("a");
  Group* gb1 = flds1->createGroup("b");
  Group* gc1 = flds1->createGroup("c");

  Group* gc2 = flds2->createGroup("c");  // Note: flds2 groups added in
                                         // different order
  Group* gb2 = flds2->createGroup("b");
  Group* ga2 = flds2->createGroup("a");

  ga1->createViewScalar("i0", 1);
  gb1->createViewScalar("f0", 100.0f);
  gc1->createViewScalar("d0", 3000.00);
  gc1->createViewScalar("d1", 6000.00);
  gc1->createViewScalar("d2", 9000.00);

  ga2->createViewScalar("i0", 1);
  gb2->createViewScalar("f0", 100.0f);
  gc2->createViewScalar("d2", 9000.00);  // Note: views of gc2 added in
                                         // different order
  gc2->createViewScalar("d1", 6000.00);
  gc2->createViewScalar("d0", 3000.00);

  // Groups were created identically, so should be equivalent.
  EXPECT_TRUE(flds1->isEquivalentTo(flds2));

  // Add something extra to flds2, making them not equivalent.
  gc2->createViewAndAllocate("extra", DataType::c_double());

  EXPECT_FALSE(flds1->isEquivalentTo(flds2));

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_load_all_protocols)
{
  // Note: This test relies on re-wiring conduit error handlers
  DataStore::setConduitSLICMessageHandlers();

  const std::string file_path_base("sidre_save_load_all_protocols.");
  DataStore ds;

  Group* flds = ds.getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  Group* gc = flds->createGroup("c");
  int ndata = 10;

  // prep a tree that can exactly restored by all
  // i/o protocols.
  // Specially, use int64 and float64 b/c the
  // json i/o case uses those types for parsed integers
  // and floating point numbers.

  ga->createViewScalar<conduit::int64>("i0", 100);
  ga->createViewScalar<conduit::float64>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", DataType::int64(ndata));
  conduit::int64* data_ptr = gc->getView("int10")->getArray();
  for(int i = 0; i < ndata; ++i)
  {
    data_ptr[i] = (conduit::int64)i;
  }

  // show the source tree
  SLIC_INFO("Source tree");
  ds.print();

  //
  // test all protocols
  //
  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    SLIC_INFO("Testing protocol: " << protocol);
    const std::string file_path = file_path_base + protocol;
    // save using current protocol
    EXPECT_TRUE(ds.getRoot()->save(file_path, protocol));

    DataStore ds_load;
    EXPECT_TRUE(ds_load.getRoot()->load(file_path, protocol));
    EXPECT_FALSE(ds.getConduitErrorOccurred());
    ds.appendToConduitErrors("test: dummy error");
    EXPECT_TRUE(ds.getConduitErrorOccurred());
    ds.clearConduitErrors();
    EXPECT_FALSE(ds.getConduitErrorOccurred());

    SLIC_INFO("Tree from protocol: " << protocol);
    // show the result
    ds_load.print();

    Group* ds_load_root = ds_load.getRoot();
    // check that the sidre hierarchy is equiv
    EXPECT_TRUE(ds.getRoot()->isEquivalentTo(ds_load_root));

    // check that the values are the same
    EXPECT_EQ(ds_load_root->getView("fields/a/i0")->getData<conduit::int64>(), 100);
    EXPECT_NEAR(ds_load_root->getView("fields/a/d0")->getData<conduit::float64>(), 3000.00, 1e-12);
    EXPECT_EQ(ds_load_root->getView("fields/b/s0")->getString(), std::string("foo"));

    conduit::int64* load_data_ptr = ds_load_root->getView("fields/c/int10")->getData();
    for(int j = 0; j < ndata; j++)
    {
      EXPECT_EQ(data_ptr[j], load_data_ptr[j]);
    }
  }

  // restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
TEST(sidre_group, fail_save_all_protocols)
{
  const std::string file_path_base("sidre_fail_save_all_protocols.");
  DataStore ds;

  Group* flds = ds.getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  Group* gc = flds->createGroup("c");
  int ndata = 10;

  // prep a tree that can exactly restored by all
  // i/o protocols.
  // Specially, use int64 and float64 b/c the
  // json i/o case uses those types for parsed integers
  // and floating point numbers.

  ga->createViewScalar<conduit::int64>("i0", 100);
  ga->createViewScalar<conduit::float64>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", DataType::int64(ndata));
  conduit::int64* data_ptr = gc->getView("int10")->getArray();
  for(int i = 0; i < ndata; ++i)
  {
    data_ptr[i] = (conduit::int64)i;
  }

  // show the source tree
  SLIC_INFO("Source tree");
  ds.print();

  //
  // test all protocols
  //
  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  for(const auto& protocol : protocols)
  {
    SLIC_INFO("Testing fail to save or load protocol: " << protocol);
    const std::string save_file_path = file_path_base + "save." + protocol;
    // create a directory named file_path
    axom::utilities::filesystem::makeDirsForPath(save_file_path);
    // fail to save, since there's a directory with the same name
    EXPECT_FALSE(ds.getRoot()->save(save_file_path, protocol));
    EXPECT_TRUE(ds.getConduitErrorOccurred());
    std::string errstring = ds.getConduitErrors();
    EXPECT_GT(errstring.length(), 0);

    DataStore ds_load;
    const std::string load_file_path = file_path_base + "load." + protocol;
    // fail to load, since there's a directory, not a data file
    EXPECT_FALSE(ds_load.getRoot()->load(load_file_path, protocol));
    EXPECT_TRUE(ds.getConduitErrorOccurred());
    errstring = ds.getConduitErrors();
    EXPECT_GT(errstring.length(), 0);
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_load_preserve_contents)
{
  // Note: This test relies on re-wiring conduit error handlers
  DataStore::setConduitSLICMessageHandlers();

  const std::string file_path_tree0("sidre_save_load_preserve_contents.tree0.");
  const std::string file_path_tree1("sidre_save_load_preserve_contents.tree1.");
  DataStore ds;

  Group* tree0 = ds.getRoot()->createGroup("tree0");

  Group* ga = tree0->createGroup("a");
  Group* gb = tree0->createGroup("b");
  Group* gc = tree0->createGroup("c");
  int ndata = 10;

  // Prepare a tree that can be exactly restored by all I/O protocols.
  // Specifically, use int64 and float64 because the json I/O case
  // uses those types for parsed integers and floating point numbers.

  ga->createViewScalar<conduit::int64>("i0", 100);
  ga->createViewScalar<conduit::float64>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", DataType::int64(ndata));
  conduit::int64* data_ptr = gc->getView("int10")->getArray();
  for(int i = 0; i < ndata; ++i)
  {
    data_ptr[i] = (conduit::int64)i;
  }

  const std::vector<std::string>& protocols = Group::getValidIOProtocols();
  std::string groupname;
  for(const auto& protocol : protocols)
  {
    std::string file_path0 = file_path_tree0 + protocol;
    EXPECT_TRUE(tree0->save(file_path0, protocol));

    Group* tree1 = tree0->createGroup("tree1");

    Group* gx = tree1->createGroup("x");
    Group* gy = tree1->createGroup("y");
    Group* gz = tree1->createGroup("z");

    gx->createViewAndAllocate("int20", DataType::int64(ndata * 2));
    conduit::int64* data_ptr20 = gx->getView("int20")->getArray();
    for(int i = 0; i < ndata * 2; ++i)
    {
      data_ptr20[i] = (conduit::int64)(i * 2);
    }
    gy->createViewScalar<conduit::int64>("i0", 400);
    gz->createViewScalar<conduit::float64>("d0", 17.00);

    std::string file_path1 = file_path_tree1 + protocol;

    EXPECT_TRUE(tree1->save(file_path1, protocol));

    // show the source tree
    SLIC_INFO("Source tree");
    ds.print();

    DataStore ds_load;
    Group* loadtree0 = ds_load.getRoot()->createGroup("tree0");
    EXPECT_TRUE(loadtree0->load(file_path0, protocol, false, groupname));
    EXPECT_TRUE(loadtree0->load(file_path1, protocol, true, groupname));
    loadtree0->rename(groupname);

    SLIC_INFO("Tree from protocol: " << protocol);
    // show the result
    ds_load.print();

    Group* ds_load_root = ds_load.getRoot();

    // check that the values are the same
    EXPECT_EQ(ds_load_root->getView("tree1/a/i0")->getData<conduit::int64>(), 100);
    EXPECT_NEAR(ds_load_root->getView("tree1/a/d0")->getData<conduit::float64>(), 3000.00, 1e-12);
    EXPECT_EQ(ds_load_root->getView("tree1/b/s0")->getString(), std::string("foo"));
    EXPECT_EQ(ds_load_root->getView("tree1/y/i0")->getData<conduit::int64>(), 400);
    EXPECT_NEAR(ds_load_root->getView("tree1/z/d0")->getData<conduit::float64>(), 17.00, 1e-12);

    conduit::int64* load_data_ptr = ds_load_root->getView("tree1/c/int10")->getData();
    for(int j = 0; j < ndata; j++)
    {
      EXPECT_EQ(data_ptr[j], load_data_ptr[j]);
    }
    load_data_ptr = ds_load_root->getView("tree1/x/int20")->getData();
    for(int j = 0; j < ndata * 2; j++)
    {
      EXPECT_EQ(data_ptr20[j], load_data_ptr[j]);
    }

    // Destroy the group so the name can be reused by the next protocol
    tree0->destroyGroup("tree1");
  }

  // restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_layout_protocols)
{
  // Tests calls to save using sidre_layout_json and conduit_layout_json
  // protocols. No loads are attempted as these protocols are used for
  // save only.

  DataStore::setConduitSLICMessageHandlers();

  const std::string file_path_base("sidre_save_layout_protocols_");
  DataStore ds;

  Group* flds = ds.getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  Group* gc = flds->createGroup("c");
  int ndata = 10;

  ga->createViewScalar<conduit::int64>("i0", 100);
  ga->createViewScalar<conduit::float64>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", DataType::int64(ndata));
  conduit::int64* data_ptr = gc->getView("int10")->getArray();
  for(int i = 0; i < ndata; ++i)
  {
    data_ptr[i] = (conduit::int64)i;
  }

  int extdatasize = 20;
  std::vector<conduit::float64> extdata(extdatasize);
  gc->createView("ext20", FLOAT64_ID, extdatasize)->setExternalDataPtr(&extdata[0]);
  for(int i = 0; i < extdatasize; ++i)
  {
    extdata[i] = (conduit::float64)(i / 3.0);
  }

  gc->createView("empty_view");

  std::string file_path = file_path_base + "sidre_layout_json";
  EXPECT_TRUE(ds.getRoot()->save(file_path, "sidre_layout_json"));
  file_path = file_path_base + "conduit_layout_json";
  EXPECT_TRUE(ds.getRoot()->save(file_path, "conduit_layout_json"));

  //restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
TEST(sidre_group, import_conduit)
{
  conduit::Node input;

  input["fields/a/i0"] = (conduit::int64)(100);
  input["fields/a/d0"] = (conduit::float64)(3000.00);
  input["fields/b/s0"] = "foo";

  int ndata = 10;
  std::vector<conduit::int64> ivec(ndata);
  input["fields/c/int10"].set(&ivec[0], ndata);
  conduit::int64_array iarray = input["fields/c/int10"].as_int64_array();
  for(int i = 0; i < ndata; ++i)
  {
    iarray[i] = (conduit::int64)i;
  }

  DataStore ds;

  EXPECT_TRUE(ds.getRoot()->importConduitTree(input));

  EXPECT_EQ(ds.getRoot()->getView("fields/a/i0")->getData<conduit::int64>(), 100);
  EXPECT_NEAR(ds.getRoot()->getView("fields/a/d0")->getData<conduit::float64>(), 3000.00, 1e-12);
  EXPECT_EQ(ds.getRoot()->getView("fields/b/s0")->getString(), std::string("foo"));

  conduit::int64* sidre_data_ptr = ds.getRoot()->getView("fields/c/int10")->getData();
  for(int j = 0; j < ndata; j++)
  {
    EXPECT_EQ(iarray[j], sidre_data_ptr[j]);
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, import_conduit_external)
{
  conduit::Node input;

  input["fields/a/i0"] = (conduit::int64)(100);
  input["fields/a/d0"] = (conduit::float64)(3000.00);
  input["fields/b/s0"] = "foo";

  int ndata = 10;
  std::vector<conduit::int64> ivec(ndata);
  input["fields/c/int10"].set(&ivec[0], ndata);
  conduit::int64_array iarray = input["fields/c/int10"].as_int64_array();
  for(int i = 0; i < ndata; ++i)
  {
    iarray[i] = (conduit::int64)i;
  }

  conduit::Node& list_item0 = input["list"].append();
  list_item0.set((conduit::int64)(12));
  conduit::Node& list_item1 = input["list"].append();

  std::vector<conduit::float64> fvec(ndata);
  list_item1.set(&fvec[0], ndata);
  conduit::float64_array farray = list_item1.as_float64_array();
  for(int i = 0; i < ndata; ++i)
  {
    farray[i] = (conduit::float64)(-i);
  }

  DataStore ds;

  //Zero copy of array data
  EXPECT_TRUE(ds.getRoot()->importConduitTreeExternal(input));

  EXPECT_EQ(ds.getRoot()->getView("fields/a/i0")->getData<conduit::int64>(), 100);
  EXPECT_NEAR(ds.getRoot()->getView("fields/a/d0")->getData<conduit::float64>(), 3000.00, 1e-12);
  EXPECT_EQ(ds.getRoot()->getView("fields/b/s0")->getString(), std::string("foo"));

  //Scalar and string Views are never external.
  EXPECT_FALSE(ds.getRoot()->getView("fields/a/i0")->isExternal());
  EXPECT_FALSE(ds.getRoot()->getView("fields/a/d0")->isExternal());
  EXPECT_FALSE(ds.getRoot()->getView("fields/b/s0")->isExternal());

  //The View holding an array is external.
  EXPECT_TRUE(ds.getRoot()->getView("fields/c/int10")->isExternal());

  conduit::int64* sidre_data_ptr = ds.getRoot()->getView("fields/c/int10")->getData();
  for(int j = 0; j < ndata; j++)
  {
    EXPECT_EQ(iarray[j], sidre_data_ptr[j]);
  }

  //Change a value in the original conduit array, then test that it's changed
  //in the Sidre external view.
  if(ndata > 3)
  {
    iarray[3] += 10;
    EXPECT_EQ(iarray[3], sidre_data_ptr[3]);
  }

  //The pointers should be the same addresses as the import treated the
  //array as an external pointer.
  EXPECT_EQ((void*)sidre_data_ptr, iarray.data_ptr());

  //Test that the list group was correctly imported. The list group created
  //above has two child views, one scalar and one external array.
  Group* list = ds.getRoot()->getGroup("list");
  for(auto& view : list->views())
  {
    if(view.isScalar())
    {
      EXPECT_EQ(view.getData<conduit::int64>(), 12);
    }
    else
    {
      //If not a scalar, it must be the external array.
      EXPECT_TRUE(view.isExternal());

      conduit::float64* sidre_flt_ptr = view.getData();
      for(int j = 0; j < ndata; j++)
      {
        EXPECT_NEAR(farray[j], sidre_flt_ptr[j], 1.0e-12);
      }

      //Change a value in the original conduit array, then test that it's
      //changed in the Sidre external view.
      if(ndata > 5)
      {
        farray[5] *= 5.9;
        EXPECT_NEAR(farray[5], sidre_flt_ptr[5], 1.0e-12);
      }

      //Check address equality, as the external view points to the storage of
      //the original array.
      EXPECT_EQ((void*)sidre_flt_ptr, farray.data_ptr());
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, import_conduit_lists)
{
  conduit::Node input;

  //Give input some child nodes that are not lists, some that are.
  input["fields/a/i0"] = (conduit::int64)(100);
  input["fields/a/d0"] = (conduit::float64)(3000.00);
  input["fields/b/s0"] = "foo";

  int ndata = 10;
  std::vector<conduit::int64> ivec(ndata);
  input["fields/c/int10"].set(&ivec[0], ndata);
  conduit::int64_array iarray = input["fields/c/int10"].as_int64_array();
  for(int i = 0; i < ndata; ++i)
  {
    iarray[i] = (conduit::int64)i;
  }

  conduit::Node& list_item0 = input["list"].append();
  list_item0.set((conduit::int64)(12));
  conduit::Node& list_item1 = input["list"].append();
  list_item1.set((conduit::float64)(75.75));
  conduit::Node& list_item2 = input["list"].append();
  list_item2.set("test_str");
  conduit::Node& list_item3 = input["list"].append();
  list_item3["val1"] = (conduit::int64)(2);
  list_item3["val2"] = (conduit::float64)(4.0);

  DataStore ds;

  EXPECT_TRUE(ds.getRoot()->importConduitTree(input));

  {
    EXPECT_EQ(ds.getRoot()->getView("fields/a/i0")->getData<conduit::int64>(), 100);
    EXPECT_NEAR(ds.getRoot()->getView("fields/a/d0")->getData<conduit::float64>(), 3000.00, 1e-12);
    EXPECT_EQ(ds.getRoot()->getView("fields/b/s0")->getString(), std::string("foo"));

    conduit::int64* sidre_data_ptr = ds.getRoot()->getView("fields/c/int10")->getData();

    for(int j = 0; j < ndata; j++)
    {
      EXPECT_EQ(iarray[j], sidre_data_ptr[j]);
    }

    Group* list = ds.getRoot()->getGroup("list");
    for(IndexType idx = list->getFirstValidGroupIndex(); indexIsValid(idx);
        idx = list->getNextValidGroupIndex(idx))
    {
      Group* child = list->getGroup(idx);
      EXPECT_EQ(child->getView("val1")->getData<conduit::int64>(), 2);
      EXPECT_NEAR(child->getView("val2")->getData<conduit::float64>(), 4.0, 1e-12);
    }

    for(IndexType idx = list->getFirstValidViewIndex(); indexIsValid(idx);
        idx = list->getNextValidViewIndex(idx))
    {
      View* view = list->getView(idx);
      const conduit::Schema& schema = view->getSchema();
      if(schema.dtype().is_int64())
      {
        EXPECT_EQ(view->getData<conduit::int64>(), 12);
      }
      if(schema.dtype().is_float64())
      {
        EXPECT_NEAR(view->getData<conduit::float64>(), 75.75, 1e-12);
      }
      if(schema.dtype().is_string())
      {
        EXPECT_EQ(view->getString(), std::string("test_str"));
      }
    }
  }

#if defined(AXOM_USE_HDF5)
  std::string lists_file = "lists.hdf5";
  std::string lists_protocol = "sidre_hdf5";
#else
  std::string lists_file = "lists.json";
  std::string lists_protocol = "sidre_json";
#endif

  EXPECT_TRUE(ds.getRoot()->save(lists_file, lists_protocol));

  DataStore load_ds;
  EXPECT_TRUE(load_ds.getRoot()->load(lists_file, lists_protocol));

  {
    EXPECT_EQ(load_ds.getRoot()->getView("fields/a/i0")->getData<conduit::int64>(), 100);
    EXPECT_NEAR(load_ds.getRoot()->getView("fields/a/d0")->getData<conduit::float64>(), 3000.00, 1e-12);
    EXPECT_EQ(load_ds.getRoot()->getView("fields/b/s0")->getString(), std::string("foo"));

    conduit::int64* sidre_data_ptr = load_ds.getRoot()->getView("fields/c/int10")->getData();

    for(int j = 0; j < ndata; j++)
    {
      EXPECT_EQ(iarray[j], sidre_data_ptr[j]);
    }

    Group* list = load_ds.getRoot()->getGroup("list");
    for(IndexType idx = list->getFirstValidGroupIndex(); indexIsValid(idx);
        idx = list->getNextValidGroupIndex(idx))
    {
      Group* child = list->getGroup(idx);
      EXPECT_EQ(child->getView("val1")->getData<conduit::int64>(), 2);
      EXPECT_NEAR(child->getView("val2")->getData<conduit::float64>(), 4.0, 1e-12);
    }

    for(IndexType idx = list->getFirstValidViewIndex(); indexIsValid(idx);
        idx = list->getNextValidViewIndex(idx))
    {
      View* view = list->getView(idx);
      const conduit::Schema& schema = view->getSchema();
      if(schema.dtype().is_int64())
      {
        EXPECT_EQ(view->getData<conduit::int64>(), 12);
      }
      if(schema.dtype().is_float64())
      {
        EXPECT_NEAR(view->getData<conduit::float64>(), 75.75, 1e-12);
      }
      if(schema.dtype().is_string())
      {
        EXPECT_EQ(view->getString(), std::string("test_str"));
      }
    }
  }
}

//------------------------------------------------------------------------------

inline int pointerToAllocatorID(const void* ptr)
{
#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  if(rm.hasAllocator(const_cast<void*>(ptr)))
  {
    umpire::Allocator allocator = rm.getAllocator(const_cast<void*>(ptr));
    return allocator.getId();
  }
#endif
  AXOM_UNUSED_VAR(ptr);
  return axom::getDefaultAllocatorID();
}

TEST(sidre_group, import_conduit_into_mem_space)
{
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_DEVICE)
  const int hostAllocId = axom::detail::getAllocatorID<axom::MemorySpace::Host>();
  const int devAllocId = axom::detail::getAllocatorID<axom::MemorySpace::Device>();
#else
  // Memory space testing is trivial without device memory.
  const int hostAllocId = axom::detail::getAllocatorID<axom::MemorySpace::Dynamic>();
  const int devAllocId = hostAllocId;
#endif

  DataStore ds;
  Group* root = ds.getRoot();

  constexpr IndexType N = 5;
  int32_t* devArray = axom::allocate<std::int32_t>(N, devAllocId);
  auto dtype = conduit::DataType::int32(N);

  conduit::Node node;
  node["origOnHost"].set_dtype(dtype);
  node["origOnDev"].set_external(dtype, devArray);
  EXPECT_EQ(pointerToAllocatorID(node["origOnHost"].data_ptr()), hostAllocId);
  EXPECT_EQ(pointerToAllocatorID(node["origOnDev"].data_ptr()), devAllocId);

  /*
    Regardless of where the source data is, the destination
    data should be in the destination group's allocator id.
  */

  {
    // Import host and device data into host memory
    const int destAllocId = hostAllocId;
    Group* destGroup = root->createGroup("destGroupWithHostData");
    destGroup->setDefaultAllocator(destAllocId);

    destGroup->importConduitTree(node);
    EXPECT_EQ(destGroup->getDefaultAllocatorID(), destAllocId);
    EXPECT_TRUE(destGroup->hasView("origOnHost"));
    EXPECT_TRUE(destGroup->hasView("origOnDev"));

    const auto* viewWithHostData = destGroup->getView("origOnHost");
    int allocIdForHostData = pointerToAllocatorID(viewWithHostData->getVoidPtr());
    EXPECT_EQ(allocIdForHostData, destAllocId);

    const auto* viewWithDevData = destGroup->getView("origOnDev");
    int allocIdForDevData = pointerToAllocatorID(viewWithDevData->getVoidPtr());
    EXPECT_EQ(allocIdForDevData, destAllocId);
  }

  {
    // Import host and device data into device memory
    const int destAllocId = devAllocId;
    Group* destGroup = root->createGroup("destGroupWithDeviceData");
    destGroup->setDefaultAllocator(destAllocId);

    destGroup->importConduitTree(node);
    EXPECT_EQ(destGroup->getDefaultAllocatorID(), destAllocId);
    EXPECT_TRUE(destGroup->hasView("origOnHost"));
    EXPECT_TRUE(destGroup->hasView("origOnDev"));

    const auto* viewWithHostData = destGroup->getView("origOnHost");
    int allocIdForHostData = pointerToAllocatorID(viewWithHostData->getVoidPtr());
    EXPECT_EQ(allocIdForHostData, destAllocId);

    const auto* viewWithDevData = destGroup->getView("origOnDev");
    int allocIdForDevData = pointerToAllocatorID(viewWithDevData->getVoidPtr());
    EXPECT_EQ(allocIdForDevData, destAllocId);
  }

  axom::deallocate(devArray);
}

//------------------------------------------------------------------------------
// getDataInfo()
//------------------------------------------------------------------------------
// Local variables and methods to avoid redundant code
namespace
{
// Variables used to set check values
IndexType num_groups_chk = 0;
IndexType num_views_chk = 0;
IndexType num_views_empty_chk = 0;
IndexType num_views_buffer_chk = 0;
IndexType num_views_external_chk = 0;
IndexType num_views_scalar_chk = 0;
IndexType num_views_string_chk = 0;
IndexType num_bytes_assoc_with_views_chk = 0;
IndexType num_bytes_external_chk = 0;
IndexType num_bytes_in_buffers_chk = 0;

// Variables used to pull test values from Conduit node
IndexType num_groups = 0;
IndexType num_views = 0;
IndexType num_views_empty = 0;
IndexType num_views_buffer = 0;
IndexType num_views_external = 0;
IndexType num_views_scalar = 0;
IndexType num_views_string = 0;
IndexType num_bytes_assoc_with_views = 0;
IndexType num_bytes_external = 0;
IndexType num_bytes_in_buffers = 0;

void resetChkVals()
{
  num_groups_chk = 0;
  num_views_chk = 0;
  num_views_empty_chk = 0;
  num_views_buffer_chk = 0;
  num_views_external_chk = 0;
  num_views_scalar_chk = 0;
  num_views_string_chk = 0;
  num_bytes_assoc_with_views_chk = 0;
  num_bytes_external_chk = 0;
  num_bytes_in_buffers_chk = 0;
}

void pullTestVals(const conduit::Node& n)
{
  num_groups = n["num_groups"].value();
  num_views = n["num_views"].value();
  num_views_empty = n["num_views_empty"].value();
  num_views_buffer = n["num_views_buffer"].value();
  num_views_external = n["num_views_external"].value();
  num_views_scalar = n["num_views_scalar"].value();
  num_views_string = n["num_views_string"].value();
  num_bytes_assoc_with_views = n["num_bytes_assoc_with_views"].value();
  num_bytes_external = n["num_bytes_external"].value();
  num_bytes_in_buffers = n["num_bytes_in_buffers"].value();
}

}  // end anonymous namespace

TEST(sidre_group, get_data_info)
{
  //
  // Check 0: Ha-ha check -- empty DataStore
  //
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  conduit::Node n0;
  root->getDataInfo(n0);

  resetChkVals();
  num_groups_chk += 1;  // count root Group

  pullTestVals(n0);

  EXPECT_EQ(num_groups, num_groups_chk);
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_empty, num_views_empty_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_views_external, num_views_external_chk);
  EXPECT_EQ(num_views_scalar, num_views_scalar_chk);
  EXPECT_EQ(num_views_string, num_views_string_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_external, num_bytes_external_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  //
  // Check 1: Add "A" Group with Views. Then, non-recursive check from root
  //          Group; i.e., don't count "A" Group
  //

  // Add some buffers, 'cause it's fun....
  ds->createBuffer(INT_ID, 10)->allocate();  // NOTE: Buffer unused
                                             // in this scope
  Buffer* buff2 = ds->createBuffer(DOUBLE_ID, 10)->allocate();

  // array for external views
  int extdata[20];

  Group* gp_A = root->createGroup("A");

  View* view_A1 = gp_A->createViewAndAllocate("dat_A1", INT_ID, 5);
  View* view_A2 = gp_A->createView("dat_A2", DOUBLE_ID, 5)->attachBuffer(buff2);
  View* view_A3 = gp_A->createView("ext_A3", INT_ID, 10, &extdata);

  conduit::Node n1;
  root->getDataInfo(n1, false /* not recursive */);

  resetChkVals();
  num_groups_chk += 1;  // count root Group

  pullTestVals(n1);

  EXPECT_EQ(num_groups, num_groups_chk);
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_empty, num_views_empty_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_views_external, num_views_external_chk);
  EXPECT_EQ(num_views_scalar, num_views_scalar_chk);
  EXPECT_EQ(num_views_string, num_views_string_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_external, num_bytes_external_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  //
  // Check 2: Recursive from root Group. Count "A" Group this time.
  //

  resetChkVals();
  num_groups_chk += 1;  // count root Group

  // "A" Group and Views
  num_groups_chk += 1;
  num_views_chk = 3;         // "dat_A1", "dat_A2", "ext_A3" Views
  num_views_buffer_chk = 2;  // "dat_A1", "dat_A2" Views
  num_bytes_assoc_with_views_chk += view_A1->getTotalBytes();
  num_bytes_assoc_with_views_chk += view_A2->getTotalBytes();
  num_views_external_chk = 1;  // "ext_A3" View
  num_bytes_external_chk += view_A3->getTotalBytes();
  num_bytes_in_buffers_chk += view_A1->getTotalBytes();
  num_bytes_in_buffers_chk += buff2->getTotalBytes();  // "dat_A2" View

  conduit::Node n2;
  root->getDataInfo(n2);

  pullTestVals(n2);

  EXPECT_EQ(num_groups, num_groups_chk);
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_empty, num_views_empty_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_views_external, num_views_external_chk);
  EXPECT_EQ(num_views_scalar, num_views_scalar_chk);
  EXPECT_EQ(num_views_string, num_views_string_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_external, num_bytes_external_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  //
  // Check 3: Add "B" and "C" Groups. Recursive check from "B" Group,
  //          counting "B" and "C" Groups
  //

  resetChkVals();

  // "B" Group and Views
  Group* gp_B = root->createGroup("B");
  num_groups_chk += 1;

  View* view_B1 = gp_B->createViewAndAllocate("dat_B1", INT_ID, 5);
  num_views_chk += 1;
  num_views_buffer_chk += 1;
  num_bytes_assoc_with_views_chk += view_B1->getTotalBytes();
  num_bytes_in_buffers_chk += view_B1->getTotalBytes();

  View* view_B2 = gp_B->createView("dat_B2")->attachBuffer(buff2)->apply(DOUBLE_ID, 5, 5);
  num_views_chk += 1;
  num_views_buffer_chk += 1;
  num_bytes_assoc_with_views_chk += view_B2->getTotalBytes();
  num_bytes_in_buffers_chk += buff2->getTotalBytes();

  // "C" Group and Views
  Group* gp_C = gp_B->createGroup("C");
  num_groups_chk += 1;

  View* view_C1 = gp_C->createViewScalar("val_C1", 3);
  num_views_chk += 1;
  num_views_scalar_chk += 1;
  num_bytes_assoc_with_views_chk += view_C1->getTotalBytes();

  View* view_C2 = gp_C->createViewString("val_C2", "Homer");
  num_views_chk += 1;
  num_views_string_chk += 1;
  num_bytes_assoc_with_views_chk += view_C2->getTotalBytes();

  // "D" Group and Views
  Group* gp_D = gp_C->createGroup("D");
  num_groups_chk += 1;

  View* view_D1 = gp_D->createView("ext_D1", INT_ID, 10, extdata + 10);
  num_views_chk += 1;
  num_views_external_chk += 1;
  num_bytes_external_chk += view_D1->getTotalBytes();

  View* view_D2 = gp_D->createView("dat_D2");
  num_views_chk += 1;
  num_views_empty_chk += 1;

  conduit::Node n3;
  gp_B->getDataInfo(n3);

  pullTestVals(n3);

  EXPECT_EQ(num_groups, num_groups_chk);
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_empty, num_views_empty_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_views_external, num_views_external_chk);
  EXPECT_EQ(num_views_scalar, num_views_scalar_chk);
  EXPECT_EQ(num_views_string, num_views_string_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_external, num_bytes_external_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  //
  // Check 4: Recursive check from root Group. Count entire Group hierarchy
  //

  resetChkVals();
  num_groups_chk += 1;  // count root Group

  // "A" Group and Views
  num_groups_chk += 1;        // "A" Group
  num_views_chk += 3;         // "dat_A1", "dat_A2", "ext_A3" Views
  num_views_buffer_chk += 2;  // "dat_A1", "dat_A2" Views
  num_bytes_assoc_with_views_chk += view_A1->getTotalBytes();
  num_bytes_assoc_with_views_chk += view_A2->getTotalBytes();
  num_views_external_chk += 1;  // "ext_A3" View
  num_bytes_external_chk += view_A3->getTotalBytes();
  num_bytes_in_buffers_chk += view_A1->getTotalBytes();
  num_bytes_in_buffers_chk += buff2->getTotalBytes();  // "dat_A2" View

  // "B" Group and Views
  num_groups_chk += 1;
  num_views_chk += 2;         // "dat_B1", "dat_B2" Views
  num_views_buffer_chk += 2;  // "dat_B1", "dat_B2" Views
  num_bytes_assoc_with_views_chk += view_B1->getTotalBytes();
  num_bytes_assoc_with_views_chk += view_B2->getTotalBytes();
  num_bytes_in_buffers_chk += view_B1->getTotalBytes();
  // Note: we do not count view_B2 Buffer bytes since it shares a Buffer
  //       with view_A2 -- only count Buffer bytes once!!

  // "C" Group and Views
  num_groups_chk += 1;
  num_views_chk += 2;         // "val_C1", "val_C2" Views
  num_views_scalar_chk += 1;  // "val_C1" View
  num_views_string_chk += 1;  // "val_C2" View
  num_bytes_assoc_with_views_chk += view_C1->getTotalBytes();
  num_bytes_assoc_with_views_chk += view_C2->getTotalBytes();

  // "D" Group and Views
  num_groups_chk += 1;
  num_views_chk += 2;           // "ext_D1", "dat_D2" Views
  num_views_external_chk += 1;  // "ext_D1" View
  num_views_empty_chk += 1;     // "dat_D2" View
  num_bytes_external_chk += view_D1->getTotalBytes();

  conduit::Node n4;
  root->getDataInfo(n4);

  pullTestVals(n4);

  EXPECT_EQ(num_groups, num_groups_chk);
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_empty, num_views_empty_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_views_external, num_views_external_chk);
  EXPECT_EQ(num_views_scalar, num_views_scalar_chk);
  EXPECT_EQ(num_views_string, num_views_string_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_external, num_bytes_external_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  //
  // Miscellaneous checks to ensure expectations are met.
  //
  // NOTE: These checks only modify check data and query data values
  //       that are relevant to the tests.
  //

  //
  // Check 5: Destroy "dat_A1" View in Group "A", but leave data intact
  //          (orphaned Buffer!)
  //

  IndexType num_buffers_datastore_chk = ds->getNumBuffers();
  IndexType num_buffer_bytes_datastore_chk = ds->getTotalAllocatedBytesInBuffers();

  num_bytes_assoc_with_views_chk -= view_A1->getTotalBytes();
  num_bytes_in_buffers_chk -= view_A1->getTotalBytes();

  gp_A->destroyView("dat_A1");

  num_views_chk -= 1;
  num_views_buffer_chk -= 1;

  IndexType num_buffers_datastore = ds->getNumBuffers();
  IndexType num_buffer_bytes_datastore = ds->getTotalAllocatedBytesInBuffers();

  conduit::Node n5;
  root->getDataInfo(n5);

  pullTestVals(n5);

  // Checks for Group data query
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  // Checks for DataStore Buffers
  EXPECT_EQ(num_buffers_datastore, num_buffers_datastore_chk);
  EXPECT_EQ(num_buffer_bytes_datastore, num_buffer_bytes_datastore_chk);

  //
  // Check 6: Destroy "dat_A2" View in Group "A" and attempt to destroy its
  //          data, but its data remains because it shares its buffer with
  //          another View
  //

  num_buffers_datastore_chk = ds->getNumBuffers();
  num_buffer_bytes_datastore_chk = ds->getTotalAllocatedBytesInBuffers();

  num_bytes_assoc_with_views_chk -= view_A2->getTotalBytes();
  // Note: num_bytes_in_buffers_chk doesn't change b/c Views Buffer is
  //       still attached to another View

  gp_A->destroyViewAndData("dat_A2");

  num_views_chk -= 1;
  num_views_buffer_chk -= 1;

  num_buffers_datastore = ds->getNumBuffers();
  num_buffer_bytes_datastore = ds->getTotalAllocatedBytesInBuffers();

  conduit::Node n6;
  root->getDataInfo(n6);

  pullTestVals(n6);

  // Checks for Group data query
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  // Checks for DataStore Buffers
  EXPECT_EQ(num_buffers_datastore, num_buffers_datastore_chk);
  EXPECT_EQ(num_buffer_bytes_datastore, num_buffer_bytes_datastore_chk);

  //
  // Check 7: Describe and allocate data for empty "dat_D2" View in Group "D".
  //          Check data in subtree rooted at root is what we expect.
  //

  view_D2->allocate(DOUBLE_ID, 5);
  num_views_empty_chk -= 1;

  num_views_buffer_chk += 1;
  num_bytes_assoc_with_views_chk += view_D2->getTotalBytes();
  num_bytes_in_buffers_chk += view_D2->getTotalBytes();

  num_buffers_datastore_chk += 1;
  num_buffer_bytes_datastore_chk += view_D2->getTotalBytes();

  num_buffers_datastore = ds->getNumBuffers();
  num_buffer_bytes_datastore = ds->getTotalAllocatedBytesInBuffers();

  conduit::Node n7;
  root->getDataInfo(n7);

  pullTestVals(n7);

  // Checks for Group data query
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_empty, num_views_empty_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  // Checks for DataStore Buffers
  EXPECT_EQ(num_buffers_datastore, num_buffers_datastore_chk);
  EXPECT_EQ(num_buffer_bytes_datastore, num_buffer_bytes_datastore_chk);

  //
  // Check 8: Destroy Group "C" and check data in subtree rooted at
  //          root to make sure changes are correct.
  //

  num_groups_chk -= 2;
  num_views_chk -= 4;

  num_views_buffer_chk -= 1;
  num_bytes_assoc_with_views_chk -= view_D2->getTotalBytes();
  num_bytes_in_buffers_chk -= view_D2->getTotalBytes();

  num_views_scalar_chk -= 1;
  num_bytes_assoc_with_views_chk -= view_C1->getTotalBytes();

  num_views_string_chk -= 1;
  num_bytes_assoc_with_views_chk -= view_C2->getTotalBytes();

  num_views_external_chk -= 1;
  num_bytes_external_chk -= view_D1->getTotalBytes();

  gp_B->destroyGroup("C");

  num_buffers_datastore = ds->getNumBuffers();
  num_buffer_bytes_datastore = ds->getTotalAllocatedBytesInBuffers();

  conduit::Node n8;
  root->getDataInfo(n8);

  pullTestVals(n8);

  // Checks for Group data query
  EXPECT_EQ(num_groups, num_groups_chk);
  EXPECT_EQ(num_views, num_views_chk);
  EXPECT_EQ(num_views_empty, num_views_empty_chk);
  EXPECT_EQ(num_views_buffer, num_views_buffer_chk);
  EXPECT_EQ(num_views_external, num_views_external_chk);
  EXPECT_EQ(num_views_scalar, num_views_scalar_chk);
  EXPECT_EQ(num_views_string, num_views_string_chk);
  EXPECT_EQ(num_bytes_assoc_with_views, num_bytes_assoc_with_views_chk);
  EXPECT_EQ(num_bytes_external, num_bytes_external_chk);
  EXPECT_EQ(num_bytes_in_buffers, num_bytes_in_buffers_chk);

  // Checks for DataStore Buffers
  EXPECT_EQ(num_buffers_datastore, num_buffers_datastore_chk);
  EXPECT_EQ(num_buffer_bytes_datastore, num_buffer_bytes_datastore_chk);

  delete ds;
}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_UMPIRE

class UmpireTest : public ::testing::TestWithParam<int>
{
public:
  void SetUp() override
  {
    allocID = GetParam();
    root = ds.getRoot();
  }

  void TearDown() override { axom::setDefaultAllocator(umpire::resource::Host); }

  static constexpr int SIZE = 100;
  DataStore ds;
  Group* root;
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  int allocID;
};

TEST_P(UmpireTest, root_default_allocator)
{
  axom::setDefaultAllocator(allocID);

  DataStore dsPrime;
  ASSERT_EQ(dsPrime.getRoot()->getDefaultAllocatorID(), allocID);
}

TEST_P(UmpireTest, get_set_allocator)
{
  int defaultAllocatorID = axom::getDefaultAllocatorID();
  ASSERT_EQ(root->getDefaultAllocator().getId(), defaultAllocatorID);
  ASSERT_EQ(root->getDefaultAllocatorID(), defaultAllocatorID);

  root->setDefaultAllocator(allocID);
  defaultAllocatorID = axom::getDefaultAllocatorID();
  ASSERT_EQ(root->getDefaultAllocatorID(), allocID);

  root->setDefaultAllocator(rm.getInstance().getAllocator(defaultAllocatorID));
  ASSERT_EQ(root->getDefaultAllocator().getId(), defaultAllocatorID);
  ASSERT_EQ(root->getDefaultAllocatorID(), defaultAllocatorID);
}

//------------------------------------------------------------------------------
TEST_P(UmpireTest, allocate)
{
  {
    View* view = root->createViewAndAllocate("v", INT_ID, SIZE, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    IndexType shape[] = {1, SIZE, 1};
    View* view = root->createViewWithShapeAndAllocate("v", INT_ID, 3, shape, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    DataType dtype = conduit::DataType::default_dtype(INT_ID);
    dtype.set_number_of_elements(SIZE);
    View* view = root->createViewAndAllocate("v", dtype, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }
}

//------------------------------------------------------------------------------
TEST_P(UmpireTest, allocate_default)
{
  root->setDefaultAllocator(allocID);

  {
    View* view = root->createViewAndAllocate("v", INT_ID, SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    IndexType shape[] = {1, SIZE, 1};
    View* view = root->createViewWithShapeAndAllocate("v", INT_ID, 3, shape);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    DataType dtype = conduit::DataType::default_dtype(INT_ID);
    dtype.set_number_of_elements(SIZE);
    View* view = root->createViewAndAllocate("v", dtype);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }
}

const int allocators[] = {axom::getUmpireResourceAllocatorID(umpire::resource::Host)
  #ifdef AXOM_USE_GPU

    #ifdef UMPIRE_ENABLE_PINNED
                            ,
                          axom::getUmpireResourceAllocatorID(umpire::resource::Pinned)
    #endif

    #ifdef UMPIRE_ENABLE_DEVICE
                            ,
                          axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    #endif

    #ifdef UMPIRE_ENABLE_CONST
                            ,
                          axom::getUmpireResourceAllocatorID(umpire::resource::Constant)
    #endif

    #ifdef UMPIRE_ENABLE_UM
                            ,
                          axom::getUmpireResourceAllocatorID(umpire::resource::Unified)
    #endif

  #endif /* defined(AXOM_USE_GPU) */
};

INSTANTIATE_TEST_SUITE_P(sidre_group, UmpireTest, ::testing::ValuesIn(allocators));

#endif  // AXOM_USE_UMPIRE
