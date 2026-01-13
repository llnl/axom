# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import numpy as np

if pysidre.AXOM_USE_HDF5:
	NPROTOCOLS = 3
	PROTOCOLS = ["sidre_json", "sidre_hdf5", "json"]
else:
	NPROTOCOLS = 2
	PROTOCOLS = ["sidre_json", "json"]


# ------------------------------------------------------------------------------
# getName()
# ------------------------------------------------------------------------------
def test_get_name():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("test")

	assert grp.getName() == "test"

	grp2 = root.getGroup("foo")
	assert grp2 == None


# ------------------------------------------------------------------------------
# getPath(), getPathName()
# ------------------------------------------------------------------------------
def test_get_path_name():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	group = root.createGroup("test/a/b/c")
	grp2 = root.getGroup("test/a")
	grp3 = root.getGroup("test")

	assert root.getName() == ""
	assert root.getPath() == ""
	assert root.getPathName() == ""

	assert grp2.getName() == "a"
	assert grp2.getPath() == "test"
	assert grp2.getPathName() == "test/a"

	assert grp3.getName() == "test"
	assert grp3.getPath() == ""
	assert grp3.getPathName() == "test"

	assert group.getName() == "c"
	assert group.getPath() == "test/a/b"
	assert group.getPathName() == "test/a/b/c"


#------------------------------------------------------------------------------
# createGroup(), getGroup(), hasGroup()  with path strings
#------------------------------------------------------------------------------
def test_group_with_path():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# Test full path access when building incrementally
	group = root.createGroup("test1").createGroup("test2").createGroup("test3")
	group2 = root.getGroup("test1/test2/test3")

	assert group2 is not None
	assert group == group2

	# Test incremental access when building full path
	groupP = root.createGroup("testA/testB/testC")
	groupP2 = root.getGroup("testA").getGroup("testB").getGroup("testC")

	assert groupP2 is not None
	assert groupP == groupP2
	# test non-const getGroup() with path
	groupPParent = root.getGroup("testA/testB")
	assert groupP.getParent() == groupPParent
	assert groupP.getParent().getName() == "testB"

	# Now verify that code will not create missing groups.
	root.createGroup("testa").createGroup("testb").createGroup("testc")
	group_bada = root.getGroup("BAD/testb/testc")
	group_badb = root.getGroup("testa/BAD/testc")
	group_badc = root.getGroup("testa/testb/BAD")

	assert group_bada is None
	assert group_badb is None
	assert group_badc is None

	# Test hasGroup with paths.
	assert not root.hasGroup("BAD/testb/testc")
	assert not root.hasGroup("testa/BAD/testc")
	assert not root.hasGroup("testa/testb/BAD")

	assert root.hasGroup("test1")
	assert root.hasGroup("test1/test2")
	assert root.hasGroup("test1/test2/test3")
	group_testa = root.getGroup("testa")
	assert group_testa.hasGroup("testb")
	assert group_testa.hasGroup("testb/testc")
	assert not group_testa.hasGroup("testb/BAD")
	assert not group_testa.hasGroup("testb/testc/BAD")

	assert root.getNumGroups() == 3
	assert root.hasGroup(0)
	assert root.hasGroup(1)
	assert root.hasGroup(2)
	assert not root.hasGroup(3)
	assert not root.hasGroup(pysidre.InvalidIndex)

	testbnumgroups = group_testa.getGroup("testb").getNumGroups()
	group_cdup = group_testa.createGroup("testb/testc")

	assert group_cdup is None
	assert group_testa.getGroup("testb").getNumGroups() == testbnumgroups


#------------------------------------------------------------------------------
# createGroup(), destroyGroup()  with path strings
#------------------------------------------------------------------------------
def test_destroy_group_with_path():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# Test full path access when building incrementally
	group = root.createGroup("test1/test2/test3")

	exp_no_groups = 0
	exp_one_group = 1

	assert root.getNumGroups() == exp_one_group
	assert root.getGroup("test1").getNumGroups() == exp_one_group
	assert root.getGroup("test1/test2").getNumGroups() == exp_one_group
	assert root.getGroup("test1/test2/test3").getNumGroups() == exp_no_groups

	root.destroyGroup("test1/test2")

	assert root.getNumGroups() == exp_one_group
	assert root.getGroup("test1").getNumGroups() == exp_no_groups
	assert not root.hasGroup("test1/test2/test3")
	assert not root.hasGroup("test1/test2")

	root.destroyGroup("test1/BAD")

	assert root.getNumGroups() == exp_one_group
	assert root.getGroup("test1").getNumGroups() == exp_no_groups


# ------------------------------------------------------------------------------
# Verify getParent()
# ------------------------------------------------------------------------------
def test_get_parent():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	parent = root.createGroup("parent")
	child = parent.createGroup("child")

	assert child.getParent() == parent


# ------------------------------------------------------------------------------
# Verify getDataStore()
# ------------------------------------------------------------------------------
def test_get_datastore():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("parent")

	assert grp.getDataStore() == ds

	other_ds = grp.getDataStore()
	assert other_ds == ds


# ------------------------------------------------------------------------------
# Verify getGroup()
# ------------------------------------------------------------------------------
def test_get_group():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	parent = root.createGroup("parent")
	child = parent.createGroup("child")
	assert child.getParent() == parent

	child1 = parent.getGroup("child")
	assert child == child1

	child2 = parent.getGroup(0)
	assert child == child2

	# Check error condition
	errgrp = parent.getGroup("non-existent group")
	assert errgrp == None


# ------------------------------------------------------------------------------
# getView()
# ------------------------------------------------------------------------------
def test_get_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	parent = root.createGroup("parent")
	view1 = parent.createView("view")

	view2 = parent.getView("view")
	assert view1 == view2

	view3 = parent.getView(0)
	assert view1 == view3

	view2 = parent.getView("non-existant view")
	assert view2 == None


#------------------------------------------------------------------------------
# createView, hasView(), getView(), destroyView() with path strings
#------------------------------------------------------------------------------
def test_view_with_path():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# Test with full path access when building incrementally
	view = root.createGroup("group1").createGroup("group2").createView("view1")
	view2 = root.getView("group1/group2/view1")

	assert view2 is not None
	assert view == view2

	# Test incremental access when building with full path
	viewP = root.createView("groupA/groupB/viewA")
	viewP2 = root.getGroup("groupA").getGroup("groupB").getView("viewA")

	assert viewP2 is not None
	assert viewP == viewP2

	# Now verify that bad paths just return None, and don't create missing groups
	v_bad1 = root.getView("BAD/groupB/viewA")
	v_bad2 = root.getView("groupA/BAD/viewA")
	v_bad3 = root.getView("groupA/groupB/BAD")

	assert v_bad1 is None
	assert v_bad2 is None
	assert v_bad3 is None

	exp_no_groups = 0
	exp_one_group = 1
	exp_two_group = 2

	assert root.getNumGroups() == exp_two_group
	assert root.hasGroup("group1")
	assert root.hasGroup("groupA")
	assert root.getGroup("group1").getNumGroups() == exp_one_group
	assert root.hasGroup("group1/group2")
	assert root.getGroup("group1/group2").getNumGroups() == exp_no_groups
	assert root.getGroup("group1/group2").getNumViews() == exp_one_group
	assert root.getGroup("group1/group2").getView("view1") == view
	assert root.getGroup("group1").getView("group2/view1") == view

	assert root.getGroup("groupA").getNumGroups() == exp_one_group
	assert root.hasGroup("groupA/groupB")
	assert root.getGroup("groupA/groupB").getNumGroups() == exp_no_groups
	assert root.getGroup("groupA/groupB").getNumViews() == exp_one_group
	assert root.getGroup("groupA/groupB").getView("viewA") == viewP
	assert root.getGroup("groupA").getView("groupB/viewA") == viewP

	root.destroyView("group1/group2/view1")

	assert root.getGroup("group1/group2").getNumViews() == exp_no_groups
	assert not root.getGroup("group1/group2").hasView("view1")
	assert root.getGroup("group1/group2").getView("view1") is None
	assert not root.hasView("group1/group2/view1")
	assert root.getView("group1/group2/view1") is None

	groupA = root.getGroup("groupA")
	assert groupA.hasView("groupB/viewA")
	assert groupA.getView("groupB/viewA") == viewP
	assert root.hasView("groupA/groupB/viewA")
	assert root.getView("groupA/groupB/viewA") == viewP

	groupA.destroyView("groupB/viewA")

	assert groupA.getGroup("groupB").getNumViews() == exp_no_groups
	assert not groupA.getGroup("groupB").hasView("viewA")
	assert groupA.getGroup("groupB").getView("viewA") is None
	assert not groupA.hasView("groupB/viewA")
	assert groupA.getView("groupB/viewA") is None
	assert root.getView("groupA/groupB/viewA") is None


#------------------------------------------------------------------------------
# Verify getViewName() and getViewIndex()
#------------------------------------------------------------------------------
def test_get_view_name_index():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	parent = root.createGroup("parent")
	view1 = parent.createView("view1")
	view2 = parent.createView("view2")

	assert parent.getNumViews() == 2

	idx1 = parent.getViewIndex("view1")
	idx2 = parent.getViewIndex("view2")

	name1 = parent.getViewName(idx1)
	name2 = parent.getViewName(idx2)

	assert name1 == "view1"
	assert view1.getName() == name1

	assert name2 == "view2"
	assert view2.getName() == name2

	idx3 = parent.getViewIndex("view3")
	assert idx3 == pysidre.InvalidIndex

	name3 = parent.getViewName(idx3)
	assert name3 == ""
	assert not pysidre.nameIsValid(name3)


#------------------------------------------------------------------------------
# Verify getFirstValidGroupIndex() and getNextValidGroupIndex()
#------------------------------------------------------------------------------
def test_get_first_and_next_group_index():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	parent = root.createGroup("parent")
	group1 = parent.createGroup("group1")
	group2 = parent.createGroup("group2")
	assert parent.getNumGroups() == 2

	idx1 = parent.getFirstValidGroupIndex()
	idx2 = parent.getNextValidGroupIndex(idx1)
	idx3 = parent.getNextValidGroupIndex(idx2)

	assert idx1 == 0
	assert idx2 == 1
	assert idx3 == pysidre.InvalidIndex

	group1out = parent.getGroup(idx1)
	group2out = parent.getGroup(idx2)

	assert group1 == group1out
	assert group2 == group2out

	# Check error conditions
	emptygrp = root.createGroup("emptyGroup")
	badidx1 = emptygrp.getFirstValidGroupIndex()
	badidx2 = emptygrp.getNextValidGroupIndex(badidx1)

	assert badidx1 == pysidre.InvalidIndex
	assert badidx2 == pysidre.InvalidIndex


#------------------------------------------------------------------------------
# Verify Groups holding items in the list format
#------------------------------------------------------------------------------
def test_child_lists():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# parent is a Group in list format.
	parent = root.createGroup("parent", True)

	# Create 10 unnamed Groups as children of parent.
	for i in range(10):
		unnamed_group = parent.createUnnamedGroup()
		unnamed_group.createViewScalar("val", i)

	# Create 15 unnamed Views as children of parent.
	for i in range(15):
		if i % 3 == 0:
			unnamed_view = parent.createView("")
		elif i % 3 == 1:
			unnamed_view = parent.createViewScalar("", i * i)
		else:
			unnamed_view = parent.createViewString("", "foo")
		if not unnamed_view.isApplied():
			unnamed_view.apply(pysidre.TypeID.INT_ID, i)
			unnamed_view.allocate(pysidre.TypeID.INT_ID, i)
			vdata = unnamed_view.getDataArray()  # Returns numpy array
			for j in range(i):
				vdata[j] = j + 3

	# Create Group not in list format, show that it can't create unnamed children.
	not_list = root.createGroup("not_list", False)
	dummy_group = not_list.createUnnamedGroup()
	dummy_view = not_list.createView("")
	assert not_list.isUsingMap()
	assert not_list.getNumGroups() == 0
	assert not_list.getNumViews() == 0
	assert dummy_group is None
	assert dummy_view is None

	# Access data from unnamed Groups held by parent.
	scalars = set()
	idx = parent.getFirstValidGroupIndex()
	while pysidre.indexIsValid(idx):
		unnamed_group = parent.getGroup(idx)
		val_view = unnamed_group.getView("val")
		val = val_view.getDataInt()
		assert val >= 0 and val < 10
		scalars.add(val)
		idx = parent.getNextValidGroupIndex(idx)

	assert parent.isUsingList()
	assert len(scalars) == 10
	assert parent.hasGroup(6)
	assert not parent.hasGroup(20)

	# Destroy five of the unnamed Groups held by parent.
	idx = parent.getFirstValidGroupIndex()
	while pysidre.indexIsValid(idx):
		if idx % 2 == 1:
			parent.destroyGroup(idx)
		idx = parent.getNextValidGroupIndex(idx)

	# Add one more unnamed Group, so there should be six child Groups.
	parent.createUnnamedGroup()
	assert parent.getNumGroups() == 6

	# Access data from the unnamed Views.
	idx = parent.getFirstValidViewIndex()
	while pysidre.indexIsValid(idx):
		unnamed_view = parent.getView(idx)
		if idx % 3 == 0:
			assert unnamed_view.getTypeID() == pysidre.TypeID.INT32_ID
			num_elems = unnamed_view.getNumElements()
			assert num_elems == idx
			vdata = unnamed_view.getDataArray()
			for j in range(num_elems):
				assert vdata[j] == j + 3
		elif idx % 3 == 1:
			assert unnamed_view.isScalar()
			val = unnamed_view.getDataInt()
			assert val == idx * idx
		else:
			assert unnamed_view.isString()
			vstr = unnamed_view.getString()
			assert vstr == "foo"
		idx = parent.getNextValidViewIndex(idx)

	root.destroyGroup("parent")


#------------------------------------------------------------------------------
# Verify results with various path arguments for items in list
#------------------------------------------------------------------------------
def test_list_item_names():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# Create a group that uses the list format.
	list_test = root.createGroup("list_test", True)

	# It is recommended that all items held in a Group that uses the list
	# format be unnamed, as the names are not useful to access those items
	# from the parent Group.  Nonetheless it is allowed, and this test
	# verifies that the names are created as expected.

	# Test a group created with createUnnamedGroup, the recommended way to
	# create groups for a list.
	unnamed_group = list_test.createUnnamedGroup()
	assert unnamed_group.getName() == ""

	# Test groups created with empty string, a simple name, and a path.

	# The API says that groups cannot be created with an empty string argument,
	# so None is returned for this one.
	blank_group = list_test.createGroup("")

	# A simple name with no path syntax will be assigned to the group.
	named_group = list_test.createGroup("named")

	# With a path, the leading names are ignored and the final name is used.
	path_group_a = list_test.createGroup("testing/path")
	path_group_b = list_test.createGroup("testing/longer/path")
	found_group_a = list_test.getGroup("testing/path")
	found_group_b = list_test.getGroup("testing/longer/path")

	assert blank_group is None
	assert named_group.getName() == "named"
	assert path_group_a is None
	assert path_group_b is None
	assert found_group_a is None
	assert found_group_b is None

	# Similar tests for views

	# An empty string is the recommended way to create an unnamed view for
	# a list.
	blank_view = list_test.createView("")

	# A simple name with no path syntax will be assigned to the view.
	named_view = list_test.createView("named")

	# With a path, the leading names are ignored and the final name is used.
	path_view_a = list_test.createView("testing/path")
	path_view_b = list_test.createView("testing/longer/path")
	found_view_a = list_test.getView("testing/path")
	found_view_b = list_test.getView("testing/longer/path")

	assert blank_view.getName() == ""
	assert named_view.getName() == "named"
	assert path_view_a is None
	assert path_view_b is None
	assert found_view_a is None
	assert found_view_b is None

	root.destroyGroup("list_test")


#------------------------------------------------------------------------------
# Test Group's list format for holding contents of a vector of strings.
#------------------------------------------------------------------------------
def test_string_list():
	# Round-trip test from Python list of strings to Group and back.
	ds = pysidre.DataStore()
	root = ds.getRoot()

	str_vec = [
		"This",
		"is",
		"a",
		"vector",
		"to",
		"test",
		"strings",
		"in",
		"sidre::Group's",
		"list",
		"format"
	]

	# my_strings is a Group in list format.
	use_list_collection = True
	my_strings = root.createGroup("my_strings", use_list_collection)

	# Put strings into the Group.
	for s in str_vec:
		# The first parameter will be ignored when creating a View in a Group
		# that uses list collections, so we use the empty string
		str_view = my_strings.createViewString("", s)
		assert str_view is not None
		assert str_view.isString()

	test_vec = []

	# Get strings from the Group.
	idx = my_strings.getFirstValidViewIndex()
	while pysidre.indexIsValid(idx):
		str_view = my_strings.getView(idx)
		assert str_view is not None
		assert str_view.isString()
		vstr = str_view.getString()
		test_vec.append(vstr)
		idx = my_strings.getNextValidViewIndex(idx)

	assert str_vec == test_vec


#------------------------------------------------------------------------------
# Iterate Groups with getFirstValidGroupIndex, getNextValidGroupIndex
#------------------------------------------------------------------------------
def test_iterate_groups():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	parent = root.createGroup("parent")
	parent.createView("view1")
	parent.createView("view2")
	parent.createView("view3")
	parent.createGroup("g1")
	parent.createView("view4")
	parent.createView("view5")
	parent.createView("view6")
	parent.createGroup("g2")
	parent.createGroup("g3")
	parent.createView("view7")
	parent.createView("view8")
	parent.createView("view9")
	parent.createGroup("g4")

	groupcount = 0
	idx = parent.getFirstValidGroupIndex()
	while pysidre.indexIsValid(idx):
		groupcount += 1
		idx = parent.getNextValidGroupIndex(idx)
	assert groupcount == 4


#------------------------------------------------------------------------------
# Iterate Groups with groups()
#------------------------------------------------------------------------------
def test_iterate_groups_with_iterator():

	ds = pysidre.DataStore()
	foo_group = ds.getRoot().createGroup("foo")
	foo_group.createGroup("bar_group")
	foo_group.createGroup("bar_group/child_1")
	foo_group.createGroup("bar_group/child_2")

	foo_group.createGroup("baz_group")
	foo_group.createGroup("baz_group/child_1")
	foo_group.createGroup("baz_group/child_2")
	foo_group.createGroup("baz_group/child_3")

	foo_group.createGroup("qux_group")
	foo_group.createGroup("qux_group/child_1")

	foo_group.createView("bar_view")
	foo_group.createView("baz_view")
	foo_group.createView("qux_view")
	foo_group.createView("quux_view")

	numExpGroups = 3
	numExpViews = 4

	# iterate through groups and views of 'foo' using Python for-loop
	nGroups = 0
	for group in foo_group.groups():
		assert group.getParent() == foo_group
		assert group.getName().endswith("_group")
		nGroups += 1
	assert nGroups == numExpGroups

	nViews = 0
	for view in foo_group.views():
		assert view.getOwningGroup() == foo_group
		assert view.getName().endswith("_view")
		nViews += 1
	assert nViews == numExpViews


#------------------------------------------------------------------------------
# Verify getFirstValidViewIndex() and getNextValidIndex()
#------------------------------------------------------------------------------
def test_get_first_and_next_view_index():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	parent = root.createGroup("parent")
	view1 = parent.createView("view1")
	view2 = parent.createView("view2")

	assert parent.getNumViews() == 2

	idx1 = parent.getFirstValidViewIndex()
	idx2 = parent.getNextValidViewIndex(idx1)
	idx3 = parent.getNextValidViewIndex(idx2)
	assert idx1 == 0
	assert idx2 == 1
	assert idx3 == pysidre.InvalidIndex

	view1out = parent.getView(idx1)
	view2out = parent.getView(idx2)
	assert view1 == view1out
	assert view2 == view2out

	# Check error conditions
	emptygrp = root.createGroup("emptyGroup")
	badidx1 = emptygrp.getFirstValidViewIndex()
	badidx2 = emptygrp.getNextValidViewIndex(badidx1)

	assert badidx1 == pysidre.InvalidIndex
	assert badidx2 == pysidre.InvalidIndex


#------------------------------------------------------------------------------
# Iterate Views with getFirstValidViewIndex, getNextValidViewIndex
#------------------------------------------------------------------------------
def test_iterate_views():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	parent = root.createGroup("parent")
	parent.createView("view1")
	parent.createView("view2")
	parent.createView("view3")
	parent.createGroup("g1")
	parent.createView("view4")
	parent.createView("view5")
	parent.createView("view6")
	parent.createGroup("g2")
	parent.createGroup("g3")
	parent.createView("view7")
	parent.createView("view8")
	parent.createView("view9")
	parent.createGroup("g4")

	viewcount = 0
	idx = parent.getFirstValidViewIndex()
	while pysidre.indexIsValid(idx):
		viewcount += 1
		idx = parent.getNextValidViewIndex(idx)
	assert viewcount == 9

 
#------------------------------------------------------------------------------
# Verify getGroupName() and getGroupIndex()
#------------------------------------------------------------------------------
def test_get_group_name_index():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	parent = root.createGroup("parent")
	grp1 = parent.createGroup("grp1")
	grp2 = parent.createGroup("grp2")
	assert parent.getNumGroups() == 2

	idx1 = parent.getGroupIndex("grp1")
	idx2 = parent.getGroupIndex("grp2")

	name1 = parent.getGroupName(idx1)
	name2 = parent.getGroupName(idx2)

	assert name1 == "grp1"
	assert grp1.getName() == name1

	assert name2 == "grp2"
	assert grp2.getName() == name2

	idx3 = parent.getGroupIndex("grp3")
	assert idx3 == pysidre.InvalidIndex

	name3 = parent.getGroupName(idx3)
	assert name3 == ""
	assert not pysidre.nameIsValid(name3)


# ------------------------------------------------------------------------------
# createView()
# createViewAndAllocate()
# destroyView()
# destroyViewAndData()
# hasView()
# ------------------------------------------------------------------------------
def test_create_destroy_has_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	group = root.createGroup("parent")

	view = group.createView("view")
	assert group.getParent() == root
	assert not view.hasBuffer()
	assert group.hasView("view")
	iview = group.getViewIndex("view")
	assert iview == 0
	iview = view.getIndex()
	assert iview == 0

	# Try creating view again, should be a no-op (returns None)
	assert group.createView("view") is None

	# Create another view to make sure destroyView only destroys one view
	group.createView("viewfiller")
	assert group.getNumViews() == 2
	iviewfiller = group.getViewIndex("viewfiller")
	assert iviewfiller == 1

	group.destroyView("view")
	assert group.getNumViews() == 1
	# Check if index changed
	assert iviewfiller == group.getViewIndex("viewfiller")

	# Destroy already destroyed view. Should be a no-op, not a failure
	group.destroyView("view")
	assert group.getNumViews() == 1
	assert not group.hasView("view")

	# Try API call that specifies specific type and length
	group.createViewAndAllocate("viewWithLength1", pysidre.TypeID.INT32_ID, 50)
	iview2 = group.getViewIndex("viewWithLength1")
	assert iview == iview2  # reuse slot

	# Error condition check - try again with duplicate name, should be a no-op
	assert group.createViewAndAllocate("viewWithLength1", pysidre.TypeID.FLOAT64_ID, 50) is None
	group.destroyViewAndData("viewWithLength1")
	assert not group.hasView("viewWithLength1")

	# Should not allow negative length
	assert group.createViewAndAllocate("viewWithLengthBadLen", pysidre.TypeID.FLOAT64_ID, -1) is None

	# Try API call that specifies data type in another way
	group.createViewAndAllocate("viewWithLength2", pysidre.TypeID.FLOAT64_ID, 50)
	assert group.createViewAndAllocate("viewWithLength2", pysidre.TypeID.FLOAT64_ID, 50) is None

	# Destroy view and its buffer using index
	indx = group.getFirstValidViewIndex()
	buffer = group.getView(indx).getBuffer()
	bindx = buffer.getIndex()
	group.destroyViewAndData(indx)
	# Buffer should be destroyed from datastore
	assert ds.getBuffer(bindx) is None

	# Destroy view but not the buffer
	view = group.createViewAndAllocate("viewWithLength2", pysidre.TypeID.INT_ID, 50)
	buff = view.getBuffer()
	group.destroyView("viewWithLength2")
	assert buff.isAllocated()


#------------------------------------------------------------------------------
# createViewAndAllocate() with zero-sized array
#------------------------------------------------------------------------------
def test_create_zero_sized_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	zero_sized_view = root.createViewAndAllocate("foo", pysidre.TypeID.INT_ID, 0)
	assert zero_sized_view.isDescribed()
	assert zero_sized_view.isAllocated()


#------------------------------------------------------------------------------
# Verify createGroup(), destroyGroup(), hasGroup()
#------------------------------------------------------------------------------
def test_create_destroy_has_group():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	grp = root.createGroup("grp")
	assert grp.getParent() == root

	assert root.hasGroup("grp")

	root.destroyGroup("grp")
	assert not root.hasGroup("grp")

	grp2 = root.createGroup("grp2")
	root.destroyGroup(root.getFirstValidGroupIndex())


#------------------------------------------------------------------------------
# Test various destroy methods
#------------------------------------------------------------------------------
def test_destroy_group_and_data():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	group0 = root.createGroup("group0")
	group1 = root.createGroup("group1")
	child0 = group0.createGroup("child0")
	child1 = group0.createGroup("child1")
	child2 = group1.createGroup("child2")
	child3 = group1.createGroup("child3")
	child4 = group1.createGroup("child4")

	child0.createViewAndAllocate("intview", pysidre.TypeID.INT_ID, 15)
	foo0 = child0.createGroup("foo")
	child0.createGroup("empty")
	child0.createViewScalar("sclview", 3.14159)
	child0.createViewString("strview", "Hello world.")
	foo0.createViewAndAllocate("fooview", pysidre.TypeID.FLOAT64_ID, 12)

	int0_view = child0.getView("intview")
	int0_vals = int0_view.getDataArray()
	for i in range(15):
		int0_vals[i] = i

	fooview = foo0.getView("fooview")
	flt0_vals = fooview.getDataArray()
	for i in range(12):
		flt0_vals[i] = float(-i)

	intbuf = int0_view.getBuffer()
	fltbuf = fooview.getBuffer()

	# Store each Buffer's index for later testing.
	int_idx = intbuf.getIndex()
	flt_idx = fltbuf.getIndex()

	# Attach buffers to views in other children
	child1.createView("intview", pysidre.TypeID.INT_ID, 15, intbuf)
	foo1 = child1.createGroup("foo")
	child1.createGroup("empty")
	child1.createViewScalar("sclview", 3.14159)
	child1.createViewString("strview", "Hello world.")
	foo1.createView("fooview", pysidre.TypeID.FLOAT64_ID, 12, fltbuf)

	child2.createView("intview", pysidre.TypeID.INT_ID, 15, intbuf)
	foo2 = child2.createGroup("foo")
	child2.createGroup("empty")
	child2.createViewScalar("sclview", 3.14159)
	child2.createViewString("strview", "Hello world.")
	foo2.createView("fooview", pysidre.TypeID.FLOAT64_ID, 12, fltbuf)

	child3.createView("intview", pysidre.TypeID.INT_ID, 15, intbuf)
	foo3 = child3.createGroup("foo")
	child3.createGroup("empty")
	child3.createViewScalar("sclview", 3.14159)
	child3.createViewString("strview", "Hello world.")
	foo3.createView("fooview", pysidre.TypeID.FLOAT64_ID, 12, fltbuf)

	child4.createView("intview", pysidre.TypeID.INT_ID, 15, intbuf)
	foo4 = child4.createGroup("foo")
	child4.createGroup("empty")
	child4.createViewScalar("sclview", 3.14159)
	child4.createViewString("strview", "Hello world.")
	foo4.createView("fooview", pysidre.TypeID.FLOAT64_ID, 12, fltbuf)

	# Beginning state: 2 Buffers, each attached to 5 Views.
	assert ds.getNumBuffers() == 2
	assert intbuf.getNumViews() == 5
	assert fltbuf.getNumViews() == 5

	# Destroy "child0/foo" by path. This destroys the View that created fltbuf.
	assert child0.getNumGroups() == 2
	group0.destroyGroupAndData("child0/foo")

	# Verify that fltbuf now is attached to 4 Views and child0 has only one
	# group "empty".
	assert fltbuf.getNumViews() == 4
	assert child0.getNumGroups() == 1
	assert child0.hasGroup("empty")

	# Destroy child3 using index argument
	idx3 = group1.getGroupIndex("child3")
	group1.destroyGroupAndData(idx3)

	# intbuf and fltbuf both lose one attached View
	assert intbuf.getNumViews() == 4
	assert fltbuf.getNumViews() == 3

	# Verify Buffers' data can be accessed by other Views
	int4_vals = child4.getView("intview").getDataArray()
	for i in range(15):
		assert int4_vals[i] == i

	flt4_vals = foo4.getView("fooview").getDataArray()
	for i in range(12):
		assert abs(flt4_vals[i] - float(-i)) < 1.0e-12

	# Destroy Groups held by child1, removes "foo" and "empty" but leaves Views
	assert child1.getNumGroups() == 2
	assert child1.getNumViews() == 3
	child1.destroyGroupsAndData()
	assert child1.getNumGroups() == 0
	assert child1.getNumViews() == 3

	# That removed one more View from fltbuf but left intbuf unchanged.
	assert intbuf.getNumViews() == 4
	assert fltbuf.getNumViews() == 2

	# Destroy the entire subtree of child4
	assert child4.getNumGroups() == 2
	assert child4.getNumViews() == 3
	child4.destroyGroupSubtreeAndData()
	assert child4.getNumGroups() == 0
	assert child4.getNumViews() == 0

	# Both buffers lost one more View
	assert intbuf.getNumViews() == 3
	assert fltbuf.getNumViews() == 1

	# The View at "group1/child2/foo/fooview" is the only View still attached to fltbuf
	assert ds.getNumBuffers() == 2
	assert group1.hasView("child2/foo/fooview")
	assert group1.getView("child2/foo/fooview").hasBuffer()
	assert group1.getView("child2/foo/fooview").getBuffer() == fltbuf

	# Destroy entire subtree of group1. This detaches last View from fltbuf and destroys fltbuf
	group1.destroyGroupSubtreeAndData()
	assert ds.getNumBuffers() == 1
	assert ds.hasBuffer(int_idx)
	assert not ds.hasBuffer(flt_idx)
	assert group1.getNumViews() == 0
	assert group1.getNumGroups() == 0

	# intbuf still attached to the "intview" Views in child0 and child1
	assert intbuf.getNumViews() == 2
	assert group0.getView("child0/intview").getBuffer() == intbuf
	assert group0.getView("child1/intview").getBuffer() == intbuf

	# Destroy everything below root, remaining buffer will be destroyed
	root.destroyGroupSubtreeAndData()
	assert not ds.hasBuffer(int_idx)
	assert ds.getNumBuffers() == 0
	assert root.getNumViews() == 0
	assert root.getNumGroups() == 0


def test_group_name_collisions():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	flds = root.createGroup("fields")
	flds.createView("a")

	assert flds.hasChildView("a")

	# Attempt to create duplicate group name
	assert root.createGroup("fields") is None

	# Attempt to create duplicate view name
	assert flds.createView("a") is None

	# Attempt to create a group named the same as an existing view
	assert flds.createGroup("a") is None

	# Attempt to create a view named the same as an existing group
	assert root.createView("fields") is None

	# Create groups with unusual names/paths
	root.createGroup("here//is/path")
	root.createGroup("éch≈o/Ωd")
	root.createGroup("../group/..")

	# Print all group names
	idx = root.getFirstValidGroupIndex()
	while pysidre.indexIsValid(idx):
		print(root.getGroup(idx).getName())
		idx = root.getNextValidGroupIndex(idx)


def test_view_copy_move():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	flds = root.createGroup("fields")
	extdata = np.array([0] * 10)

	views = [None] * 6
	names = [None] * 6

	# Create views in different states
	views[0] = flds.createView("empty0")
	views[1] = flds.createView("empty1", pysidre.TypeID.INT32_ID, 10)
	views[2] = flds.createViewAndAllocate("buffer", pysidre.TypeID.INT32_ID, 10)
	views[3] = flds.createView("external", pysidre.TypeID.INT32_ID, 10)
	views[3].setExternalData(extdata)
	views[4] = flds.createViewScalar("scalar", 25)
	views[5] = flds.createViewString("string", "I am string")

	buffdata = flds.getView("buffer").getDataArray()
	for i in range(10):
		extdata[i] = i
		buffdata[i] = i + 100

	for i in range(6):
		names[i] = views[i].getName()
		assert flds.hasView(names[i])

	# Test moving a view from flds to sub1
	sub1 = flds.createGroup("sub1")

	for i in range(6):
		sub1.moveView(views[i])
		assert not flds.hasView(names[i])
		assert sub1.hasView(names[i])

		# Moving to same group is a no-op
		sub1.moveView(views[i])
		assert sub1.hasView(names[i])

	sub2 = flds.createGroup("sub2")

	for i in range(6):
		sub2.copyView(views[i])
		assert sub1.hasView(names[i])
		assert sub2.hasView(names[i])

	# Check copies
	view1 = sub1.getView("empty0")
	view2 = sub2.getView("empty0")
	assert view1 is not view2
	assert view2.isEmpty()
	assert not view2.isDescribed()
	assert not view2.isAllocated()
	assert not view2.isApplied()

	view1 = sub1.getView("empty1")
	view2 = sub2.getView("empty1")
	assert view1 is not view2
	assert view2.isEmpty()
	assert view2.isDescribed()
	assert not view2.isAllocated()
	assert not view2.isApplied()

	view1 = sub1.getView("buffer")
	view2 = sub2.getView("buffer")
	assert view1 is not view2
	assert view2.hasBuffer()
	assert view2.isDescribed()
	assert view2.isAllocated()
	assert view2.isApplied()
	assert view1.getBuffer() == view2.getBuffer()
	assert view1.getBuffer().getNumViews() == 2

	view1 = sub1.getView("external")
	view2 = sub2.getView("external")
	assert view1 is not view2
	assert view2.isExternal()
	assert view2.isDescribed()
	assert view2.isAllocated()
	assert view2.isApplied()
	assert view1.getDataInt() == view2.getDataInt()

	view1 = sub1.getView("scalar")
	view2 = sub2.getView("scalar")
	assert view1 is not view2
	assert view2.isScalar()
	assert view2.isDescribed()
	assert view2.isAllocated()
	assert view2.isApplied()
	assert view1.getDataInt() == view2.getDataInt()

	view1 = sub1.getView("string")
	view2 = sub2.getView("string")
	assert view1 is not view2
	assert view2.isString()
	assert view2.isDescribed()
	assert view2.isAllocated()
	assert view2.isApplied()
	svalue = view1.getString()
	assert svalue == "I am string"


def test_groups_move_copy():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	flds = root.createGroup("fields")

	ga = flds.createGroup("a")
	gb = flds.createGroup("b")
	gc = flds.createGroup("c")

	f0value = 100.0
	val = 101.0

	bschild = gb.createGroup("childOfB")

	ga.createView("i0").setScalar(1)
	gb.createView("f0").setScalar(f0value)
	gc.createView("d0").setScalar(3000.0)
	bschild.createView("val").setScalar(val)

	buffercount = ds.getNumBuffers()

	# Check that all sub groups exist
	assert flds.hasGroup("a")
	assert flds.hasGroup("b")
	assert flds.hasGroup("c")

	# Move "b" to a child of "sub"
	assert gb.getIndex() == 1
	assert gb.getParent() == flds
	gsub = flds.createGroup("sub")
	gb0 = gsub.moveGroup(gb)

	# gb0 is an alias to gb
	assert gb == gb0
	assert gb.getIndex() == 0
	assert gb.getParent() == gsub
	assert gb0.getNumGroups() == 1
	assert gb0.getGroup("childOfB") == bschild
	assert bschild.getNumGroups() == 0
	assert buffercount == ds.getNumBuffers()

	assert gb0.getNumViews() == 1
	assert gb0.hasChildView("f0")
	if gb0.hasChildView("f0"):
		assert gb0.getView("f0").getDataFloat() == f0value

	assert bschild.getNumViews() == 1
	assert bschild.hasChildView("val")
	if bschild.hasChildView("val"):
		assert bschild.getView("val").getDataFloat() == val

	assert flds.getNumGroups() == 3
	assert flds.hasGroup("a")
	assert flds.hasGroup("sub")
	assert flds.hasGroup("c")

	assert flds.getGroup("sub").getGroup("b") == gb

	# Verify that we can copy a group into an empty group
	containCopy = root.createGroup("containCopy")
	theCopy = containCopy.copyGroup(flds)
	assert theCopy.isEquivalentTo(flds)
	assert containCopy.getNumGroups() == 1
	assert buffercount == ds.getNumBuffers()

	# Verify that we can copy a group when there is no name clash
	anotherCopy = root.createGroup("anotherCopy")
	anotherCopy.createGroup("futureSiblingGroup")
	theOtherCopy = anotherCopy.copyGroup(flds)
	assert anotherCopy.getNumGroups() == 2
	assert theOtherCopy.isEquivalentTo(flds)
	assert buffercount == ds.getNumBuffers()

	# Verify that we cannot copy a group when there is a name clash
	otherB = containCopy.createGroup("b")
	otherB.createView("f1").setScalar(42.0)
	otherB.createGroup("Q")
	triedCopy = gsub.copyGroup(otherB)
	assert triedCopy is None
	assert gsub.getNumGroups() == 1
	assert gsub.hasChildGroup("b")
	assert buffercount == ds.getNumBuffers()

	assert gb0.getNumGroups() == 1
	assert gb0.getGroup("childOfB") == bschild
	assert gb0.getNumViews() == 1
	assert gb0.hasChildView("f0")
	if gb0.hasChildView("f0"):
		assert gb0.getView("f0").getDataFloat() == f0value


def test_group_deep_copy():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	flds = root.createGroup("fields")

	ga = flds.createGroup("a")
	gb = flds.createGroup("b")

	dval0 = 100.0
	dval1 = 301.0

	ga.createView("i0").setScalar(1)
	ga.createView("d0").setScalar(dval0)
	gb.createView("d1").setScalar(dval1)
	gb.createView("s0").setString("my string")

	# Check that all sub groups exist
	assert flds.hasGroup("a")
	assert flds.hasGroup("b")

	viewlen = 8
	ownsbuf = ga.createViewAndAllocate("ownsbuf", pysidre.TypeID.INT32_ID, viewlen)
	int_vals = ownsbuf.getDataArray()
	for i in range(viewlen):
		int_vals[i] = i + 1

	buflen = 24
	dbuff = ds.createBuffer()
	dbuff.allocate(pysidre.TypeID.FLOAT64_ID, buflen)
	buf_ptr = dbuff.getDataArray()
	for i in range(buflen):
		buf_ptr[i] = 2.0 * float(i)

	NUM_VIEWS = 4
	size = [5, 4, 10, 11]
	stride = [3, 2, 2, 1]
	offset = [2, 9, 0, 10]
	names = ["viewa", "viewb", "viewc", "viewd"]

	for i in range(NUM_VIEWS):
		ga.createView(names[i], dbuff).apply(size[i], offset[i], stride[i])

	extlen = 30
	ext_array = np.array([-1.0 * float(i) for i in range(extlen)])

	for i in range(NUM_VIEWS):
		gb.createView(names[i], ext_array).apply(pysidre.TypeID.FLOAT64_ID, size[i], offset[i], stride[i])

	deep_copy = root.createGroup("deep_copy")
	deep_copy.deepCopyGroup(flds)

	assert deep_copy.hasGroup("fields/a")
	assert deep_copy.hasGroup("fields/b")

	copy_ga = deep_copy.getGroup("fields/a")
	copy_gb = deep_copy.getGroup("fields/b")

	io_val = copy_ga.getView("i0").getDataInt()
	assert io_val == 1
	assert abs(copy_ga.getView("d0").getDataFloat() - dval0) < 1.0e-12
	assert abs(copy_gb.getView("d1").getDataFloat() - dval1) < 1.0e-12
	assert copy_gb.getView("s0").getString() == "my string"

	for i in range(NUM_VIEWS):
		assert copy_ga.hasView(names[i])
		copy_view = copy_ga.getView(names[i])
		assert copy_view.hasBuffer()

		# The deep copy creates a compact buffer in the copied View, associated only with that View.
		buffer = copy_view.getBuffer()
		assert buffer.getNumViews() == 1
		assert buffer.getNumElements() == size[i]
		assert copy_view.getOffset() == 0
		assert copy_view.getStride() == 1

		fdata = copy_view.getDataArray()
		for j in range(size[i]):
			assert abs(fdata[j] - 2.0 * (offset[i] + j * stride[i])) < 1.0e-12

	for i in range(NUM_VIEWS):
		assert copy_gb.hasView(names[i])
		copy_view = copy_gb.getView(names[i])
		assert copy_view.hasBuffer()

		# The deep copy creates a compact buffer in the copied View, associated only with that View.
		buffer = copy_view.getBuffer()
		assert buffer.getNumViews() == 1
		assert buffer.getNumElements() == size[i]
		assert copy_view.getOffset() == 0
		assert copy_view.getStride() == 1

		fdata = copy_view.getDataArray()
		for j in range(size[i]):
			assert abs(fdata[j] - (-1.0 * (offset[i] + j * stride[i]))) < 1.0e-12


def test_create_destroy_view_and_buffer2():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("grp")

	viewName1 = "viewBuffer1"
	viewName2 = "viewBuffer2"

	view1 = grp.createViewAndAllocate(viewName1, pysidre.TypeID.INT_ID, 1)
	view2 = grp.createViewAndAllocate(viewName2, pysidre.TypeID.INT_ID, 1)

	assert grp.hasView(viewName1)
	assert grp.getView(viewName1) == view1

	assert grp.hasView(viewName2)
	assert grp.getView(viewName2) == view2

	bufferId1 = view1.getBuffer().getIndex()

	grp.destroyViewAndData(viewName1)

	assert not grp.hasView(viewName1)
	assert ds.getNumBuffers() == 1

	buffer1 = ds.getBuffer(bufferId1)
	assert buffer1 is None

	view3 = grp.createView("viewBuffer3")
	grp.destroyViewsAndData()
	# should be no-op
	grp.destroyViewsAndData()


def test_create_destroy_alloc_view_and_buffer():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("grp")

	viewName1 = "viewBuffer1"
	viewName2 = "viewBuffer2"

	# Use create + alloc convenience methods
	# This one is the DataType method
	view1 = grp.createViewAndAllocate(viewName1, pysidre.TypeID.INT_ID, 10)

	assert grp.hasChildView(viewName1)
	assert grp.getView(viewName1) == view1

	v1_vals = view1.getDataArray()
	for i in range(10):
		v1_vals[i] = i

	assert view1.getNumElements() == 10
	assert view1.getTotalBytes() == 10 * 4  # Assuming int is 4 bytes

	grp.destroyViewAndData(viewName1)


def test_create_view_of_buffer_with_schema():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# Use create + alloc convenience methods
	base = root.createViewAndAllocate("base", pysidre.TypeID.INT_ID, 10)
	base_vals = base.getDataArray()
	for i in range(10):
		if i < 5:
			base_vals[i] = 10
		else:
			base_vals[i] = 20

	base_buff = base.getBuffer()

	# Create two views into this buffer
	# View for the first 5 values
	sub_a = root.createView("sub_a", base_buff)
	sub_a.apply(pysidre.TypeID.INT_ID, 5)

	sub_a_vals = sub_a.getDataArray()
	for i in range(5):
		assert sub_a_vals[i] == 10


def test_create_destroy_view_and_data():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("grp")

	view_name1 = "viewBuffer1"
	view_name2 = "viewBuffer2"

	view1 = grp.createViewAndAllocate(view_name1, pysidre.TypeID.INT32_ID, 1)
	view2 = grp.createViewAndAllocate(view_name2, pysidre.TypeID.INT32_ID, 1)

	assert grp.hasView(view_name1)
	assert grp.getView(view_name1) == view1

	assert grp.hasView(view_name2)
	assert grp.getView(view_name2) == view2

	tmpbuf = view1.getBuffer()
	bufferid1 = tmpbuf.getIndex()

	grp.destroyViewAndData(view_name1)

	assert not grp.hasView(view_name1)
	assert ds.getNumBuffers() == 1



def test_create_destroy_alloc_view_and_data():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("grp")

	view_name1 = "viewBuffer1"
	view_name2 = "viewBuffer2"

	# Use create + alloc convenience methods
	# this one is the DataType & method
	view1 = grp.createViewAndAllocate(view_name1, pysidre.TypeID.INT32_ID, 10)
	assert grp.hasView(view_name1)
	assert grp.getView(view_name1) == view1

	# TODO getData, need numpy array implementation
	v1_vals = view1.getDataArray()

	for i in range(10):
		v1_vals[i] = i

	assert view1.getNumElements() == 10
	grp.destroyViewAndData(view_name1)


def test_create_view_of_buffer_with_datatype():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# Use create + alloc convenience methods
	# this one is the DataType & method
	base = root.createViewAndAllocate("base", pysidre.TypeID.INT32_ID, 10)
	base_vals = base.getDataArray()

	base_vals[0:5] = 10
	base_vals[5:10] = 20

	base_buff = base.getBuffer()

	# Create view into this buffer
	sub_a = root.createView("sub_a", pysidre.TypeID.INT32_ID, 10, base_buff)

	sub_a_vals = root.getView("sub_a").getDataArray()

	for i in range(5):
		assert sub_a_vals[i] == 10
	for i in range(5,10):
		assert sub_a_vals[i] == 20


def test_save_restore_empty_datastore():
	file_path_base = "py_sidre_empty_datastore_"
	ds1 = pysidre.DataStore()
	root1 = ds1.getRoot()

	for i in range(NPROTOCOLS):
		file_path = file_path_base + PROTOCOLS[i]
		root1.save(file_path, PROTOCOLS[i])

	for i in range(NPROTOCOLS):
		if PROTOCOLS[i] != "sidre_hdf5":
			continue
		file_path = file_path_base + PROTOCOLS[i]
		
		ds2 = pysidre.DataStore()
		root2 = ds2.getRoot()

		root2.load(file_path, PROTOCOLS[i])

		assert ds2.getNumBuffers() == 0
		assert root2.getNumGroups() == 0
		assert root2.getNumViews() == 0


def test_save_restore_scalars_and_strings():
	file_path_base = "py_sidre_save_scalars_and_strings_"
	ds1 = pysidre.DataStore()
	root1 = ds1.getRoot()

	view = root1.createViewScalar("i0", 1)
	view = root1.createViewScalar("f0", 1.0)
	view = root1.createViewScalar("d0", 10.0)
	view = root1.createViewString("s0", "I am a string")

	for i in range(NPROTOCOLS):
		file_path = file_path_base + PROTOCOLS[i]
		root1.save(file_path, PROTOCOLS[i])

	for i in range(NPROTOCOLS):
		if PROTOCOLS[i] != "sidre_hdf5":
			continue
		file_path = file_path_base + PROTOCOLS[i]
		
		ds2 = pysidre.DataStore()
		root2 = ds2.getRoot()

		root2.load(file_path, PROTOCOLS[i])

		assert root1.isEquivalentTo(root2)

		view = root2.getView("i0")
		i0 = view.getDataInt()
		assert i0 == 1

		view = root2.getView("f0")
		f0 = view.getDataFloat()
		assert f0 == 1.0

		view = root2.getView("d0")
		d0 = view.getDataFloat()
		assert d0 == 10.0

		view = root2.getView("s0")
		s0 = view.getString()
		assert s0 == "I am a string"


def test_save_restore_external_data():
	file_path_base = "py_sidre_save_external_"

	nfoo = 10
	foo1 = np.array(range(nfoo))

	# dtype is necessary, or garbage conversion from float64 --> int64 takes place
	foo2 = np.zeros(10, dtype = int)
	foo3 = np.empty(0)
	foo4 = np.array([i+1 for i in range(10)])

	shape = np.array([10,2])
	int2d1 = np.column_stack((foo1, foo1 + nfoo))
	int2d2 = np.zeros((10,2), dtype = int)

	ds1 = pysidre.DataStore()
	root1 = ds1.getRoot()

	root1.createView("external_array", pysidre.TypeID.INT64_ID, nfoo, foo1)
	root1.createView("empty_array", pysidre.TypeID.INT64_ID, 0, foo3)
	root1.createView("external_undescribed").setExternalData(foo4)
	root1.createViewWithShape("int2d", pysidre.TypeID.INT64_ID, 2, shape, int2d1)

	for protocol in PROTOCOLS:
		file_path = file_path_base + protocol
		assert root1.save(file_path, protocol) == True

	# Now load back in
	for protocol in PROTOCOLS:
		# Only restore sidre_hdf5 protocol
		if protocol != "sidre_hdf5":
			continue
		file_path = file_path_base + protocol

		ds2 = pysidre.DataStore()
		root2 = ds2.getRoot()

		assert root2.load(file_path, protocol) == True

		# Load has the set type and size of the view.
		# Now set the external address before calling load_external
		view1 = root2.getView("external_array")
		assert view1.isExternal() == True, "external_array is external"
		assert view1.isDescribed() == True, "external_array is described"
		assert view1.getTypeID() == pysidre.TypeID.INT64_ID, "external_array get TypeId"
		assert view1.getNumElements() == nfoo, "external_array get num elements"
		view1.setExternalData(foo2)

		view2 = root2.getView("empty_array")
		assert view2.isExternal() == True, "empty_array is external"
		assert view2.isDescribed() == True, "empty_array is described"
		assert view2.getTypeID() == pysidre.TypeID.INT64_ID, "empty_array get TypeId"
		view2.setExternalData(foo3)

		view3 = root2.getView("external_undescribed")
		assert view3.isEmpty() == True, "external_undescribed is empty"
		assert view3.isDescribed() == False, "external_undescribed is undescribed"

		extents = np.zeros(7)
		view4 = root2.getView("int2d")
		assert view4.isExternal() == True, "int2d is external"
		assert view4.isDescribed() == True, "int2d is described"
		assert view4.getTypeID() == pysidre.TypeID.INT64_ID, "int2d get TypeId"
		assert view4.getNumElements() == nfoo * 2, "int2d get num elements"
		assert view4.getNumDimensions() == 2, "int2d get num dimensions"

		rank, extents = view4.getShape(7, extents)
		assert rank == 2, "int2d rank"
		assert extents[0] == nfoo
		assert extents[1] == 2
		view4.setExternalData(int2d2)

		# Read external data into views
		assert root2.loadExternalData(file_path) == True

		# Check loaded data
		assert np.array_equal(foo1,foo2), "compare foo1 foo2"

		assert np.array_equal(view1.getDataArray(), foo2)
		assert np.array_equal(view2.getDataArray(), foo3)
		assert np.array_equal(view4.getDataArray(), int2d2)

		assert np.array_equal(int2d1,int2d2)



def test_save_restore_other():
	file_path_base = "py_sidre_empty_other_"
	ndata = 10

	ds1 = pysidre.DataStore()
	root1 = ds1.getRoot()

	shape1 = np.array([ndata, 2])
	view1 = root1.createView("empty_view")
	view2 = root1.createView("empty_described", pysidre.TypeID.INT32_ID, ndata)
	view3 = root1.createViewWithShape("empty_shape", pysidre.TypeID.INT32_ID, 2, shape1)
	view4 = root1.createViewWithShapeAndAllocate("buffer_shape", pysidre.TypeID.INT32_ID, 2, shape1)

	for protocol in PROTOCOLS:
		file_path = file_path_base + protocol
		assert root1.save(file_path, protocol) == True

	# Now load back in
	for protocol in PROTOCOLS:
		# Only restore sidre_hdf5 protocol
		if protocol != "sidre_hdf5":
			continue

		file_path = file_path_base + protocol

		ds2 = pysidre.DataStore()
		root2 = ds2.getRoot()

		root2.load(file_path, protocol)

		view1 = root2.getView("empty_view")
		assert view1.isEmpty() == True, "empty_view is empty"
		assert view1.isDescribed() == False, "empty_view is described"

		view2 = root2.getView("empty_described")
		assert view2.isEmpty() == True, "empty_described is empty"
		assert view2.isDescribed() == True, "empty_described is described"
		assert view2.getTypeID() == pysidre.TypeID.INT32_ID, "empty_described get TypeID"
		assert view2.getNumElements() == ndata, "empty_described get num elements"

		view3 = root2.getView("empty_shape")
		assert view3.isEmpty() == True, "empty_shape is empty"
		assert view3.isDescribed() == True, "empty_shape is described"
		assert view3.getTypeID() == pysidre.TypeID.INT32_ID, "empty_shape get TypeID"
		assert view3.getNumElements() == ndata * 2, "empty_shape get num elements"
		shape2 = np.zeros(7)
		rank, shape2 = view3.getShape(7, shape2)
		assert rank == 2, "empty_shape rank"
		assert shape2[0] == ndata and shape2[1] == 2, "empty_shape get shape"

		view4 = root2.getView("buffer_shape")
		assert view4.hasBuffer() == True, "buffer_shape has buffer"
		assert view4.isDescribed() == True, "buffer_shape is described"
		assert view4.getTypeID() == pysidre.TypeID.INT32_ID, "buffer_shape get TypeID"
		assert view4.getNumElements() == ndata * 2, "buffer_shape get num elements"
		shape2 = np.zeros(7)
		rank, shape2 = view4.getShape(7, shape2)
		assert rank == 2, "buffer_shape rank"
		assert shape2[0] == ndata and shape2[1] == 2, "buffer_shape get shape"



def test_rename_group():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	child1 = root.createGroup("g_a")
	child2 = root.createGroup("g_b")
	child3 = root.createGroup("g_c")

	# Rename should not change the index
	assert child1.getIndex() == 0
	success = child1.rename("g_r")
	assert success
	assert child1.getName() == "g_r"
	assert child1.getIndex() == 0
	assert root.hasGroup("g_r")
	assert not root.hasGroup("g_a")

	# Try to rename to path
	success = child2.rename("fields/g_s")
	assert not success
	assert child2.getName() == "g_b"

	# Try to rename to existing group name
	success = child3.rename("g_b")
	assert not success
	assert child3.getName() == "g_c"

	# Rename root group
	assert not pysidre.indexIsValid(root.getIndex())
	assert root.getParent() == root
	assert root.getName() == ""
	root.rename("newroot")
	assert not pysidre.indexIsValid(root.getIndex())
	assert root.getParent() == root
	assert root.getName() == "newroot"


# Fortran comment - redo these, the C++ tests were heavily rewritten
def test_save_restore_simple():
	file_path = "py_out_sidre_group_save_restore_simple"
	ds = pysidre.DataStore()
	root = ds.getRoot()
	flds = root.createGroup("fields")

	ga = flds.createGroup("a")

	i0_view = ga.createViewScalar("i0", 1)

	assert root.hasGroup("fields") == True
	assert flds.hasGroup("a") == True
	assert ga.hasView("i0") == True

	root.save(file_path, "sidre_conduit_json")

	ds2 = pysidre.DataStore()
	root2 = ds2.getRoot()

	root2.load(file_path, "sidre_conduit_json")

	flds = root2.getGroup("fields")
	assert flds.hasGroup("a") == True
	ga = flds.getGroup("a")
	i0_view = ga.getView("i0")
	assert i0_view.getDataInt() == 1


# Fortran comment - redo these, the C++ tests were heavily rewritten
def test_save_restore_complex():
	file_path = "py_out_sidre_group_save_restore_complex"

	ds = pysidre.DataStore()
	root = ds.getRoot()
	flds = root.createGroup("fields")

	ga = flds.createGroup("a")
	gb = flds.createGroup("b")
	gc = flds.createGroup("c")

	ga.createViewScalar("i0", 1)
	gb.createViewScalar("f0", 100.0)
	gc.createViewScalar("d0", 3000.0)

	# Check that all sub groups exist
	assert flds.hasGroup("a") == True
	assert flds.hasGroup("b") == True
	assert flds.hasGroup("c") == True

	root.save(file_path, "sidre_conduit_json")

	ds2 = pysidre.DataStore()
	root2 = ds2.getRoot()

	root2.load(file_path, "sidre_conduit_json")

	flds = root2.getGroup("fields")

	# Check that all sub groups exist
	assert flds.hasGroup("a") == True
	assert flds.hasGroup("b") == True
	assert flds.hasGroup("c") == True

	ga = flds.getGroup("a")
	gb = flds.getGroup("b")
	gc = flds.getGroup("c")

	i0_view = ga.getView("i0")
	f0_view = gb.getView("f0")
	d0_view = gc.getView("d0")

	assert i0_view.getDataInt() == 1
	assert f0_view.getDataFloat() == 100.0
	assert d0_view.getDataFloat() == 3000.0


# Fortran - for some reason not part of main program
def test_save_load_preserve_contents():
	file_path_base0 = "py_sidre_save_preserve_contents_tree0_"
	file_path_base1 = "py_sidre_save_preserve_contents_tree1_"

	ds = pysidre.DataStore()
	root = ds.getRoot()
	tree0 = root.createGroup("tree0")

	ga = tree0.createGroup("a")
	gb = tree0.createGroup("b")
	gc = tree0.createGroup("c")

	i0_view = ga.createViewScalar("i0", 100)
	f0_view = ga.createViewScalar("f0", 3000.0)
	s0_view = gb.createViewString("s0", "foo")
	i10_view = gc.createViewAndAllocate("int10", pysidre.TypeID.INT32_ID, 10)

	v1_vals = i10_view.getDataArray()
	for i in range(10):
		v1_vals[i] = i

	for protocol in PROTOCOLS:
		# Only restore sidre_hdf5 protocol
		if protocol != "sidre_hdf5":
			continue

		file_path0 = file_path_base0 + protocol
		tree0.save(file_path0, protocol)

		tree1 = root.createGroup("tree1")

		gx = tree1.createGroup("x")
		gy = tree1.createGroup("y")
		gz = tree1.createGroup("z")

		i20_view = gx.createViewAndAllocate("int20", pysidre.TypeID.INT32_ID, 20)
		v2_vals = i20_view.getDataArray()
		for i in range(20):
			v2_vals[i] = 2 * i

		i1_view = gy.createViewScalar("i1", 400)
		f1_view = gz.createViewScalar("f1", 17.0)

		file_path1 = file_path_base1 + protocol
		assert tree1.save(file_path1, protocol)

		dsload = pysidre.DataStore()
		ldroot = dsload.getRoot()

		ldtree0 = ldroot.createGroup("tree0")
		ldtree0.load(file_path0, protocol)
		ldtree0.load(file_path1, protocol, True)
		ldtree0.rename("tree1")


		i0_view = ldroot.getView("tree1/a/i0")
		i0 = i0_view.getDataInt()
		assert i0 == 100

		f0_view = ldroot.getView("tree1/a/f0")
		f0 = f0_view.getDataFloat()
		assert f0 == 3000.0

		s0_view = ldroot.getView("tree1/b/s0")
		s0 = s0_view.getString()
		assert s0 == "foo"

		i1_view = ldroot.getView("tree1/y/i1")
		i1 = i1_view.getDataInt()
		assert i1 == 400

		f1_view = ldroot.getView("tree1/z/f1")
		f1 = f1_view.getDataFloat()
		assert f1 == 17.0

		i10_view = ldroot.getView("tree1/c/int10")
		i20_view = ldroot.getView("tree1/x/int20")

		v1_vals = i10_view.getDataArray()
		v2_vals = i20_view.getDataArray()

		for i in range(10):
			assert v1_vals[i] == i

		for i in range(20):
			assert v2_vals[i] == 2 * i

		# Delete the group so it is ready to use by the next protocol
		root.destroyGroup("tree1")
