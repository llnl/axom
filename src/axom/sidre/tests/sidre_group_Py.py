# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre

if pysidre.AXOM_USE_HDF5:
	NPROTOCOLS = 3
	PROTOCOLS = ["sidre_json", "sidre_hdf5", "json"]
else:
	NPROTOCOLS = 2
	PROTOCOLS = ["sidre_json", "json"]


# Verify getName()
def test_get_name():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("test")

	assert grp.getName() == "test"

	grp2 = root.getGroup("foo")
	assert grp2 == None


# Verify getPathName()
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


# Verify getParent()
def test_get_parent():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	parent = root.createGroup("parent")
	child = parent.createGroup("child")

	assert child.getParent() == parent


# Verify getDataStore()
def test_get_datastore():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("parent")

	assert grp.getDataStore() == ds

	const_ds = grp.getDataStore()
	assert const_ds == ds


# Verify getGroup()
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


# Verify getView()
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


# Verify getViewName() and getViewIndex()
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


# Verify getFirstValidGroupIndex() and getNextValidGroupIndex()
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


# Verify getFirstValidViewIndex() and getNextValidIndex()
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


# Verify getGroupName() and getGroupIndex()
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


# Verify createViewEmpty(), destoryView(), hasView()
def test_create_destroy_has_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	grp = root.createGroup("parent")
	view0 = grp.createView("view0")
	assert grp.getParent() == root
	assert not view0.hasBuffer()

	assert grp.hasView("view0")
	# Try creating view again, should be no-op
	view1 = grp.createView("view0")
	assert view1 == None

	grp.destroyView("view0")
	assert not grp.hasView("view0"), "grp.hasView(\"view0\")"

	# Try API call that specifies type and length
	view1 = grp.createViewAndAllocate("viewWithLength1", pysidre.TypeID.FLOAT_ID, 50)
	grp.destroyView("viewWithLength1")
	assert not grp.hasView("viewWithLength1"), "grp.hasView(\"viewWithLength1\")"

	view1 = grp.createViewAndAllocate("viewWithLengthBadLen", pysidre.TypeID.FLOAT64_ID, -1)
	assert view1 == None

	# Try API call that specifies data type in another way
	view1 = grp.createViewAndAllocate("viewWithLength2", pysidre.TypeID.FLOAT64_ID, 50)
	view2 = grp.createViewAndAllocate("viewWithLength2", pysidre.TypeID.FLOAT64_ID, 50)
	assert view2 == None
	# Destory this view using index
	grp.destroyViewAndData(grp.getFirstValidGroupIndex())


# Verify createGroup(), destroyGroup(), hasGroup()
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


def test_group_name_collisions():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	flds = root.createGroup("fields")
	view = flds.createView("a")

	assert flds.hasView("a")

	# Attempt to create duplicate group name
	badgrp = root.createGroup("fields")
	assert badgrp == None

	# Check error condition - attempt to create duplicate view name.
	view = flds.createView("a")
	assert view == None


# Restore this after copy_move is working. ATK-667 (?)
def test_view_copy_move():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	flds = root.createGroup("fields")

	i0_view = flds.createViewScalar("i0", 1)
	f0_view = flds.createViewScalar("f0", 100.0)
	d0_view = flds.createViewScalar("d0", 3000.0)

	assert flds.hasView("i0")
	assert flds.hasView("f0")
	assert flds.hasView("d0")

	# Test moving a view from flds to sub
	subgrp = flds.createGroup("sub")
	tmpview = subgrp.moveView(flds.getView("d0"))
	flds.print()
	assert not flds.hasView("d0")
	assert flds.hasGroup("sub")
	assert subgrp.hasView("d0")

	# TODO - functionality will change with numpy
	# Check data value
	# Conduit error as-is currently
	# assert tmpview.getData() == 3000.0

	# Test copying a view from flds to sub
	tmpview = subgrp.copyView(flds.getView("i0"))
	flds.print()

	assert flds.hasView("i0")
	assert subgrp.hasView("i0")

	# We expect teh actaul data points to be the same
	assert flds.getView("i0").getBuffer() == flds.getGroup("sub").getView("i0").getBuffer()


def test_groups_move_copy():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	flds = root.createGroup("fields")

	ga = flds.createGroup("a")
	gb = flds.createGroup("b")
	gc = flds.createGroup("c")

	i0_view = ga.createViewScalar("i0", 1)
	f0_view = gb.createViewScalar("f0", 100.0)
	d0_view = gc.createViewScalar("d0", 300.0)

	# Check that all sub groups exist
	assert flds.hasGroup("a")
	assert flds.hasGroup("b")
	assert flds.hasGroup("c")

	# Move "b" to a child of "sub"
	subgrp = flds.createGroup("sub")
	tmpgrp = subgrp.moveGroup(gb)

	flds.print()

	assert flds.hasGroup("a")
	assert flds.hasGroup("sub")
	assert flds.hasGroup("c")

	assert tmpgrp == gb


def test_create_destroy_view_and_data():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	grp = root.createGroup("grp")

	view_name1 = "viewBuffer1"
	view_name2 = "viewBuffer2"

	view1 = grp.createViewAndAllocate(view_name1, pysidre.TypeID.INT_ID, 1)
	view2 = grp.createViewAndAllocate(view_name2, pysidre.TypeID.INT_ID, 1)

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
	view1 = grp.createViewAndAllocate(view_name1, pysidre.TypeID.INT_ID, 10)
	assert grp.hasView(view_name1)
	assert grp.getView(view_name1) == view1

	# TODO getData, need numpy array implementation
	# v1_vals = view1.getData()


def test_create_view_of_buffer_with_datatype():
	# TODO getData, need numpy array implementation
	pass


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

		# TODO functionality will change with numpy getData()
		view = root2.getView("i0")
		i0 = view.getData()
		assert i0 == 1

		# Wrong 0 == 1.0
		# view = root2.getView("f0")
		# f0 = view.getData()
		# assert f0 == 1.0

		# Wrong 0 == 1.0
		# view = root2.getView("d0")
		# d0 = view.getData()
		# assert d0 == 10.0

		view = root2.getView("s0")
		s0 = view.getString()
		assert s0 == "I am a string"


def test_save_restore_external_data():
	# TODO getData, need numpy array implementation for external data
	pass


def test_save_restore_other():
	# TODO need numpy array for createViewWithShape
	# file_path_base = "py_sidre_empty_other_"
	# ds1 = pysidre.DataStore()
	# root1 = ds1.getRoot()
	pass


def test_rename_group():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	child1 = root.createGroup("g_a")
	child2 = root.createGroup("g_b")
	child3 = root.createGroup("g_c")

	assert child1.rename("g_r")
	assert child1.getName() == "g_r"
	assert root.hasGroup("g_r")
	assert not root.hasGroup("g_a")

	assert not child2.rename("fields/g_s")
	assert child2.getName() == "g_b"

	assert not child3.rename("g_b")
	assert child3.getName() == "g_c"


# Fortran comment - TODO - redo these, the C++ tests were heavily rewritten
#  def test_call save_restore_simple
#  def test_call save_restore_complex
# Complex requires numpy getArray() implementation
# def test_save_load_preserve_contents
# TODO requires numpy getArray() implementation
