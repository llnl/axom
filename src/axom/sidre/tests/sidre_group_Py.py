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

	other_ds = grp.getDataStore()
	assert other_ds == ds


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

	assert not flds.hasView("d0")
	assert flds.hasGroup("sub")
	assert subgrp.hasView("d0")
	assert tmpview.getDataFloat() == 3000.0

	# Test copying a view from flds to sub
	tmpview = subgrp.copyView(flds.getView("i0"))

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

	assert child1.rename("g_r")
	assert child1.getName() == "g_r"
	assert root.hasGroup("g_r")
	assert not root.hasGroup("g_a")

	assert not child2.rename("fields/g_s")
	assert child2.getName() == "g_b"

	assert not child3.rename("g_b")
	assert child3.getName() == "g_c"


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
