# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre


EMPTYVIEW = 1
BUFFERVIEW = 2
EXTERNALVIEW = 3
SCALARVIEW = 4
STRINGVIEW = 5
NOTYPE = 6

# Helper function to get state
def get_state(view):
	if view.isEmpty():
		return EMPTYVIEW
	elif view.hasBuffer():
		return BUFFERVIEW
	elif view.isExternal():
		return EXTERNALVIEW
	elif view.isScalar():
		return SCALARVIEW
	elif view.isString():
		return STRINGVIEW
	else:
		return NOTYPE


# Helper function to check values
def check_view_values(view, state, is_described, is_allocated, is_applied, length):
	dims = [0,0]

	name = view.getName()
	assert get_state(view) == state
	assert view.isDescribed() == is_described, f"{name} is described"
	assert view.isAllocated() == is_allocated, f"{name} is allocated"
	assert view.isApplied() == is_applied, f"{name} is applied"
	assert view.getNumElements() == length, f"{name} getNumElements"

	if view.isDescribed():
		assert view.getNumDimensions() == 1, f"{name} getNumDimensions"
	
		# TODO Implement numpy array for getShape()
		# assert view.getShape() == 1, f"{name} getShape"
		# assert dims[0] == length, f"{name} dims[0]"


def test_create_views():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	dv_0 = root.createViewAndAllocate("field0", pysidre.TypeID.INT_ID, 1)
	dv_1 = root.createViewAndAllocate("field1", pysidre.TypeID.INT_ID, 1)

	db_0 = dv_0.getBuffer()
	db_1 = dv_1.getBuffer()

	assert db_0.getIndex() == 0
	assert db_1.getIndex() == 1


def test_get_path_name():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	v1 = root.createView("test/a/b/v1")
	v2 = root.createView("test/v2")
	v3 = root.createView("v3")

	assert v1.getName() == "v1"
	assert v1.getPath() == "test/a/b"
	assert v1.getPathName() == "test/a/b/v1"

	assert v2.getName() == "v2"
	assert v2.getPath() == "test"
	assert v2.getPathName() == "test/v2"

	assert v3.getName() == "v3"
	assert v3.getPath() == ""
	assert v3.getPathName() == "v3"


def test_create_view_from_path():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	baz = root.createView("foo/bar/baz")
	assert root.hasGroup("foo")
	assert root.getGroup("foo").hasGroup("bar")

	bar = root.getGroup("foo").getGroup("bar")
	assert bar.hasView("baz")
	assert bar.getView("baz") == baz


def test_scalar_view():
	# Inner helper function
	def check_scalar_values(view, state, is_described, is_allocated, is_applied, typeID, length):
		dims = [0,0]

		name = view.getName()
		assert get_state(view) == state
		assert view.isDescribed() == is_described, f"{name} is described"
		assert view.isAllocated() == is_allocated, f"{name} is allocated"
		assert view.isApplied() == is_applied, f"{name} is applied"

		assert view.getTypeID() == typeID, f"{name} getTypeID"
		assert view.getNumElements() == length, f"{name} getNumElements"
		assert view.getNumDimensions() == 1, f"{name} getNumDimensions"
		
		# TODO Implement numpy array for getShape()
		# assert view.getShape() == 1, f"{name} getShape"
		# assert dims[0] == length, f"{name} dims[0]"

	ds = pysidre.DataStore()
	root = ds.getRoot()

	i1 = 1
	i0view = root.createView("i0")
	i0view.setScalar(i1)
	check_scalar_values(i0view, SCALARVIEW, True, True, True, pysidre.TypeID.INT32_ID, 1)
	i2 = i0view.getData()
	assert i1 == i2

	i1 = 2
	i1view = root.createViewScalar("i1", i1)
	check_scalar_values(i1view, SCALARVIEW, True, True, True, pysidre.TypeID.INT32_ID, 1)
	i2 = i1view.getData()
	assert i1 == i2

	s1 = "i am a string"
	s0view = root.createView("s0")
	s0view.setString(s1)
	check_scalar_values(s0view, STRINGVIEW, True, True, True, pysidre.TypeID.CHAR8_STR_ID, len(s1) + 1)
	s2 = s0view.getString()
	assert s1 == s2

	s1 = "i too am a string"
	s1view = root.createViewString("s1", s1)
	check_scalar_values(s1view, STRINGVIEW, True, True, True, pysidre.TypeID.CHAR8_STR_ID, len(s1) + 1)
	s2 = s1view.getString()
	assert s1 == s2

	# Test group access to scalars
	# TODO - (shroud specific) Group access to View data?



#def test_dealloc():

#def test_alloc_zero_items():

#def test_alloc_and_dealloc_multiview():

def test_int_buffer_from_view():
	elem_count = 10
	ds = pysidre.DataStore()
	root = ds.getRoot()

	dv = root.createViewAndAllocate("u0", pysidre.TypeID.INT32_ID, elem_count)
	data = dv.getData()

	# TODO Implement numpy to getData()



def test_int_buffer_from_view_conduit_value():
	# TODO Implement numpy to getData()
	pass


def test_int_array_multi_view():
	# TODO Implement numpy to getData() / getVoidPtr()
	pass


def test_init_int_array_multi_view():
	# TODO Implement numpy to getData() / getVoidPtr()
	pass


def test_int_array_depth_view():
	# TODO Implement numpy to getData() / getVoidPtr()
	pass


def test_int_array_view_attach_buffer():
	# TODO Implement numpy to getData() / getVoidPtr()
	pass


def test_int_array_offset_stride():
	# TODO Implement numpy to getData() / getVoidPtr()
	pass


def test_int_array_multi_view_resize():
	# TODO Implement numpy to getData() / getVoidPtr()
	pass


def test_int_array_realloc():
	# TODO Implement numpy to getData() / getVoidPtr()
	pass


def test_simple_opaque():
	# TODO Implement numpy to getData() / getVoidPtr() / createViewExternal(ptr)
	pass


def test_clear_view():
	BLEN = 10
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# Create an empty view
	view = root.createView("v_empty")
	check_view_values(view, EMPTYVIEW, False, False, False, 0)
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)

	# Describe an empty view
	view = root.createView("v_described", pysidre.TypeID.INT32_ID, BLEN)
	check_view_values(view, EMPTYVIEW, True, False, False, BLEN)
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)

	# Scalar view
	view = root.createViewScalar("v_scalar", 1)
	check_view_values(view, SCALARVIEW, True, True, True, 1)
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)

	# String view
	view = root.createViewString("v_string", "string-test")
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)

	# Allocated view, Buffer will be released
	nbuf = ds.getNumBuffers()
	view = root.createViewAndAllocate("v_allocated", pysidre.TypeID.INT32_ID, BLEN)
	check_view_values(view, BUFFERVIEW, True, True, True, BLEN)
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)
	assert ds.getNumBuffers() == nbuf

	# Undescribed buffer
	nbuf = ds.getNumBuffers()
	dbuff = ds.createBuffer()
	view = root.createView("v_undescribed_buffer", dbuff)
	check_view_values(view, BUFFERVIEW, False, False, False, 0)
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)
	assert ds.getNumBuffers() == nbuf

	# Explicit buffer attached to two views
	dbuff = ds.createBuffer()
	dbuff.allocate(pysidre.TypeID.INT32_ID, BLEN)
	nbuf = ds.getNumBuffers()
	assert dbuff.getNumViews() == 0

	vother = root.createView("v_other", pysidre.TypeID.INT32_ID, BLEN)
	view = root.createView("v_buffer", pysidre.TypeID.INT32_ID, BLEN)
	vother.attachBuffer(dbuff)
	assert dbuff.getNumViews() == 1
	view.attachBuffer(dbuff)
	assert dbuff.getNumViews() == 2

	check_view_values(view, BUFFERVIEW, True, True, True, BLEN)
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)

	assert ds.getNumBuffers() == nbuf
	assert dbuff.getNumViews() == 1

	# External View
	# TODO Implement numpy to createView(ptr)

	# Opaque view
	# TODO Implement numpy to createView(ptr)
