# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import numpy as np

NUM_BYTES_INT_32 = 4
NUM_BYTES_DOUBLE = 8

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
	dims = np.array([0,0])

	name = view.getName()
	assert get_state(view) == state
	assert view.isDescribed() == is_described, f"{name} is described"
	assert view.isAllocated() == is_allocated, f"{name} is allocated"
	assert view.isApplied() == is_applied, f"{name} is applied"
	assert view.getNumElements() == length, f"{name} getNumElements"

	if view.isDescribed():
		assert view.getNumDimensions() == 1, f"{name} getNumDimensions"
	
		ndims, dims = view.getShape(1, dims)
		assert ndims == 1, f"{name} getShape"
		assert dims[0] == length, f"{name} dims[0]"


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
		dims = np.array([0,0])

		name = view.getName()
		assert get_state(view) == state
		assert view.isDescribed() == is_described, f"{name} is described"
		assert view.isAllocated() == is_allocated, f"{name} is allocated"
		assert view.isApplied() == is_applied, f"{name} is applied"

		assert view.getTypeID() == typeID, f"{name} getTypeID"
		assert view.getNumElements() == length, f"{name} getNumElements"
		assert view.getNumDimensions() == 1, f"{name} getNumDimensions"
		
		ndims, dims = view.getShape(1, dims)
		assert ndims == 1, f"{name} getShape"
		assert dims[0] == length, f"{name} dims[0]"

	ds = pysidre.DataStore()
	root = ds.getRoot()

	i1 = 1
	i0view = root.createView("i0")
	i0view.setScalar(i1)
	check_scalar_values(i0view, SCALARVIEW, True, True, True, pysidre.TypeID.INT32_ID, 1)
	i2 = i0view.getDataInt()
	assert i1 == i2

	i1 = 2
	i1view = root.createViewScalar("i1", i1)
	check_scalar_values(i1view, SCALARVIEW, True, True, True, pysidre.TypeID.INT32_ID, 1)
	i2 = i1view.getDataInt()
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


# Test functions that are commented out in Fortran's "main" method, but unimplemented.
#def test_dealloc():

#def test_alloc_zero_items():

#def test_alloc_and_dealloc_multiview():

def test_int_buffer_from_view():
	elem_count = 10
	ds = pysidre.DataStore()
	root = ds.getRoot()

	dv = root.createViewAndAllocate("u0", pysidre.TypeID.INT32_ID, elem_count)
	data = dv.getDataArray()

	for i in range(elem_count):
		data[i] = i * i

	dv.print()

	assert dv.getNumElements() == elem_count
	assert dv.getTotalBytes() == NUM_BYTES_INT_32 * elem_count


def test_int_array_multi_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	dbuff = ds.createBuffer(pysidre.TypeID.INT32_ID, 10)

	dbuff.allocate()
	data = dbuff.getDataArray()

	for i in range(10):
		data[i] = i

	dbuff.print()

	dv_e = root.createView("even", dbuff)
	dv_o = root.createView("odd", dbuff)
	assert dbuff.getNumViews() == 2, f"{dbuff.getNumViews()} == 2"

	# NOTE - The data/data ptr here is different from C++ implementation.
	#        Buffer representations are the same.
	# "even"
	# C++ : [0,1,2,3,4,5,6,7,8,9]
	# python: [0,2,4,6,8]
	# "odd"
	# C++ : [1,2,3,4,5,6,7,8,9, <garbage val>]
	# python: [1,3,5,7,9]

	dv_e.apply(5,0,2)
	dv_e.print()
	data_e = dv_e.getDataArray()

	dv_o.apply(5,1,2)
	dv_o.print()
	data_o = dv_o.getDataArray()

	for i in range(5):
		assert data_e[i] % 2 == 0
		assert data_o[i] % 2 == 1


def test_init_int_array_multi_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	dbuff = ds.createBuffer()
	dbuff.allocate(pysidre.TypeID.INT32_ID, 10)

	data = dbuff.getDataArray()
	for i in range(10):
		data[i] = i
	dbuff.print()

	dv_e = root.createView("even", dbuff)
	dv_o = root.createView("odd", dbuff)
	assert dbuff.getNumViews() == 2, f"{dbuff.getNumViews()} == 2"

	# NOTE - The data/data ptr here is different from C++ implementation.
	#        Buffer representations are the same.
	# "even"
	# C++ : [0,1,2,3,4,5,6,7,8,9]
	# python: [0,2,4,6,8]
	# "odd"
	# C++ : [1,2,3,4,5,6,7,8,9, <garbage val>]
	# python: [1,3,5,7,9]

	dv_e.apply(5,0,2)
	dv_e.print()
	data_e = dv_e.getDataArray()

	dv_o.apply(5,1,2)
	dv_o.print()
	data_o = dv_o.getDataArray()

	for i in range(5):
		assert data_e[i] % 2 == 0
		assert data_o[i] % 2 == 1



def test_int_array_depth_view():
	ds = pysidre.DataStore()
	depth_nelems = 10
	total_nelems = 4 * depth_nelems

	dbuff = ds.createBuffer(pysidre.TypeID.INT32_ID, total_nelems)

	# Get access to our root data Group
	root = ds.getRoot()

	# Allocate buffer to hold data for 4 "depth" views
	dbuff.allocate()

	data = dbuff.getDataArray()

	for i in range(total_nelems):
		data[i] = i / depth_nelems

	dbuff.print()
	assert dbuff.getNumElements() == 4 * depth_nelems

	# Create 4 "depth" views and apply offsets into buffer
	view0 = root.createView("depth_view_0", dbuff)
	view0.apply(depth_nelems, 0 * depth_nelems)

	view1 = root.createView("depth_view_1", dbuff)
	view1.apply(depth_nelems, 1 * depth_nelems)

	view2 = root.createView("depth_view_2", dbuff)
	view2.apply(depth_nelems, 2 * depth_nelems)

	view3 = root.createView("depth_view_3", dbuff)
	view3.apply(depth_nelems, 3 * depth_nelems)

	assert dbuff.getNumViews() == 4

	view0.print()
	view1.print()
	view2.print()
	view3.print()

	# Check values in depth views
	data = view0.getDataArray()
	assert np.all(data == 0)

	data = view1.getDataArray()
	assert np.all(data == 1)

	data = view2.getDataArray()
	assert np.all(data == 2)

	data = view3.getDataArray()
	assert np.all(data == 3)


def test_int_array_view_attach_buffer():
	# Create our main data store
	ds = pysidre.DataStore()

	# Get access to our root data Group
	root = ds.getRoot()

	field_nelems = 10

	# Create 2 "field" views with type and # elems
	elem_count = 0
	field0 = root.createView("field0", pysidre.TypeID.INT32_ID, field_nelems)
	elem_count = elem_count + field0.getNumElements()
	print(f"elem_count field0 {elem_count}")
	field1 = root.createView("field1", pysidre.TypeID.INT32_ID, field_nelems)
	elem_count = elem_count + field1.getNumElements()
	print(f"elem_count field1 {elem_count}")

	assert elem_count == 2 * field_nelems

	# Create buffer to hold data for all fields and allocate
	dbuff = ds.createBuffer(pysidre.TypeID.INT32_ID, elem_count)
	dbuff.allocate()

	assert dbuff.getNumElements() == elem_count

	# Initialize buffer data for testing below
	data = dbuff.getDataArray()
	for i in range(elem_count):
		data[i] = i / field_nelems

	dbuff.print()

	# Attach field views to buffer and apply offsets into buffer
	offset0 = 0 * field_nelems
	field0.attachBuffer(dbuff)
	field0.apply(field_nelems, offset0)

	offset1 = 1 * field_nelems
	field1.attachBuffer(dbuff)
	field1.apply(field_nelems, offset1)

	assert dbuff.getNumViews() == 2

	# Print field views...
	field0.print()
	field1.print()

	# Check values in field views
	data = field0.getDataArray()
	assert np.size(data) == field_nelems
	assert np.all(data == 0)
	assert field0.getOffset() == offset0

	data = field1.getDataArray()
	assert np.size(data) == field_nelems
	assert np.all(data == 1)
	assert field1.getOffset() == offset1



def test_int_array_offset_stride():
	# create our main data store
	ds = pysidre.DataStore()

	# get access to our root data Group
	root = ds.getRoot()

	field_nelems = 20
	field0 = root.createViewAndAllocate("field0", pysidre.TypeID.DOUBLE_ID, field_nelems)
	assert field0.getNumElements() == field_nelems
	assert field0.getBytesPerElement() == NUM_BYTES_DOUBLE
	assert field0.getTotalBytes() == NUM_BYTES_DOUBLE * field_nelems
	assert field0.getOffset() == 0
	assert field0.getStride() == 1

	dbuff = field0.getBuffer()

	# Initialize buffer data for testing below
	data = dbuff.getDataArray()

	for i in range(field_nelems):
		data[i] = (i+1) * 1.001

	dbuff.print()

	# Create two more views into field0's buffer and test stride and offset
	v1_nelems = 3
	v1_stride = 3
	v1_offset = 2
	view1 = root.createView("offset_stride1", dbuff)
	view1.apply(v1_nelems, v1_offset, v1_stride)
	data1 = view1.getDataArray()

	v2_nelems = 3
	v2_stride = 3
	v2_offset = 3
	view2 = root.createView("offset_stride2", dbuff)
	view2.apply(v2_nelems, v2_offset, v2_stride)
	data2 = view2.getDataArray()

	v3_nelems = 5
	v3_stride = 1
	v3_offset = 12
	view3 = root.createView("offset_stride3", dbuff)
	view3.apply(v3_nelems, v3_offset, v3_stride)
	data3 = view3.getDataArray()

	assert view1.getNumElements() == v1_nelems
	assert view1.getBytesPerElement() == NUM_BYTES_DOUBLE
	assert view1.getOffset() == v1_offset
	assert view1.getStride() == v1_stride
	assert view1.getTotalBytes() == NUM_BYTES_DOUBLE * (1 + (v1_stride * (v1_nelems - 1)))
	assert data1[0] == data[v1_offset]

	assert view2.getNumElements() == v2_nelems
	assert view2.getBytesPerElement() == NUM_BYTES_DOUBLE
	assert view2.getOffset() == v2_offset
	assert view2.getStride() == v2_stride
	assert view2.getTotalBytes() == NUM_BYTES_DOUBLE * (1 + (v2_stride * (v2_nelems - 1)))
	assert data2[0] == data[v2_offset]

	assert view3.getNumElements() == v3_nelems
	assert view3.getBytesPerElement() == NUM_BYTES_DOUBLE
	assert view3.getOffset() == v3_offset
	assert view3.getStride() == v3_stride
	assert view3.getTotalBytes() == NUM_BYTES_DOUBLE * (1 + (v3_stride * (v3_nelems - 1)))
	assert data3[0] == data[v3_offset]

	# Test stride and offset against other types of views
	other = root.createGroup("other")
	view1 = other.createView("key_empty")
	assert view1.getOffset() == 0
	assert view1.getStride() == 1
	assert view1.getNumElements() == 0
	assert view1.getBytesPerElement() == 0

	# Opaque - not described
	view1 = other.createView("key_opaque", data)
	assert view1.getOffset() == 0
	assert view1.getStride() == 1
	assert view1.getNumElements() == 0
	assert view1.getBytesPerElement() == 0
	assert view1.getTotalBytes() == 0

	view1 = other.createViewString("key_str", "val_str")
	assert view1.getOffset() == 0
	assert view1.getStride() == 1
	assert view1.getBytesPerElement() == 1

	view1 = other.createViewScalar("key_int", 5)
	int_data = view1.getDataInt()
	assert int_data == 5
	assert view1.getOffset() == 0
	assert view1.getStride() == 1
	assert view1.getNumElements() == 1
	assert view1.getBytesPerElement() == NUM_BYTES_INT_32
	assert view1.getTotalBytes() == NUM_BYTES_INT_32


def test_int_array_multi_view_resize():
	# This example creates a 4 * 10 buffer of ints,
	# and 4 views that point the 4 sections of 10 ints

	# We then create a new buffer to support 4*12 ints
	# and 4 views that point into them

	# after this we use the old buffers to copy the values
	# into the new views

	# Create our main data store
	ds = pysidre.DataStore()

	# Get access to our root data Group
	root = ds.getRoot()

	# Create a group to hold the "old" or data we want to copy
	r_old = root.createGroup("r_old")

	# Create a view to hold the base buffer and allocate
	# we will create 4 sub views of this array
	base_old = r_old.createViewAndAllocate("base_data", pysidre.TypeID.INT32_ID, 40)


	# Init the buff with values that align with the 4 subsections
	data = base_old.getDataArray()
	data[0:10] = 1
	data[10:20] = 2
	data[20:30] = 3
	data[30:40] = 4

	# Setup our 4 views
	buff_old = base_old.getBuffer()
	r0_old = r_old.createView("r0", buff_old)
	r1_old = r_old.createView("r1", buff_old)
	r2_old = r_old.createView("r2", buff_old)
	r3_old = r_old.createView("r3", buff_old)

	# Each view is offset by 10
	offset = 0
	r0_old.apply(10, offset)
	offset += 10
	r1_old.apply(10, offset)
	offset += 10
	r2_old.apply(10, offset)
	offset += 10
	r3_old.apply(10, offset)

	# Check that our views actually point to the expected data
	r0_old_data = r0_old.getDataArray()
	r1_old_data = r1_old.getDataArray()
	r2_old_data = r2_old.getDataArray()
	r3_old_data = r3_old.getDataArray()

	for i in range(10):
		assert r0_old_data[i] == 1
		assert r1_old_data[i] == 2
		assert r2_old_data[i] == 3
		assert r3_old_data[i] == 4

	# Create a group to hold the "new" or data we want to copy into
	r_new = root.createGroup("r_new")

	# Create a view to hold the base buffer
	base_new = r_new.createView("base_data")
	base_new.allocate(pysidre.TypeID.INT32_ID, 48)
	base_new_data = base_new.getDataArray()
	for i in range(48):
		base_new_data[i] = 0

	buff_new = base_new.getBuffer()
	buff_new.print()

	# Create the 4 sub views of this array
	r0_new = r_new.createView("r0", buff_new)
	r1_new = r_new.createView("r1", buff_new)
	r2_new = r_new.createView("r2", buff_new)
	r3_new = r_new.createView("r3", buff_new)

	# Apply views to r0, r1, r2, r3
	# Each view is offset by 12
	offset = 0
	r0_new.apply(12, offset)
	offset += 12
	r1_new.apply(12, offset)
	offset += 12
	r2_new.apply(12, offset)
	offset += 12
	r3_new.apply(12, offset)

	"""
	# Note - this requires getNode() (Conduit Node)
	! update r2 as an example first
	call buff_new%print()
	call r2_new%print()

	! copy the subset of value
	r2_new->getNode().update(r2_old->getNode())
	call r2_new%print()
	call buff_new%print()


	! check pointer values
	int * r2_new_ptr = (int *) SIDRE_view_get_data_pointer(r2_new)

	for(int i=0  i<10  i++)
	{
	EXPECT_EQ(r2_new_ptr[i], 3)
	}

	for(int i=10  i<12  i++)
	{
	EXPECT_EQ(r2_new_ptr[i], 0)     ! assumes zero-ed alloc
	}


	! update the other views
	r0_new->getNode().update(r0_old->getNode())
	r1_new->getNode().update(r1_old->getNode())
	r3_new->getNode().update(r3_old->getNode())

	call buff_new%print()
	"""


def test_int_array_realloc():
	# Create our main data store
	ds = pysidre.DataStore()

	# Get access to our root data Group
	root = ds.getRoot()

	a1 = root.createViewAndAllocate("a1", pysidre.TypeID.DOUBLE_ID, 5)
	a2 = root.createViewAndAllocate("a2", pysidre.TypeID.DOUBLE_ID, 5)

	a1_data = a1.getDataArray()
	a2_data = a2.getDataArray()

	assert np.size(a1_data) == 5
	assert np.size(a2_data) == 5

	a1_data[0:5] = 5.0
	a2_data[0:5] = -5.0

	a1.reallocate(10)
	a2.reallocate(15)

	a1_data = a1.getDataArray()
	a2_data = a2.getDataArray()

	assert np.size(a1_data) == 10
	assert np.size(a2_data) == 15

	assert np.all(a1_data[0:5] == 5.0)
	assert np.all(a2_data[0:5] == -5.0)

	a1_data[5:10] = 10.0
	a2_data[5:10] = -10.0
	a2_data[10:15] = -15.0

	assert np.all(a1_data[5:10] == 10.0)
	assert np.all(a2_data[5:10] == -10.0)
	assert np.all(a2_data[10:15] == -15.0)

	ds.print()


def test_simple_opaque():
	# Create our main data store
	ds = pysidre.DataStore()

	# Get access to our root data Group
	root = ds.getRoot()

	src_data = np.array([42])
	opq_view = root.createView("my_opaque", src_data)

	# External data is held in the view, not buffer
	assert ds.getNumBuffers() == 0

	assert opq_view.isExternal() == True
	assert opq_view.isApplied() == False
	assert opq_view.isOpaque() == True
	assert opq_view.getTypeID() == pysidre.TypeID.NO_TYPE_ID

	# Apply type to get data
	opq_view.apply(pysidre.TypeID.INT32_ID, 1)
	opq_data = opq_view.getDataArray()
	assert opq_data[0] == 42

	ds.print()


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
	ext_data = np.array(BLEN)
	view = root.createView("v_external", pysidre.TypeID.INT32_ID, BLEN, ext_data)
	check_view_values(view, EXTERNALVIEW, True, True, True, BLEN)
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)


	# Opaque view
	ext_data = np.array(BLEN)
	view = root.createView("v_opaque", ext_data)
	check_view_values(view, EXTERNALVIEW, False, True, False, 0)
	view.clear()
	check_view_values(view, EMPTYVIEW, False, False, False, 0)


