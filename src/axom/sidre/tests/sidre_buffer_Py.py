# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import numpy as np

NUM_BYTES_INT_32 = 4

def test_create_buffers():
	ds = pysidre.DataStore()
	assert ds.getNumBuffers() == 0

	dbuff_0 = ds.createBuffer()
	assert ds.getNumBuffers() == 1
	assert dbuff_0.getIndex() == 0

	dbuff_1 = ds.createBuffer()
	assert ds.getNumBuffers() == 2
	assert dbuff_1.getIndex() == 1

	# Destroy by index
	ds.destroyBuffer(0)
	assert ds.getNumBuffers() == 1

	dbuff_0b = ds.createBuffer()
	assert ds.getNumBuffers() == 2
	assert dbuff_0b.getIndex() == 0

	ds.print()


def test_alloc_buffer_for_int_array():
	ds = pysidre.DataStore()
	dbuff = ds.createBuffer()
	elem_count = 10

	dbuff.allocate(pysidre.TypeID.INT32_ID, elem_count)

	# Should be a warning and no-op, buffer is already allocated, we don't want
	# to re-allocate and leak memory.
	dbuff.allocate()

	assert dbuff.getTypeID() == pysidre.TypeID.INT32_ID
	assert dbuff.getNumElements() == elem_count
	assert dbuff.getTotalBytes() == NUM_BYTES_INT_32 * elem_count

	data = dbuff.getDataArray()

	assert type(data[0]) == np.int32
	
	for i in range(elem_count):
		data[i] = i * i

	for i in range(elem_count):
		assert data[i] == i * i

	dbuff.print()
	ds.print()


def test_init_buffer_for_int_array():
	elem_count = 10

	ds = pysidre.DataStore()
	dbuff = ds.createBuffer()

	dbuff.allocate(pysidre.TypeID.INT32_ID, elem_count)

	assert dbuff.getTypeID() == pysidre.TypeID.INT32_ID
	assert dbuff.getNumElements() == elem_count
	assert dbuff.getTotalBytes() == NUM_BYTES_INT_32 * elem_count

	data = dbuff.getDataArray()

	assert type(data[0]) == np.int32

	for i in range(elem_count):
		data[i] = i * i

	for i in range(elem_count):
		assert data[i] == i * i

	dbuff.print()
	ds.print()


def test_realloc_buffer():
	orig_elem_count = 5
	mod_elem_count = 10

	ds = pysidre.DataStore()
	dbuff = ds.createBuffer()

	dbuff.allocate(pysidre.TypeID.INT32_ID, orig_elem_count)

	assert dbuff.getTypeID() == pysidre.TypeID.INT32_ID
	assert dbuff.getNumElements() == orig_elem_count
	assert dbuff.getTotalBytes() == NUM_BYTES_INT_32 * orig_elem_count

	data = dbuff.getDataArray()

	for i in range(orig_elem_count):
		data[i] = orig_elem_count

	for i in range(orig_elem_count):
		assert data[i] == orig_elem_count

	dbuff.reallocate(mod_elem_count)

	assert dbuff.getTypeID() == pysidre.TypeID.INT32_ID
	assert dbuff.getNumElements() == mod_elem_count
	assert dbuff.getTotalBytes() == NUM_BYTES_INT_32 * mod_elem_count

	data = dbuff.getDataArray()

	assert type(data[0]) == np.int32

	for i in range(orig_elem_count,mod_elem_count):
		data[i] = mod_elem_count

	for i in range(0,mod_elem_count):
		value = orig_elem_count
		if i > 4:
			value = mod_elem_count
		assert data[i] == value