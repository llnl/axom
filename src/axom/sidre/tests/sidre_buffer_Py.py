# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre

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

	dbuff.allocate(pysidre.TypeID.INT_ID, 10)
	dbuff.allocate()

    # Type changed to INT32_ID? Why?
	assert dbuff.getTypeID() == pysidre.TypeID.INT32_ID
	assert dbuff.getNumElements() == 10
	assert dbuff.getTotalBytes() == NUM_BYTES_INT_32 * 10

	# TODO Requires nanobind's numpy support to get access to underlying data instead of CapsuleType
	data_ptr = dbuff.getVoidPtr()
	
	# for i in range(10):
		# data_ptr[i] = i * i

	dbuff.print()
	ds.print()


def test_init_buffer_for_int_array():
	ds = pysidre.DataStore()
	dbuff = ds.createBuffer()

	dbuff.allocate(pysidre.TypeID.INT_ID, 10)

    # Type changed to INT32_ID? Why?
	assert dbuff.getTypeID() == pysidre.TypeID.INT32_ID
	assert dbuff.getNumElements() == 10
	assert dbuff.getTotalBytes() == NUM_BYTES_INT_32 * 10

	# TODO Requires nanobind's numpy support to get access to underlying data instead of CapsuleType
	data_ptr = dbuff.getVoidPtr()
	
	# for i in range(10):
		# data_ptr[i] = i * i

	dbuff.print()
	ds.print()


def test_realloc_buffer():
	ds = pysidre.DataStore()
	dbuff = ds.createBuffer()

	dbuff.allocate(pysidre.TypeID.INT_ID, 5)

    # Type changed to INT32_ID? Why?
	assert dbuff.getTypeID() == pysidre.TypeID.INT32_ID
	assert dbuff.getNumElements() == 5
	assert dbuff.getTotalBytes() == NUM_BYTES_INT_32 * 5

	# TODO Requires nanobind's numpy support to get access to underlying data instead of CapsuleType
	data_ptr = dbuff.getVoidPtr()

	# for i in range(5):
		# data_ptr[i] = 5

	# for i in range(5):
		# assert data_ptr[i] == 5

	dbuff.reallocate(10)

    # Type changed to INT32_ID? Why?
	assert dbuff.getTypeID() == pysidre.TypeID.INT32_ID
	assert dbuff.getNumElements() == 10
	assert dbuff.getTotalBytes() == NUM_BYTES_INT_32 * 10

	# TODO Requires nanobind's numpy support to get access to underlying data instead of CapsuleType
	data_ptr = dbuff.getVoidPtr()

	# for i in range(5,10):
		# data_ptr[i] = 10

	# for i in range(0,10):
		# value = 5
		# if i > 4:
		# 	value = 10
		# assert data_ptr[i] == value