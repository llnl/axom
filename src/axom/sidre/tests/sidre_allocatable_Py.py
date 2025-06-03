# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import numpy as np

# Allocate array via python
# Register with datastore then 
# Query metadata using datastore API.
def test_external_allocatable_int():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	iarray = np.array(range(1,11))

	view = root.createView("iarray", iarray)
	view.apply(pysidre.TypeID.INT64_ID, 10)

	assert view.isExternal() == True
	assert view.getTypeID() == pysidre.TypeID.INT64_ID
	assert view.getNumElements() == len(iarray)
	assert view.getNumDimensions() == 1

	extents = np.zeros(7)
	rank,extents = view.getShape(7, extents)
	assert rank == 1
	assert extents[0] == np.size(iarray)

	ipointer = view.getDataArray()
	assert np.array_equal(ipointer, iarray)


def test_external_allocatable_int_3d():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# create 3D numpy array
	iarray = np.empty((2, 3, 4), dtype=int)

	for i in range(2):
	    for j in range(3):
	        for k in range(4):
	            iarray[i, j, k] = (i+1)*100 + (j+1)*10 + (k+1)
	view = root.createView("iarray", iarray)
	view.apply(pysidre.TypeID.INT64_ID, 3, np.array([2,3,4]))

	assert view.isExternal() == True
	assert view.getTypeID() == pysidre.TypeID.INT64_ID
	assert view.getNumElements() == np.size(iarray)
	assert view.getNumDimensions() == 3

	extents = np.zeros(7)
	rank,extents = view.getShape(7, extents)
	assert rank == 3
	assert extents[0] == iarray.shape[0]
	assert extents[1] == iarray.shape[1]
	assert extents[2] == iarray.shape[2]

	ipointer = view.getDataArray()
	assert np.array_equal(ipointer, iarray)


# register a static (non-allocatable) array with the datastore as external view
def test_external_static_int():
	# may not be relevant for Python?
	# TODO
	pass


# check other types
def test_external_allocatable_double():
	# TODO
	pass


# Datastore owns a multi-dimension array.
def test_datastore_int_3d():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	extents = [2,3,4]

	view = root.createViewWithShapeAndAllocate("iarray", pysidre.TypeID.INT32_ID, 3, extents)

	# TODO return a python numpy array
	# data = view.getData_Int()
	# print(data)
	# print(type(data))

	# typeID = view.getTypeID()
	# assert typeID == pysidre.TypeID.INT32_ID

	# num_elements = view.getNumElements()
	# assert num_elements == len(data)

	# rank = view.getNumDimensions()
	# assert rank == 3

	# rank = view.getShape(7, extents)
	# assert rank == 3

