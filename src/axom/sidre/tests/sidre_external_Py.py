# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import numpy as np

def test_create_external_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	length = 11

	idata = np.array(range(length))
	print(f"PYTHON SIDE: idata type is {type(idata[0])}")
	print(f"PYTHON SIDE: idata is {idata}")

	ddata = np.array([x * 2.0 for x in range(length)])
	print(f"PYTHON SIDE: ddata type is {type(ddata[0])}")
	print(f"PYTHON SIDE: ddata is {ddata}")

	iview = root.createView("idata", idata)
	iview.apply(pysidre.TypeID.INT64_ID, length)
	iview.print()

	dview = root.createView("ddata", ddata)
	dview.apply(pysidre.TypeID.FLOAT64_ID, length)
	dview.print()

	assert root.getNumViews() == 2

	idata_chk = iview.getDataArray()
	assert len(idata_chk) == length
	assert np.array_equal(idata_chk, idata)

	ddata_chk = dview.getDataArray()
	assert len(ddata_chk) == length
	assert np.array_equal(ddata_chk, ddata)


# External numpy array via python
# Register with datastore then
# Query metadata using datastore API.
def test_external_int():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	iarray = np.array(range(1,11))

	view = root.createView("iarray", iarray)
	view.apply(pysidre.TypeID.INT64_ID, 10)

	assert view.isExternal() == True
	assert view.getTypeID() == pysidre.TypeID.INT64_ID
	assert view.getNumElements() == np.size(iarray)
	assert view.getNumDimensions() == 1

	extents = np.zeros(7)
	rank,extents = view.getShape(7, extents)
	assert rank == 1
	assert extents[0] == np.size(iarray)

	ipointer = view.getDataArray()
	assert np.array_equal(ipointer, iarray)


def test_external_int_3d():
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


# check other types
def test_external_float():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	darray = np.array([(i + 0.5) for i in range(1,11)])
	view = root.createView("darray", darray)
	view.apply(pysidre.TypeID.FLOAT64_ID, 10)

	assert view.getTypeID() == pysidre.TypeID.FLOAT64_ID
	assert view.getNumElements() == np.size(darray)

	dpointer = view.getDataArray()
	assert np.array_equal(dpointer, darray)


# Datastore owns a multi-dimension array.
def test_datastore_int_3d():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	extents_in = [2,3,4]

	view = root.createViewWithShapeAndAllocate("iarray", pysidre.TypeID.INT32_ID, 3, extents_in)

	ipointer = view.getDataArray()

	assert view.getTypeID() == pysidre.TypeID.INT32_ID
	assert view.getNumElements() == np.size(ipointer)
	assert view.getNumDimensions() == 3
	assert view.getNumDimensions() == ipointer.ndim

	extents = np.zeros(7)
	rank,extents = view.getShape(7, extents)
	assert rank == 3
	assert extents[0] == ipointer.shape[0]
	assert extents[1] == ipointer.shape[1]
	assert extents[2] == ipointer.shape[2]

	# Reshape as 1D using shape
	extents_in[0] = np.size(ipointer)
	view.apply(pysidre.TypeID.INT32_ID, 1, np.array([extents_in[0]]))
	assert view.getNumElements() == np.size(ipointer)

	# Reshape as 1D using length
	view.apply(pysidre.TypeID.INT32_ID, extents_in[0])
	assert view.getNumElements() == np.size(ipointer)


# Corresponding fortran test implementation needs to be fixed
# def test_save_load_external_view():
