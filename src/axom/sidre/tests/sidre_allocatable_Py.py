# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre

# Allocate array via python
# Register with datastore then 
# Query metadata using datastore API.
def test_external_allocatable_int():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# TODO Requires nanobind's numpy support to get access to underlying data instead of CapsuleType
	# TODO Requires creating a new function, createArrayView, that performs a createView followed by an apply()(???)
	# https://github.com/LLNL/axom/blob/307465db9c9f4f6653e9460954b53adb02d41e4a/src/axom/sidre/interface/c_fortran/genfsidresplicer.py#L80-L106

	# create numpy array
	# view = root.createArrayView("array", array)


def test_external_allocatable_int_3d():
	ds = pysidre.DataStore()
	root = ds.getRoot()
	# TODO

	# create 3D numpy array
	# view = root.createArrayView("array", array)

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

