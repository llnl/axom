# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import numpy as np

# This example file is based on sidre_createdatastore.cpp.
# The python interface is a work-in-progress, and does not yet support
# all the features in the C++ source.

def create_datastore(region):
	ds = pysidre.DataStore()
	root = ds.getRoot()

	# Create two attributes
	ds.createAttributeScalar("vis", 0)
	ds.createAttributeScalar("restart", 1)

	# Create group children of root group
	state = root.createGroup("state")
	nodes = root.createGroup("nodes")
	fields = root.createGroup("fields")

	# Populate "state" group
	state.createViewScalar("cycle", 25)
	state.createViewScalar("time", 1.2562e-2)
	state.createViewString("name", "sample_20171206_a")

	N = 16
	nodecount = N * N * N
	eltcount = (N - 1) * (N - 1) * (N - 1)

	"""
	Populate "nodes" group

	"x", "y", and "z" are three views into a shared Sidre buffer object that
	holds 3 * nodecount doubles.  These views might describe the location of
	each node in a 16 x 16 x 16 hexahedron mesh.  Each view is described by
	number of elements, offset, and stride into that data.
	"""
	buff = ds.createBuffer(pysidre.TypeID.DOUBLE_ID, 3 * nodecount).allocate()
	nodes.createView("x", buff).apply(pysidre.TypeID.DOUBLE_ID, nodecount, 0, 3)
	nodes.createView("y", buff).apply(pysidre.TypeID.DOUBLE_ID, nodecount, 1, 3)
	nodes.createView("z", buff).apply(pysidre.TypeID.DOUBLE_ID, nodecount, 2, 3)

	"""
	Populate "fields" group

	"temp" is a view into a buffer that is not shared with another View.
	In this case, the data Buffer is allocated directly through the View
	object.  Likewise with "rho."  Both Views have the default offset (0)
	and stride (1).  These Views could point to data associated with
	each of the 15 x 15 x 15 hexahedron elements defined by the nodes above.
	"""
	temp = fields.createViewAndAllocate("temp", pysidre.TypeID.DOUBLE_ID, eltcount)
	rho = fields.createViewAndAllocate("rho", pysidre.TypeID.DOUBLE_ID, eltcount)

	# Explicitly set values for the "vis" Attribute on the "temp" and "rho" buffers.
	temp.setAttributeScalar("vis", 1)
	rho.setAttributeScalar("vis", 1)

	# The "fields" Group also contains a child Group "ext" which holds a pointer
	# to an externally owned integer array.  Although Sidre does not own the
	# data, the data can still be described to Sidre.
	ext = fields.createGroup("ext")

	# numpy of int region has been passed in as a function argument.  As with "temp"
	# and "rho", view "region" has default offset and stride.
	ext.createView("region", region).apply(pysidre.TypeID.INT_ID, eltcount)

	return ds


def access_datastore(ds):
	# Retrieve Group pointers
	root = ds.getRoot()
	state = root.getGroup("state")
	nodes = root.getGroup("nodes")
	fields = root.getGroup("fields")

	# Accessing a Group that is not there gives None
	# Requesting a nonexistent View also gives None
	goofy = root.getGroup("goofy")
	if goofy == None:
		print("No such group: goofy")
	else:
		print("Something is very wrong!")

	# Access items in "state" group
	cycle = state.getView("cycle").getDataInt()
	time = state.getView("time").getDataFloat()
	name = state.getView("name").getString()

	# Access some items in "nodes" and "fields" groups
	y = nodes.getView("y").getDataArray()
	ystride = nodes.getView("y").getStride()
	temp = fields.getView("temp").getDataArray()
	region = fields.getView("ext/region").getDataArray()

	# Nudge the 3rd node, adjust temp and region of the 3rd element
	y[2 * ystride] += 0.0032
	temp[2] *= 1.0021
	region[2] = 6

	return ds

def iterate_datastore(ds):
	fill_line = "=" * 80
	print(fill_line)

	# iterate through the attributes in ds
	print("The datastore has the following attributes:")
	for attr in ds.attributes():
		print(f"* [{attr.getIndex()}] '{attr.getName()}' of type "
			  f"{attr.getTypeID()} "

			  # Requires conduit::Node information
			  # f"and default value: {attr.getDefaultNodeRef().to_yaml()}\n"
			  )

	# iterate through the buffers in ds
	print(fill_line)
	print("The datastore has the following buffers:")
	for buff in ds.buffers():
		print(f"* [{buff.getIndex()}] "
			  f"{'Allocated' if buff.isAllocated() else 'Unallocated'} buffer with "
			  f"{buff.getNumElements()} elements of type {buff.getTypeID()} with "
			  f"{buff.getNumViews()} views")
	print(fill_line)

	# iterate through the groups of the root group
	print("The root group has the following groups:")
	for grp in ds.getRoot().groups():
		print(f"* [{grp.getIndex()}] '{grp.getName()}' with "
			  f"{grp.getNumGroups()} groups and {grp.getNumViews()} views")
	print(fill_line)

	# iterate through the views of the 'state' group
	print("The 'state' group has the following views:")
	for view in ds.getRoot().getGroup("state").views():
		print(f"* [{view.getIndex()}] '{view.getName()}' -- "
			  f"{'Allocated' if view.isAllocated() else 'Unallocated'} view of type "
			  f"{view.getTypeID()} and {view.getNumElements()} elements")
	print(fill_line)


if __name__=="__main__":
	region = np.zeros(3375, dtype = int)
	ds = create_datastore(region)
	access_datastore(ds)
	iterate_datastore(ds)
