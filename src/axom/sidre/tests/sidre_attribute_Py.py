# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import numpy as np

# Global attribute values, used by multiple tests
g_name_color = "color"
g_color_none = "white"
g_color_red = "red"
g_color_blue = "blue"

g_name_animal = "animal"
g_animal_none = "human"
g_animal_cat = "cat"
g_animal_dog = "dog"

g_namea = "a"
g_nameb = "b"

g_name_dump = "dump"
g_dump_no = 0
g_dump_yes = 1

g_name_size = "size"
g_size_small = 1.2
g_size_medium = 2.3
g_size_large = 3.4

# Python equivalent of nullptr
g_attr_null = None

if pysidre.AXOM_USE_HDF5:
	g_nprotocols = 3
	g_protocols = ["sidre_json", "sidre_hdf5", "json"]
else:
	g_nprotocols = 2
	g_protocols = ["sidre_json", "json"]

g_protocol_saves_attributes = {
	"sidre_hdf5",
	"sidre_conduit_json",
	"sidre_json",
	"sidre_layout_json"
}

# ------------------------------------------------------------------------------
# Create attribute in a Datastore
def test_create_attr():
	print("Some warnings are expected in the 'create_attr' test")

	ds = pysidre.DataStore()

	nattrs = ds.getNumAttributes()
	assert nattrs == 0

	has_index = ds.hasAttribute(0)
	assert not has_index
	has_name = ds.hasAttribute(g_name_color)
	assert not has_name

	# Create string attribute
	color = ds.createAttributeString(g_name_color, g_color_none)
	assert color is not None
	assert color.getTypeID() == pysidre.TypeID.CHAR8_STR_ID

	attr_index = color.getIndex()
	assert attr_index == 0

	nattrs = ds.getNumAttributes()
	assert nattrs == 1

	has_name = ds.hasAttribute(g_name_color)
	assert has_name
	has_index = ds.hasAttribute(0)
	assert has_index

	# Try to change default to a different type.
	# Check template of setDefaultScalar.
	ok = color.setDefaultScalar(1)
	assert not ok
	ok = color.setDefaultScalar(3.14)
	assert not ok

	# Change to legal values.
	ok = color.setDefaultString("unknown")
	assert ok
	ok = color.setDefaultString("string")  # non-const string
	assert ok

	attr = ds.getAttribute(g_name_color)
	assert attr == color
	attrc = ds.getAttribute(g_name_color)
	assert attrc == color

	attr = ds.getAttribute(0)
	assert attr == color
	attrc = ds.getAttribute(0)
	assert attrc == color

	ds.destroyAttribute(color)
	nattrs = ds.getNumAttributes()
	assert nattrs == 0
	has_name = ds.hasAttribute(g_name_color)
	assert not has_name
	# At this point color points to deallocated memory

	# Create additional attributes
	dump = ds.createAttributeScalar(g_name_dump, g_dump_no)
	assert dump is not None

	attr_index = dump.getIndex()
	assert attr_index == 0

	size = ds.createAttributeScalar(g_name_size, g_size_small)
	assert size is not None

	attr_index = size.getIndex()
	assert attr_index == 1

	nattrs = ds.getNumAttributes()
	assert nattrs == 2

	ok = dump.setDefaultScalar(1)
	assert ok
	# try to change default to a different type
	ok = dump.setDefaultString(g_name_dump)
	assert not ok
	ok = dump.setDefaultString("yes")
	assert not ok
	ok = dump.setDefaultScalar(3.1415)
	assert not ok

	ds.destroyAllAttributes()
	nattrs = ds.getNumAttributes()
	assert nattrs == 0


def test_view_attr():
	print("Some warnings are expected in the 'view_attr' test")

	ds = pysidre.DataStore()

	# Create all attributes for DataStore
	attr_color = ds.createAttributeString(g_name_color, g_color_none)
	assert attr_color is not None

	attr_animal = ds.createAttributeString(g_name_animal, g_animal_none)
	assert attr_animal is not None

	root = ds.getRoot()

	# ----------------------------------------
	# Set the first attribute in a Group
	grp1 = root.createGroup("grp1")
	view1a = grp1.createView(g_namea)
	assert view1a is not None

	assert not view1a.hasAttributeValue(g_attr_null)
	assert not view1a.hasAttributeValue(attr_color)

	# Check values of unset attributes
	out1x = view1a.getAttributeString(attr_color)
	assert out1x == g_color_none

	out1y = view1a.getAttributeString(attr_animal)
	assert out1y == g_animal_none

	ok = view1a.setAttributeString(attr_color, g_color_red)
	assert ok

	assert view1a.hasAttributeValue(attr_color)

	out = view1a.getAttributeString(attr_color)
	assert out == g_color_red

	# reset attribute value
	ok = view1a.setAttributeString(attr_color, g_color_blue)
	assert ok

	out1b = view1a.getAttributeString(attr_color)
	assert out1b == g_color_blue

	# Check second, unset attribute. Should be default value
	assert not view1a.hasAttributeValue(attr_animal)
	out1d = view1a.getAttributeString(attr_animal)
	assert out1d == g_animal_none

	# Now set second attribute
	ok = view1a.setAttributeString(attr_animal, g_animal_dog)
	assert ok

	out1c = view1a.getAttributeString(attr_animal)
	assert out1c == g_animal_dog

	# ----------------------------------------
	# Set the second attribute in a Group
	grp2 = root.createGroup("grp2")

	view2a = grp2.createView(g_namea)
	assert view2a is not None

	assert not view2a.hasAttributeValue(attr_color)
	assert not view2a.hasAttributeValue(attr_animal)

	ok = view2a.setAttributeString(attr_animal, g_animal_dog)
	assert ok

	assert not view2a.hasAttributeValue(attr_color)
	assert view2a.hasAttributeValue(attr_animal)

	out2a = view2a.getAttributeString(attr_animal)
	assert out2a == g_animal_dog

	# Get the first, unset, attribute
	out2b = view2a.getAttributeString(attr_color)
	assert out2b == g_color_none

	# Now set first attribute
	ok = view2a.setAttributeString(attr_color, g_color_red)
	assert ok

	assert view2a.hasAttributeValue(attr_color)
	assert view2a.hasAttributeValue(attr_animal)

	out2c = view2a.getAttributeString(attr_color)
	assert out2c == g_color_red

	# Try to get a scalar from string
	novalue = view2a.getAttributeScalarInt(attr_color)
	assert novalue == 0

	# ----------------------------------------
	# Set attribute on second View in a Group
	grp3 = root.createGroup("grp3")
	view3a = grp3.createView(g_namea)
	assert view3a is not None
	view3b = grp3.createView(g_nameb)
	assert view3b is not None

	ok = view3b.setAttributeString(attr_animal, g_animal_dog)
	assert ok

	assert not view3b.hasAttributeValue(attr_color)
	assert view3b.hasAttributeValue(attr_animal)

	out3a = view3b.getAttributeString(attr_animal)
	assert out3a == g_animal_dog

	# ----------------------------------------
	# Moving a view should preserve attributes
	grp4 = root.createGroup("grp4")

	grp4.moveView(view3b)

	out4a = view3b.getAttributeString(attr_animal)
	assert out4a == g_animal_dog

	# Create an attribute which will be destroyed
	view3a.setAttributeString(attr_animal, g_animal_dog)

	grp3.destroyView(g_namea)
	grp4.destroyView(g_nameb)


def test_view_int_and_double():
	print("Some warnings are expected in the 'view_int_and_double' test")

	ds = pysidre.DataStore()

	# Create all attributes for DataStore
	attr_dump = ds.createAttributeScalar(g_name_dump, g_dump_no)
	assert attr_dump is not None
	assert attr_dump.getTypeID() == pysidre.TypeID.INT32_ID

	attr_size = ds.createAttributeScalar(g_name_size, g_size_small)
	assert attr_size is not None
	assert attr_size.getTypeID() == pysidre.TypeID.FLOAT64_ID

	root = ds.getRoot()

	# ----------------------------------------
	# Create a View
	grp1 = root.createGroup("grp1")
	view1a = grp1.createView(g_namea)
	assert view1a is not None

	# Get default values
	dump = view1a.getAttributeScalarInt(attr_dump)
	assert dump == g_dump_no

	size = view1a.getAttributeScalarFloat(attr_size)
	assert size == g_size_small

	# Set values
	ok = view1a.setAttributeScalar(attr_dump, g_dump_yes)
	assert ok
	dump = -1  # clear value
	dump = view1a.getAttributeScalarInt(attr_dump)
	assert dump == g_dump_yes

	ok = view1a.setAttributeScalar(attr_size, g_size_medium)
	assert ok
	size = 0.0  # clear value
	size = view1a.getAttributeScalarFloat(attr_size)
	assert size == g_size_medium

	# Set values with incorrect types
	ok = view1a.setAttributeScalar(attr_dump, g_size_small)
	assert not ok
	ok = view1a.setAttributeString(attr_dump, g_namea)
	assert not ok
	ok = view1a.setAttributeString(attr_dump, "g_namea")
	assert not ok

	# Try to get a string from a scalar
	nostr = view1a.getAttributeString(attr_dump)
	 # In Python, might return None or empty string
	assert nostr is None or nostr == ""

	i = -1
	i = view1a.getAttributeScalarInt(g_attr_null)
	assert i == 0


def test_set_default():
	print("Some warnings are expected in the 'set_default' test")

	ds = pysidre.DataStore()

	# Create all attributes for DataStore
	attr_dump = ds.createAttributeScalar(g_name_dump, g_dump_no)
	assert attr_dump is not None
	assert attr_dump.getTypeID() == pysidre.TypeID.INT32_ID

	attr_size = ds.createAttributeScalar(g_name_size, g_size_small)
	assert attr_size is not None
	assert attr_size.getTypeID() == pysidre.TypeID.FLOAT64_ID

	root = ds.getRoot()

	# ----------------------------------------
	# Create a View
	grp1 = root.createGroup("grp1")
	view1a = grp1.createView(g_namea)
	assert view1a is not None

	# reset unset attribute 1
	assert not view1a.hasAttributeValue(attr_dump)

	ok = view1a.setAttributeToDefault(attr_dump)
	assert ok

	assert not view1a.hasAttributeValue(attr_dump)

	# Set value
	ok = view1a.setAttributeScalar(attr_dump, g_dump_yes)
	assert ok
	assert view1a.hasAttributeValue(attr_dump)

	# reset set attribute 1
	ok = view1a.setAttributeToDefault(attr_dump)
	assert ok
	assert not view1a.hasAttributeValue(attr_dump)

	# reset unset attribute 2
	assert not view1a.hasAttributeValue(attr_size)

	ok = view1a.setAttributeToDefault(attr_size)
	assert ok

	assert not view1a.hasAttributeValue(attr_size)

	# Check errors
	ok = view1a.setAttributeToDefault(g_attr_null)
	assert not ok


# ------------------------------------------------------------------------------
#  get attribute as Conduit::Node

# Requires conduit::Node information
"""
def test_as_node():
	print("Some warnings are expected in the 'as_node' test")

	ds = pysidre.DataStore()

	# Create attributes for DataStore
	attr_color = ds.createAttributeString(g_name_color, g_color_none)
	assert attr_color is not None

	attr_dump = ds.createAttributeScalar(g_name_dump, g_dump_no)
	assert attr_dump is not None

	root = ds.getRoot()

	# ----------------------------------------
	# Set the first attribute in a Group
	grp1 = root.createGroup("grp1")
	view1a = grp1.createView(g_namea)
	assert view1a is not None

	ok = view1a.setAttributeString(attr_color, g_color_red)
	assert ok

	node1 = view1a.getAttributeNodeRef(attr_color)
	assert node1.as_string() == g_color_red

	node2 = view1a.getAttributeNodeRef(attr_dump)
	assert node2.as_int() == g_dump_no

	node3 = view1a.getAttributeNodeRef(g_attr_null)
	assert node3.schema().dtype().is_empty()
"""


def test_overloads():
	print("Some warnings are expected in the 'overloads' test")

	ds = pysidre.DataStore()

	# Create string and scalar attributes
	attr_color = ds.createAttributeString(g_name_color, g_color_none)
	assert attr_color is not None
	icolor = attr_color.getIndex()
	assert icolor == 0

	attr_dump = ds.createAttributeScalar(g_name_dump, g_dump_no)
	assert attr_dump is not None
	idump = attr_dump.getIndex()
	assert idump == 1

	assert attr_color == ds.getAttribute(g_name_color)
	assert attr_color == ds.getAttribute(icolor)

	# ----------------------------------------
	root = ds.getRoot()
	view = root.createView("view1")

	# string
	ok = view.setAttributeString(attr_color, g_color_red)
	assert ok
	ok = view.setAttributeString(icolor, g_color_red)
	assert ok
	ok = view.setAttributeString(g_name_color, g_color_red)
	assert ok

	attr1a = view.getAttributeString(attr_color)
	assert attr1a == g_color_red
	attr2a = view.getAttributeString(icolor)
	assert attr2a == g_color_red
	attr3a = view.getAttributeString(g_name_color)
	assert attr3a == g_color_red

	# scalar
	ok = view.setAttributeScalar(attr_dump, g_dump_yes)
	assert ok
	ok = view.setAttributeScalar(idump, g_dump_yes)
	assert ok
	ok = view.setAttributeScalar(g_name_dump, g_dump_yes)
	assert ok

	attr1b = view.getAttributeScalarInt(attr_dump)
	assert attr1b == g_dump_yes
	attr2b = view.getAttributeScalarInt(idump)
	assert attr2b == g_dump_yes
	attr3b = view.getAttributeScalarInt(g_name_dump)
	assert attr3b == g_dump_yes

	assert view.getAttributeScalarInt(attr_dump) == g_dump_yes
	assert view.getAttributeScalarInt(idump) == g_dump_yes
	assert view.getAttributeScalarInt(g_name_dump) == g_dump_yes

	# Requires conduit::Node information
	"""
	node1 = view.getAttributeNodeRef(attr_dump)
	assert node1.as_int() == g_dump_yes
	node2 = view.getAttributeNodeRef(idump)
	assert node2.as_int() == g_dump_yes
	node3 = view.getAttributeNodeRef(g_name_dump)
	assert node3.as_int() == g_dump_yes
	"""

	assert view.hasAttributeValue(attr_dump)
	assert view.hasAttributeValue(idump)
	assert view.hasAttributeValue(g_name_dump)

	ok = view.setAttributeToDefault(attr_dump)
	assert ok
	ok = view.setAttributeToDefault(idump)
	assert ok
	ok = view.setAttributeToDefault(g_name_dump)
	assert ok

	# Attribute no longer set
	assert not view.hasAttributeValue(attr_dump)
	assert not view.hasAttributeValue(idump)
	assert not view.hasAttributeValue(g_name_dump)

	# Check some errors
	assert view.getAttributeScalarInt(g_attr_null) == 0
	assert view.getAttributeScalarInt(pysidre.InvalidIndex) == 0
	assert view.getAttributeScalarInt("noname") == 0


def test_loop_attributes():
	ds = pysidre.DataStore()

	# Create attributes for DataStore
	color = ds.createAttributeString(g_name_color, g_color_none)
	icolor = color.getIndex()
	assert icolor == 0

	dump = ds.createAttributeScalar(g_name_dump, g_dump_no)
	idump = dump.getIndex()
	assert idump == 1

	size = ds.createAttributeScalar(g_name_size, g_size_small)
	isize = size.getIndex()
	assert isize == 2

	# Loop over attribute indices in DataStore
	idx1 = ds.getFirstValidAttributeIndex()
	assert idx1 == 0
	idx2 = ds.getNextValidAttributeIndex(idx1)
	assert idx2 == 1
	idx3 = ds.getNextValidAttributeIndex(idx2)
	assert idx3 == 2
	idx4 = ds.getNextValidAttributeIndex(idx3)
	assert idx4 == pysidre.InvalidIndex
	idx5 = ds.getNextValidAttributeIndex(idx4)
	assert idx5 == pysidre.InvalidIndex

	# ----------------------------------------
	root = ds.getRoot()

	# set all attributes
	view1 = root.createView("view1")
	view1.setAttributeString(color, g_color_red)
	view1.setAttributeScalar(dump, g_dump_yes)
	view1.setAttributeScalar(size, g_size_large)

	idx1 = view1.getFirstValidAttrValueIndex()
	assert idx1 == 0
	idx2 = view1.getNextValidAttrValueIndex(idx1)
	assert idx2 == 1
	idx3 = view1.getNextValidAttrValueIndex(idx2)
	assert idx3 == 2
	idx4 = view1.getNextValidAttrValueIndex(idx3)
	assert idx4 == pysidre.InvalidIndex

	# set first attribute
	view2 = root.createView("view2")
	view2.setAttributeString(color, g_color_red)

	idx1 = view2.getFirstValidAttrValueIndex()
	assert idx1 == 0
	idx2 = view2.getNextValidAttrValueIndex(idx1)
	assert idx2 == pysidre.InvalidIndex

	# set last attribute
	view3 = root.createView("view3")
	view3.setAttributeScalar(size, g_size_large)

	idx1 = view3.getFirstValidAttrValueIndex()
	assert idx1 == 2
	idx2 = view3.getNextValidAttrValueIndex(idx1)
	assert idx2 == pysidre.InvalidIndex

	# set first and last attributes
	view4 = root.createView("view4")
	view4.setAttributeString(color, g_color_red)
	view4.setAttributeScalar(size, g_size_large)

	idx1 = view4.getFirstValidAttrValueIndex()
	assert idx1 == 0
	idx2 = view4.getNextValidAttrValueIndex(idx1)
	assert idx2 == 2
	idx3 = view4.getNextValidAttrValueIndex(idx2)
	assert idx3 == pysidre.InvalidIndex

	# no attributes
	view5 = root.createView("view5")

	idx1 = view5.getFirstValidAttrValueIndex()
	assert idx1 == pysidre.InvalidIndex
	idx2 = view5.getNextValidAttrValueIndex(idx1)
	assert idx2 == pysidre.InvalidIndex


def test_iterate_attributes():
	ds = pysidre.DataStore()

	# Create attributes for DataStore
	color = ds.createAttributeString(g_name_color, g_color_none)
	icolor = color.getIndex()
	assert icolor == 0

	dump = ds.createAttributeScalar(g_name_dump, g_dump_no)
	idump = dump.getIndex()
	assert idump == 1

	size = ds.createAttributeScalar(g_name_size, g_size_small)
	isize = size.getIndex()
	assert isize == 2

	for attr in ds.attributes():
		idx = attr.getIndex()
		if idx == 0:
			assert attr.getName() == g_name_color
		elif idx == 1:
			assert attr.getName() == g_name_dump
		elif idx == 2:
			assert attr.getName() == g_name_size
		else:
			raise AssertionError(f"Unexpected attribute: {{id:{idx}, name:'{attr.getName()}'}}")

	# Destroy all attributes and check for empty
	ds.destroyAllAttributes()
	assert ds.getNumAttributes() == 0

	for attr in ds.attributes():
		raise AssertionError(f"Expected no attributes, but found: {{id:{attr.getIndex()}, name:'{attr.getName()}'}}")


def test_save_attributes():
	# Setup test data
	idata = np.zeros(5, dtype=int)
	file_path_base = "sidre_attribute_datastore_"

	ds1 = pysidre.DataStore()
	root1 = ds1.getRoot()

	# Create attributes for DataStore
	color = ds1.createAttributeString(g_name_color, g_color_none)
	assert color is not None

	dump = ds1.createAttributeScalar(g_name_dump, g_dump_no)
	assert dump is not None

	size = ds1.createAttributeScalar(g_name_size, g_size_small)
	assert size is not None

	assert ds1.getNumAttributes() == 3

	# empty
	view1a = root1.createView("empty")
	view1a.setAttributeString(color, "color-empty")
	view1a.setAttributeScalar(dump, g_dump_yes)
	view1a.setAttributeScalar(size, g_size_small)

	# buffer
	view1b = root1.createViewAndAllocate("buffer", pysidre.TypeID.INT_ID, 5)
	bdata = view1b.getDataArray()
	view1b.setAttributeString(color, "color-buffer")
	view1b.setAttributeScalar(size, g_size_medium)

	# external
	view1c = root1.createView("external", pysidre.TypeID.INT_ID, 5, idata)
	view1c.setAttributeScalar(size, g_size_large)

	# scalar
	view1d = root1.createViewScalar("scalar", 1)
	view1d.setAttributeString(color, "color-scalar")

	# string
	view1e = root1.createViewString("string", "value")
	view1e.setAttributeString(color, "color-string")

	# empty without attributes
	root1.createView("empty-no-attributes")

	for i in range(5):
		idata[i] = i
		bdata[i] = i

	# ----------------------------------------

	# Save using all protocols
	for i in range(g_nprotocols):
		file_path = file_path_base + g_protocols[i]
		root1.save(file_path, g_protocols[i])

	# Delete ds1 reference to simulate closing the file
	del ds1

	# ----------------------------------------

	# Only restore sidre_hdf5 protocol
	for i in range(g_nprotocols):
		if g_protocols[i] != "sidre_hdf5":
			continue

		file_path = file_path_base + g_protocols[i]

		ds2 = pysidre.DataStore()
		root2 = ds2.getRoot()

		root2.load(file_path, g_protocols[i])
		assert ds2.getNumAttributes() == 3

		# Asserts requires conduit::Node information
		# Check available attributes
		attr_color = ds2.getAttribute(g_name_color)
		#assert attr_color.getDefaultNodeRef().as_string() == g_color_none

		attr_dump = ds2.getAttribute(g_name_dump)
		#assert attr_dump.getDefaultNodeRef().as_int() == g_dump_no

		attr_size = ds2.getAttribute(g_name_size)
		#assert attr_size.getDefaultNodeRef().as_double() == g_size_small

		# Check attributes assigned to Views

		view2a = root2.getView("empty")
		assert view2a.hasAttributeValue(g_name_color)
		assert view2a.hasAttributeValue(g_name_dump)
		assert view2a.hasAttributeValue(g_name_size)
		assert view2a.getAttributeString(attr_color) == "color-empty"
		assert view2a.getAttributeScalarInt(attr_dump) == g_dump_yes
		assert view2a.getAttributeScalarFloat(attr_size) == g_size_small

		view2b = root2.getView("buffer")
		assert view2b.hasAttributeValue(g_name_color)
		assert not view2b.hasAttributeValue(g_name_dump)
		assert view2b.hasAttributeValue(g_name_size)
		assert view2b.getAttributeString(attr_color) == "color-buffer"
		assert view2b.getAttributeScalarFloat(attr_size) == g_size_medium

		view2c = root2.getView("external")
		assert not view2c.hasAttributeValue(g_name_color)
		assert not view2c.hasAttributeValue(g_name_dump)
		assert view2c.hasAttributeValue(g_name_size)
		assert view2c.getAttributeScalarFloat(attr_size) == g_size_large

		view2d = root2.getView("scalar")
		assert view2d.hasAttributeValue(g_name_color)
		assert not view2d.hasAttributeValue(g_name_dump)
		assert not view2d.hasAttributeValue(g_name_size)
		assert view2d.getAttributeString(attr_color) == "color-scalar"

		view2e = root2.getView("string")
		assert view2e.hasAttributeValue(g_name_color)
		assert not view2e.hasAttributeValue(g_name_dump)
		assert not view2e.hasAttributeValue(g_name_size)
		assert view2e.getAttributeString(attr_color) == "color-string"

		view2f = root2.getView("empty-no-attributes")
		assert not view2f.hasAttributeValue(g_name_color)
		assert not view2f.hasAttributeValue(g_name_dump)
		assert not view2f.hasAttributeValue(g_name_size)

		del ds2


def test_save_by_attribute():
	idata = np.zeros(5, dtype=int)
	jdata = np.zeros(5, dtype=int)
	file_path_base = "sidre_attribute_by_attribute_"

	ds1 = pysidre.DataStore()
	root1 = ds1.getRoot()

	# Create attributes for DataStore
	dump = ds1.createAttributeScalar(g_name_dump, g_dump_no)
	assert dump is not None

	# scalar
	root1.createViewScalar("view1", 1).setAttributeScalar(dump, g_dump_yes)
	root1.createViewScalar("view2", 2)

	# Create a deep path with and without attribute
	root1.createViewScalar("grp1a/grp1b/view3", 3)
	root1.createViewScalar("grp2a/view4", 4) # make sure empty "views" not saved
	root1.createViewScalar("grp2a/grp2b/view5", 5).setAttributeScalar(dump, g_dump_yes)
	root1.createView("view6", pysidre.TypeID.INT32_ID, 5, idata).setAttributeScalar(dump, g_dump_yes)
	root1.createView("grp3a/grp3b/view7", pysidre.TypeID.INT32_ID, 5, jdata)

	for i in range(5):
		idata[i] = i
		jdata[i] = i + 10

	# ----------------------------------------

	# Save with attribute filter
	for i in range(g_nprotocols):
		file_path = file_path_base + g_protocols[i]
		root1.save(file_path, g_protocols[i], dump)

	# Delete ds1 reference to simulate closing the file
	del ds1

	# ----------------------------------------

	# Only restore sidre_hdf5 protocol
	for i in range(g_nprotocols):
		if g_protocols[i] != "sidre_hdf5":
			continue

		file_path = file_path_base + g_protocols[i]
		ds2 = pysidre.DataStore()
		root2 = ds2.getRoot()
		root2.load(file_path, g_protocols[i])

		# Only views with the dump attribute should exist
		assert root2.hasView("view1")
		assert not root2.hasView("view2")
		assert not root2.hasView("grp1a/grp1b/view3")
		assert not root2.hasView("grp2a/view4")
		assert root2.hasView("grp2a/grp2b/view5")
		assert root2.hasView("view6")
		assert not root2.hasView("grp3a/grp3b/view7")

		del ds2


def test_save_load_group_with_attributes_new_ds():

	for protocol in g_protocols:
		ext = "hdf5" if protocol.endswith("hdf5") else "json"
		filename = f"saveFile_{protocol}.{ext}"

		# Set up first datastore and save to disk
		ds1 = pysidre.DataStore()
		ds1.createAttributeScalar("attr", 10)
		ds1.createAttributeString(g_name_color, g_color_none)

		gr1 = ds1.getRoot().createGroup("gr")
		gr1.createViewScalar("scalar1", 1)
		gr1.createViewScalar("scalar2", 2).setAttributeString(g_name_color, g_color_red)
		gr1.createViewScalar("scalar3", 3).setAttributeString(g_name_color, g_color_blue)

		assert ds1.getNumAttributes() == 2
		assert (pysidre.TypeID.INT32_ID == ds1.getAttribute("attr").getTypeID() or
				pysidre.TypeID.INT64_ID == ds1.getAttribute("attr").getTypeID())
		assert pysidre.TypeID.CHAR8_STR_ID == ds1.getAttribute(g_name_color).getTypeID()

		assert not gr1.getView("scalar1").hasAttributeValue(g_name_color)

		assert gr1.getView("scalar2").hasAttributeValue(g_name_color)
		assert gr1.getView("scalar2").getAttributeString(g_name_color) == g_color_red

		assert gr1.getView("scalar3").hasAttributeValue(g_name_color)
		assert gr1.getView("scalar3").getAttributeString(g_name_color) == g_color_blue

		gr1.save(filename, protocol)

		# Check if protocol supports saving attributes
		if protocol not in g_protocol_saves_attributes:
			print(f"Skipping attribute load tests for protocol '{protocol}' -- "
				   "it doesn't support saving attributes")
			continue

		# Load second datastore from saved data
		ds2 = pysidre.DataStore()
		gr2 = ds2.getRoot().createGroup("gr")
		gr2.load(filename, protocol)

		assert ds2.getNumAttributes() == 2
		assert (pysidre.TypeID.INT32_ID == ds2.getAttribute("attr").getTypeID() or
				pysidre.TypeID.INT64_ID == ds2.getAttribute("attr").getTypeID())
		assert pysidre.TypeID.CHAR8_STR_ID == ds2.getAttribute(g_name_color).getTypeID()

		assert gr2.hasView("scalar1")
		assert not gr2.getView("scalar1").hasAttributeValue(g_name_color)

		assert gr2.hasView("scalar2")
		assert gr2.getView("scalar2").hasAttributeValue(g_name_color)
		assert gr2.getView("scalar2").getAttributeString(g_name_color) == g_color_red

		assert gr2.hasView("scalar3")
		assert gr2.getView("scalar3").hasAttributeValue(g_name_color)
		assert gr2.getView("scalar3").getAttributeString(g_name_color) == g_color_blue

		# Compare attributes between ds1 and ds2
		assert ds1.getNumAttributes() == ds2.getNumAttributes()
		assert ds1.getAttribute(g_name_color).getName() == ds2.getAttribute(g_name_color).getName()
		assert ds1.getAttribute(g_name_color).getTypeID() == ds2.getAttribute(g_name_color).getTypeID()
		
		# Requires conduit::Node information
		# assert ds1.getAttribute(g_name_color).getDefaultNodeRef().to_string() == ds2.getAttribute(g_name_color).getDefaultNodeRef().to_string()


def test_save_load_group_with_attributes_same_ds():

	for protocol in g_protocols:
		ext = "hdf5" if protocol.endswith("hdf5") else "json"
		filename = f"saveFile_{protocol}.{ext}"

		print(f"Checking attribute save/load w/ protocol '{protocol}' using file '{filename}'")

		# Create the DataStore and attributes
		ds = pysidre.DataStore()
		ds.createAttributeScalar("attr", 10)
		ds.createAttributeString(g_name_color, g_color_none)

		# Attach some attributes to views
		gr = ds.getRoot().createGroup("gr")
		gr.createViewScalar("scalar1", 1)
		gr.createViewScalar("scalar2", 2).setAttributeString(g_name_color, g_color_red)
		gr.createViewScalar("scalar3", 3).setAttributeString(g_name_color, g_color_blue)

		assert ds.getNumAttributes() == 2

		assert not gr.getView("scalar1").hasAttributeValue(g_name_color)

		assert gr.getView("scalar2").hasAttributeValue(g_name_color)
		assert gr.getView("scalar2").getAttributeString(g_name_color) == g_color_red

		assert gr.getView("scalar3").hasAttributeValue(g_name_color)
		assert gr.getView("scalar3").getAttributeString(g_name_color) == g_color_blue

		gr.save(filename, protocol)

		# Check if protocol supports saving attributes
		if protocol not in g_protocol_saves_attributes:
			print(f"Skipping attribute load tests for protocol '{protocol}' -- "
				   "it doesn't support saving attributes")
			continue

		# First load: check that things are still as expected after loading
		gr.load(filename, protocol)
		assert ds.getNumAttributes() == 2
		assert not gr.getView("scalar1").hasAttributeValue(g_name_color)
		assert gr.getView("scalar2").hasAttributeValue(g_name_color)
		assert gr.getView("scalar2").getAttributeString(g_name_color) == g_color_red
		assert gr.getView("scalar3").hasAttributeValue(g_name_color)
		assert gr.getView("scalar3").getAttributeString(g_name_color) == g_color_blue

		# Second load: check that changes to attributes get overwritten when loading
		ds.destroyAttribute("attr")
		gr.getView("scalar1").setAttributeString(g_name_color, g_color_red)
		gr.getView("scalar2").setAttributeToDefault(g_name_color)
		gr.getView("scalar3").setAttributeString(g_name_color, g_color_red)

		assert ds.getNumAttributes() == 1

		# Reload group from file; this should revert changes
		gr.load(filename, protocol)

		# Check that things are reverted after loading
		assert ds.getNumAttributes() == 2
		assert not gr.getView("scalar1").hasAttributeValue(g_name_color)
		assert gr.getView("scalar2").hasAttributeValue(g_name_color)
		assert gr.getView("scalar2").getAttributeString(g_name_color) == g_color_red
		assert gr.getView("scalar3").hasAttributeValue(g_name_color)
		assert gr.getView("scalar3").getAttributeString(g_name_color) == g_color_blue