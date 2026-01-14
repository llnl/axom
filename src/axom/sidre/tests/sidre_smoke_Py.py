# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
from conduit import Node

# Python automatically calls destructor during garbage collection
def test_create_datastore():
	ds = pysidre.DataStore()
	assert True

def test_valid_invalid():
	ds = pysidre.DataStore()

	idx = 3
	assert idx != pysidre.InvalidIndex

	name = "foo"
	assert pysidre.nameIsValid(name)

	root = ds.getRoot()
	assert root.getGroupName(idx) == pysidre.InvalidName
	assert root.getGroupIndex(name) == pysidre.InvalidIndex

	n = Node()
	n["foo"] = "bar"
	print("sdfdsfasdfasdgadgapsdgadsgpojsagp")

	print("python Node after ")
	n = ds.myFunc(n)
	print(n)
	print("sdfdsfasdfasdgadgapsdgadsgpojsagp")

def test_conduit_in_sidre_smoke():
	# make sure we are importing conduit module ok.
	n = Node()
	n["field"] = 100
	assert n["field"] == 100
	# print("sdfdsfasdfasdgadgapsdgadsgpojsagp")
	# print(type(n))
	# print("sdfdsfasdfasdgadgapsdgadsgpojsagp")
