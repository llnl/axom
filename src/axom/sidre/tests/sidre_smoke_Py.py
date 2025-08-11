# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre

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