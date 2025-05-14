# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre
import ctypes

def test_create_external_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	length = 11

	idata = list(range(length))
	ddata = [x * 2.0 for x in range(length)]

	# Cannot construct a types.CapsuleType object?
	iview = root.createView("idata", idata)
	iview.apply(pysidre.TypeID.INT_ID, length)

	assert True

# Corresponding fortran test implementation needs to be fixed
# def test_save_load_external_view():
