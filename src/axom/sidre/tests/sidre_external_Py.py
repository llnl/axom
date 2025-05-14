# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre

def test_create_external_view():
	ds = pysidre.DataStore()
	root = ds.getRoot()

	length = 11

	idata = list(range(length))
	ddata = [x * 2.0 for x in range(length)]

    # TODO - Cannot pass python list as void *, need special nanobind numpy handling
	# iview = root.createView("idata", idata)
	# iview.apply(pysidre.TypeID.INT_ID, length)

	# dview = root.createView("ddata", ddata)
	# dview.apply(pysidre.TypeID.DOUBLE_ID, length)

	# assert root.getNumViews() == 2
	pass



	assert True

# Corresponding fortran test implementation needs to be fixed
# def test_save_load_external_view():
