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
	dview.apply(pysidre.TypeID.DOUBLE_ID, length)
	dview.print()

	assert root.getNumViews() == 2

	idata_chk = iview.getDataIntPtr()
	assert len(idata_chk) == length
	assert np.array_equal(idata_chk, idata)

	ddata_chk = dview.getDataDoublePtr()
	assert len(ddata_chk) == length
	assert np.array_equal(ddata_chk, ddata)

# Corresponding fortran test implementation needs to be fixed
# def test_save_load_external_view():
