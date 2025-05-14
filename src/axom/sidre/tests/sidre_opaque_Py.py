import pysidre

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import pysidre

from enum import Enum

#------------------------------------------------------------------------------
# Some simple types and functions used in tests
#------------------------------------------------------------------------------
class Centering(Enum):
	Zone = 1
	Node = 2
	UnknownCentering = 3

class DType(Enum):
	Double = 1
	Int = 2
	UnknownType = 3

class Extent:
	def __init__(self, lo, hi):
		self.m_ilo = lo
		self.m_ihi = hi

	def getNumPts(self, cent):
		if cent == Centering.Zone:
			return self.m_ihi - self.m_ilo + 1
		elif cent == Centering.Node:
			return self.m_ihi - self.m_ilo + 2
		else:
			return -1

class MeshVar:
	def __init__(self, cent, dtype, depth=1):
		self.m_cent = cent
		self.m_type = dtype
		self.m_depth = depth

	def getNumVals(self, extent):
		return extent.getNumPts(self.m_cent) * self.m_depth


def test_basic_inout():
	ihi_val = 9
	ds = pysidre.DataStore()
	root = ds.getRoot()
	problem_gp = root.createGroup("problem")

	ext = Extent(0, ihi_val)

	# You can't do this, even WITH nanobind's special numpy array support
	# (no support for Arbitrary Python objects):
	# https://nanobind.readthedocs.io/en/latest/ndarray.html#limitations-related-to-dtypes
	# ext_view = problem_gp.createView("ext", ext)


def test_meshvar_test():
	# Workflow similarly involves creating a View of arbitrary type
	pass