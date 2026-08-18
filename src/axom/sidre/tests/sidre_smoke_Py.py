# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import axom.sidre as sidre
from conduit import Node


# Python automatically calls destructor during garbage collection
def test_create_datastore():
    ds = sidre.DataStore()
    assert True


def test_valid_invalid():
    ds = sidre.DataStore()

    idx = 3
    assert idx != sidre.InvalidIndex

    name = "foo"
    assert sidre.nameIsValid(name)

    root = ds.getRoot()
    assert root.getGroupName(idx) == sidre.InvalidName
    assert root.getGroupIndex(name) == sidre.InvalidIndex


def test_conduit_in_sidre_smoke():
    # make sure we can import conduit module
    n = Node()
    n["field"] = 100
    assert n["field"] == 100
