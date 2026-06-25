# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Tests for the deprecated 'pysidre' compatibility shim.

The Sidre bindings moved from a top-level 'pysidre' module to the 'axom.sidre' package.
'pysidre' survives as a deprecation shim that re-exports 'axom.sidre' and warns on import.
These tests check that the import keeps working, warns once, and exposes the same objects as 'axom.sidre'.
"""

import importlib
import sys
import warnings


def _fresh_import_pysidre():
    """Import 'pysidre' with a clean module cache so its import-time
    DeprecationWarning is (re)emitted deterministically."""
    sys.modules.pop("pysidre", None)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        module = importlib.import_module("pysidre")
    deprecations = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    return module, deprecations


def test_pysidre_import_warns_once():
    _module, deprecations = _fresh_import_pysidre()
    assert len(deprecations) == 1
    assert "axom.sidre" in str(deprecations[0].message)


def test_pysidre_reexports_axom_sidre():
    import axom.sidre as sidre

    pysidre, _ = _fresh_import_pysidre()

    # Core symbols resolve, and to the *same* objects as axom.sidre.
    assert pysidre.DataStore is sidre.DataStore
    assert pysidre.InvalidIndex == sidre.InvalidIndex
    assert pysidre.__version__ == sidre.__version__


def test_pysidre_datastore_roundtrip():
    pysidre, _ = _fresh_import_pysidre()
    ds = pysidre.DataStore()
    root = ds.getRoot()
    grp = root.createGroup("via_shim")
    assert root.hasGroup("via_shim")
    assert grp.getName() == "via_shim"
