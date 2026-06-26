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

import pytest


def _fresh_import_pysidre():
    """Import 'pysidre' with a clean module cache so its import-time
    DeprecationWarning is (re)emitted deterministically."""
    sys.modules.pop("pysidre", None)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        module = importlib.import_module("pysidre")
    deprecations = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    return module, deprecations


def _clear_axom_imports():
    # These tests swap between the real staged package and synthetic packages
    # under tmp_path; cached modules would otherwise bypass sys.path changes.
    for name in list(sys.modules):
        if name == "axom" or name.startswith("axom.") or name == "pysidre":
            sys.modules.pop(name, None)


def _sidre_init_source():
    _clear_axom_imports()
    import axom.sidre as sidre

    # Exercise the installed package initializer verbatim instead of keeping a
    # test-local copy of its import-error handling logic.
    with open(sidre.__file__, "r", encoding="utf-8") as sidre_init:
        return sidre_init.read()


def _write_fake_axom_sidre(tmp_path, monkeypatch, sidre_init_source, sidre_extension_source=None):
    # Build a minimal axom.sidre package layout. Leaving _sidre absent models a
    # component-disabled install; adding _sidre.py models a discoverable module
    # whose loader/import body fails.
    package_root = tmp_path / "axom"
    sidre_root = package_root / "sidre"
    sidre_root.mkdir(parents=True)
    (package_root / "__init__.py").write_text("", encoding="utf-8")
    (sidre_root / "__init__.py").write_text(sidre_init_source, encoding="utf-8")
    if sidre_extension_source is not None:
        (sidre_root / "_sidre.py").write_text(sidre_extension_source, encoding="utf-8")

    _clear_axom_imports()
    monkeypatch.syspath_prepend(str(tmp_path))


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


def test_axom_sidre_missing_extension_gets_component_message(tmp_path, monkeypatch):
    _write_fake_axom_sidre(tmp_path, monkeypatch, _sidre_init_source())

    with pytest.raises(ImportError) as caught:
        importlib.import_module("axom.sidre")

    assert "The 'axom.sidre' extension module ('_sidre') is not available" in str(caught.value)
    assert "AXOM_ENABLE_SIDRE=ON" in str(caught.value)


def test_axom_sidre_loader_import_error_is_not_masked(tmp_path, monkeypatch):
    # A discoverable _sidre that raises ImportError represents loader failures
    # such as missing shared libraries; those errors must remain actionable.
    _write_fake_axom_sidre(
        tmp_path,
        monkeypatch,
        _sidre_init_source(),
        "raise ImportError('libsidre_dependency_missing')\n",
    )

    with pytest.raises(ImportError) as caught:
        importlib.import_module("axom.sidre")

    assert "libsidre_dependency_missing" in str(caught.value)
    assert "extension module ('_sidre') is not available" not in str(caught.value)
