# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
"""Import-behavior tests for the 'axom.sidre' package.

These check that 'axom.sidre' produces an actionable ImportError when the
compiled '_sidre' extension is absent (a component-disabled install), and that
a genuine loader failure (a discoverable '_sidre' that itself raises
ImportError) is surfaced rather than masked by the component-missing message.
"""

import importlib
from pathlib import Path
import sys

import pytest


def _clear_axom_imports():
    # These tests swap between the real staged package and synthetic packages
    # under tmp_path; cached modules would otherwise bypass sys.path changes.
    for name in list(sys.modules):
        if name == "axom" or name.startswith("axom."):
            sys.modules.pop(name, None)


def _sidre_init_source():
    # Exercise the staged package initializer verbatim instead of keeping a
    # test-local copy of its import-error handling logic. Locate the file on
    # sys.path without importing axom.sidre; re-importing the real nanobind
    # extension after removing it from sys.modules can abort in some builds.
    for entry in sys.path:
        sidre_init = Path(entry) / "axom" / "sidre" / "__init__.py"
        if sidre_init.is_file():
            return sidre_init.read_text(encoding="utf-8")

    raise RuntimeError("Could not locate axom.sidre.__init__.py on sys.path")


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
