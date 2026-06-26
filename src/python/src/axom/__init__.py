# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""Python bindings for `LLNL Axom <https://github.com/LLNL/axom>`_.

Axom is a CS infrastructure library for high-performance computing applications.
Each bound Axom component is exposed as a submodule of this ``axom`` package (for example :mod:`axom.sidre`).
A submodule is importable only when the corresponding component was enabled in the underlying Axom build.
Importing a component that was not built raises :class:`ImportError` with a message naming the missing component.

The set of submodules present in a given installation therefore mirrors
the ``AXOM_ENABLE_<COMPONENT>`` configuration of the Axom build the bindings were compiled against.
"""

from importlib import metadata as _metadata

# ``axom`` is a regular package (it ships this ``__init__.py``), not an
# implicit namespace package.  All bound components install into this single
# package directory from one Axom build; mixing components from different
# builds is unsupported (see the build-id discussion in the bindings design
# notes).

__all__ = ["__version__"]


try:
    __version__ = _metadata.version("axom")
except _metadata.PackageNotFoundError:
    # CMake build-tree staging does not create Python distribution metadata.
    # Component modules, e.g. axom.sidre, still expose the C++ Axom version.
    __version__ = "0+unknown"
