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

# ``axom`` is a regular package (it ships this ``__init__.py``), not an
# implicit namespace package.  All bound components install into this single
# package directory from one Axom build; mixing components from different
# builds is unsupported (see the build-id discussion in the bindings design
# notes).

__all__ = ["__version__"]


def _discover_version() -> str:
    """Return the Axom version string.

    The version is owned by the C++ build (``AXOM_VERSION_FULL`` in  ``axom/config.hpp``)
    and surfaced on each extension module's ``__version__`` attribute.
    We read it from the ``sidre`` extension when present so there is a single source of truth.
    If no component extension is importable (an unusual, effectively content-free install) 
    we fall back to a sentinel rather than failing the package import.
    """
    try:
        from axom.sidre import _sidre  # noqa: WPS433 (local import is intentional)

        return _sidre.__version__
    except Exception:  # pragma: no cover - defensive; see docstring
        return "0+unknown"


__version__ = _discover_version()
