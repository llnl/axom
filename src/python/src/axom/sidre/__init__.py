# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""Python bindings for Axom's Sidre component.

This package re-exports the compiled ``axom.sidre._sidre`` extension module.
The extension is present only when Axom was configured with Sidre and Python bindings enabled.
If it is missing, importing :mod:`axom.sidre` raises a :class:`ImportError` that names the component,
rather than surfacing an opaque loader error.
"""

try:
    from . import _sidre
except ImportError as exc:  # pragma: no cover - exercised only in partial installs
    raise ImportError(
        "The 'axom.sidre' extension module ('_sidre') is not available in this "
        "installation. It is built only when Axom is configured with the Sidre "
        "component and Python bindings enabled "
        "(AXOM_ENABLE_SIDRE=ON together with nanobind). Rebuild Axom with those "
        "options, or install a build that includes them, to use axom.sidre."
    ) from exc

# Re-export the extension's public surface so ``axom.sidre.DataStore`` etc.
# resolve directly on this package.  ``_sidre.__all__`` is not defined by the
# nanobind module, so fall back to a filtered ``dir()`` that drops dunders and
# the private extension handle itself.
__version__ = _sidre.__version__

__all__ = [_name for _name in dir(_sidre) if not _name.startswith("_")]

globals().update({_name: getattr(_sidre, _name) for _name in __all__})

# ``__version__`` is conventionally public but intentionally excluded from the
# wildcard surface above (it starts with an underscore); expose it explicitly.
__all__.append("__version__")
