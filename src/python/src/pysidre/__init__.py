# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""Deprecated compatibility shim for the former top-level ``pysidre`` module.

Axom's Sidre Python bindings used to install as a bare top-level extension module named ``pysidre``.
They now live in the :mod:`axom.sidre` package.
This shim re-exports :mod:`axom.sidre` under the old name so that existing ``import pysidre`` code keeps working,
and emits a single :class:`DeprecationWarning` on import.

The shim will be removed in the future, and  code should ``import axom.sidre`` directly.
"""

import warnings as _warnings

_warnings.warn(
    "'pysidre' is deprecated and will be removed in a future Axom release; "
    "import 'axom.sidre' instead.",
    DeprecationWarning,
    stacklevel=2,
)

# Re-export everything axom.sidre exposes, under the legacy module name.
from axom.sidre import *  # noqa: F401,F403 (intentional re-export)
from axom.sidre import __all__ as _sidre_all
from axom.sidre import __version__  # noqa: F401

__all__ = list(_sidre_all)
