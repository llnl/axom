# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

# Type stub for the ``axom.sidre`` package.
#
# At runtime ``__init__.py`` re-exports the compiled ``axom.sidre._sidre``
# extension's public surface dynamically (via ``globals().update(...)``), which
# a static type checker cannot follow. This stub mirrors that re-export
# statically: ``from ._sidre import *`` pulls the typed declarations from the
# generated, adjacent ``_sidre.pyi`` so that ``axom.sidre.DataStore`` etc.
# resolve for mypy/pyright. Keep this in sync with the re-export logic in
# ``__init__.py``; the runtime module is the source of truth.

from ._sidre import *  # noqa: F401,F403

# ``__version__`` is conventionally public but excluded from the wildcard
# surface (it starts with an underscore), so re-export it explicitly here,
# matching ``__init__.py``.
__version__: str
