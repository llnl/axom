[comment]: # (#################################################################)
[comment]: # (Copyright Lawrence Livermore National Security, LLC and other)
[comment]: # (Axom Project Contributors. See top-level LICENSE and COPYRIGHT)
[comment]: # (files for dates and other details.)
[comment]: #
[comment]: # (# SPDX-License-Identifier: BSD-3-Clause)
[comment]: # (#################################################################)

# Axom Python package source

This directory (`src/python/`) holds the canonical source of Axom's Python package, `axom`.
It is consumed by two independent build paths that must produce the same on-disk layout:

1. **The CMake build (in tree).** When Axom is configured with a component's Python bindings enabled (currently Sidre), 
   the build stages this tree into the build directory and installs it into a `site-packages`-shaped prefix.
   See `src/axom/sidre/CMakeLists.txt`: it copies the files below into `${PROJECT_BINARY_DIR}/python/` (so the build tree is import-ready)
   and installs them under `AXOM_PYTHON_MODULE_INSTALL_PREFIX`. 
   The compiled extension (`_sidre`) and its type stub are emitted into this layout by the build; they are not checked in.

2. **[planned] The pip/uv wheel (out of tree).** A thin, binding-only wheel built with scikit-build-core
   will treat this directory as its package root (`wheel.packages = ["src/axom", "src/pysidre"]` in a sibling `pyproject.toml`),
   compiling the binding translation unit against an already-installed Axom. 


## Layout

This is a standard "src layout" Python project root:

```
src/python/
  README.md                     <- this file
  src/
    axom/                       <- the 'axom' namespace package (regular package)
      __init__.py               <- package version (sourced from the extension)
      py.typed                  <- PEP 561 marker (typed package)
      sidre/
        __init__.py             <- re-exports the compiled 'axom.sidre._sidre'
        (_sidre.<tag>.so)       <- compiled extension, produced by the build
        (_sidre.pyi)            <- type stub, produced by the build
    pysidre/
      __init__.py               <- deprecation shim re-exporting 'axom.sidre'
```

Parenthesized entries are build products and are intentionally not in the repository.

Each bound Axom component installs as a submodule of the `axom` package
(`axom.sidre`, and later `axom.quest`, `axom.primal`, ...).
A submodule is importable only when its component was enabled in the underlying Axom build.

## What goes here vs. what does not

- **Here:** importable pure-Python sources that are part of the installed package:
  package `__init__.py` files, the `py.typed` marker, and any future pure-Python helpers or shims.
- **Not here:** the C++ binding code (each component's nanobind translation unit lives with that component,
  e.g. `src/axom/sidre/nanobind_sidre.cpp`),
  generated artifacts (the `.so` and `.pyi` are produced by the build),
  and tests/examples (those live under the component, e.g. `src/axom/sidre/tests/*_Py.py`).

## Notes

- These files are installed verbatim (no template substitution). They contain no CMake-configured values.
- A `pyproject.toml` for the standalone wheel is not present yet. We will add it in the future when we add the wheel. 
  Until then this directory is consumed only by the CMake build.
- End-user instructions for installing and importing the bindings currently live in the Sidre user guide's "Python interface" page.
