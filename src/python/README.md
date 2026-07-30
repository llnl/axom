[comment]: # (#################################################################)
[comment]: # (Copyright Lawrence Livermore National Security, LLC and other)
[comment]: # (Axom Project Contributors. See top-level LICENSE and COPYRIGHT)
[comment]: # (files for dates and other details.)
[comment]: #
[comment]: # (# SPDX-License-Identifier: BSD-3-Clause)
[comment]: # (#################################################################)

# Axom Python package source

This directory (`src/python/`) holds the canonical source of Axom's Python package.
It is consumed by two independent build paths that must produce the same on-disk layout:

1. **The CMake build (in tree).** When Axom is configured with a component's Python bindings enabled (currently Sidre),
   the build stages this tree into the build directory and installs it into a `site-packages`-shaped prefix.
   See `src/axom/sidre/CMakeLists.txt`: it copies the files below into `${PROJECT_BINARY_DIR}/python/`
   (so the build tree is import-ready) and installs them under `AXOM_PYTHON_MODULE_INSTALL_PREFIX`.
   The compiled extension (`_sidre`) and its type stub are emitted into this layout by the build; they are not checked in.

2. **The pip/uv wheel (out of tree).** A thin, binding-only wheel built with scikit-build-core
   (the `pyproject.toml` and `CMakeLists.txt` beside this file). It treats this directory as its
   package root (`wheel.packages = ["src/axom"]`) and compiles the binding translation unit against
   an already-installed Axom.

This file discusses contributor-facing concerns. Installing and using the bindings is documented
in the Sidre user guide's "Python interface" page (`src/axom/sidre/docs/sphinx/python_interface.rst`).

## Build paths at a glance

Both build paths install the same package layout and should expose the same Python API:

- the **same** binding translation unit (`src/axom/sidre/nanobind_sidre.cpp`)
- under the **same** nanobind domain (`NB_DOMAIN axom`)
- with the **same** pure-Python tree from this directory.

They differ in where the Axom C++ libraries come from:

- **In-tree CMake build.** Axom's normal CMake build compiles the C++ libraries,
  builds `_sidre` in the same build tree, stages the package under `<build>/python/`,
  and installs it under `AXOM_PYTHON_MODULE_INSTALL_PREFIX`.
- **Thin pip/uv wheel.** The scikit-build-core project in this directory consumes
  an already-installed Axom via `find_package(axom CONFIG REQUIRED)` and compiles
  only the Python binding module against that install.

The wheel deliberately does not build Axom, Conduit, HDF5, RAJA, Umpire, MPI, or other TPLs.
Those come from the CMake/spack side. Conduit is also not listed as a Python dependency
because `axom.sidre` must import the Python module from the same Conduit build that Axom links.

## Layout

This is a standard "src layout" Python project root:

```
src/python/
  README.md                     <- this file
  pyproject.toml                <- scikit-build-core project for the pip/uv wheel
  CMakeLists.txt                <- wheel build: finds an installed Axom, builds the extension
  src/
    axom/                       <- the 'axom' regular package
      __init__.py               <- top-level package metadata
      py.typed                  <- PEP 561 marker (typed package)
      sidre/
        __init__.py             <- re-exports the compiled 'axom.sidre._sidre'
        __init__.pyi            <- package stub; re-exports '_sidre.pyi' for type checkers
        (_sidre.<tag>.so)       <- compiled extension, produced by the build
        (_sidre.pyi)            <- type stub, produced by the build
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

## Wheel build reference

The wheel compiles Axom's Python binding against an existing Axom install.
It is specific to that install and host-config; it is not repaired with `auditwheel`
and is not intended for PyPI.

Use an absolute `AXOM_DIR` pointing at the directory containing `axom-config.cmake`,
normally `$AXOM_INSTALL/lib/cmake`:

```bash
uv build --wheel -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake" src/python
```

The underlying CMake package variable is `axom_DIR`.
`AXOM_DIR` is accepted as an Axom-conventional alias.
Do not use `CMAKE_PREFIX_PATH` for `uv build` or `uv pip install` since scikit-build-core
uses it internally for the isolated build environment.

Conduit is found through `axom-config.cmake` in the normal case.
Add `Conduit_DIR` only if Axom's recorded Conduit package path no longer resolves:

```bash
uv build --wheel \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake" \
  -C cmake.define.Conduit_DIR="$CONDUIT_INSTALL/lib/cmake/conduit" \
  src/python
```

Add `CONDUIT_PYTHON_MODULE_DIR` only if Conduit's Python package path is not
recorded by Conduit's CMake config:

```bash
uv build --wheel \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake" \
  -C cmake.define.CONDUIT_PYTHON_MODULE_DIR="$CONDUIT_INSTALL/lib/pythonX.Y/site-packages" \
  src/python
```

For MPI-enabled Axom installs, pass the same compiler and MPI wrapper family
used by the Axom build. Copy these from the Axom build's `CMakeCache.txt` or
host-config:

```bash
uv build --wheel \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake" \
  -C cmake.define.CMAKE_C_COMPILER="$AXOM_C_COMPILER" \
  -C cmake.define.CMAKE_CXX_COMPILER="$AXOM_CXX_COMPILER" \
  -C cmake.define.MPI_C_COMPILER="$AXOM_MPI_C_COMPILER" \
  -C cmake.define.MPI_CXX_COMPILER="$AXOM_MPI_CXX_COMPILER" \
  src/python
```

Build from the source tree that produced the install. The build compares the
wheel metadata version with the installed Axom version and fails if they differ.

The wheel also installs development helpers:

```bash
axom-python-config --host-config  # path to axom/share/axom-python-host-config.cmake
axom-python-config --env-script   # path to axom/share/axom-python-env.sh
```

Use the host-config to seed downstream CMake projects with the same Axom,
Conduit, compiler, MPI and Python settings used by the wheel:

```bash
cmake -C "$(axom-python-config --host-config)" -S <source-dir> -B <build-dir>
```

### Developer loop (editable, rebuild-on-import)

nanobind's recommended editable flow rebuilds the extension automatically when
you re-import it after editing the binding source. Rebuild-on-import is a
scikit-build-core *experimental* feature (`editable.rebuild=true`) and may change;
if it misbehaves, reinstall the editable wheel to force a rebuild:

```bash
uv pip install nanobind 'scikit-build-core[pyproject]'
uv pip install -e src/python --no-build-isolation \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake" \
  -C build-dir=build/py -C editable.rebuild=true
(cd "$(mktemp -d)" && uv run --project "$OLDPWD" \
   pytest -o python_files='*_Py.py' "$OLDPWD/src/axom/sidre/tests/")
```

Note that Axom's Python tests are named `*_Py.py`, which pytest's default `python_files` patterns do not match
and several tests write output files into the current directory, so we run them from a scratch directory.

### Stable ABI (abi3)

By default the wheel is tagged for the exact CPython that built it.
Opt into a single abi3 wheel that serves every CPython >= 3.12 on the machine
by passing both flags together (the CMake option makes nanobind build the limited-API module;
the scikit-build-core setting sets the wheel tag, and the two must agree):

```bash
uv build --wheel \
  -C cmake.define.AXOM_PYTHON_STABLE_ABI=ON \
  -C wheel.py-api=cp312 \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake" \
  src/python
```

Below Python 3.12 nanobind silently builds a non-stable module, so only enable this on a 3.12+ interpreter;
the build fails with an explicit message if the interpreter does not provide `Development.SABIModule`,
rather than quietly producing a mislabelled wheel. Stable ABI relaxes the Python-version coupling,
not the toolchain coupling: an abi3 wheel is still specific to the host-config it was built against.
Free-threaded (`abi3t`) wheels are not built today; scikit-build-core 1.0+ can emit those tags once the
bindings and Conduit run under a free-threaded interpreter.

### Package metadata and extras

Wheel metadata is static, but whether the underlying Axom is an MPI build is a build-time choice,
so the wheel cannot force MPI dependencies at install time.
The `mpi` extra declares `mpi4py`, and the `test` extra declares `pytest`.
Runtime dependencies intentionally stay minimal: `numpy` is required,
while Conduit's Python module is exposed by the generated `conduit.pth` file.

## Notes

- This directory doubles as the root of the pip/uv wheel project (`pyproject.toml` + `CMakeLists.txt`),
  which reuses these files.
