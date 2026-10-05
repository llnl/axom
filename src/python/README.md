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

2. **The pip/uv wheel.** The `pyproject.toml` and `CMakeLists.txt` beside this file
   use scikit-build-core to compile the bindings against an installed Axom.
   CMake finds that install with `find_package(axom CONFIG REQUIRED)`, and
   `wheel.packages = ["src/axom"]` includes the Python sources.

Both builds compile `src/axom/sidre/nanobind_sidre.cpp` with `NB_DOMAIN axom`
and install the same Python package. The wheel requires an existing installation
of Axom and its third-party libraries. Its `conduit` Python module must come
from the Conduit build that Axom links.

This README covers building and packaging. For installation and usage, see the
Sidre user guide's [Python interface](../axom/sidre/docs/sphinx/python_interface.rst).

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
      config.py                 <- locates the wheel-generated helpers (axom-python-config)
      py.typed                  <- PEP 561 marker (typed package)
      sidre/
        __init__.py             <- re-exports the compiled 'axom.sidre._sidre'
        __init__.pyi            <- package stub; re-exports '_sidre.pyi' for type checkers
        (_sidre.<tag>.so)       <- compiled extension, produced by the build
        (_sidre.pyi)            <- type stub, produced by the build
```

Parenthesized entries are build products and are intentionally not in the repository.

The wheel build also generates `axom/share/axom-python-host-config.cmake` and
`axom/share/axom-python-env.sh`. In a CMake installation these files are absent,
so `axom.config.has_wheel_config()` returns `False` and the path accessors raise
`FileNotFoundError`.

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

The wheel depends on the Axom install and host-config used to build it.
It does not bundle shared libraries with `auditwheel` and is not intended for PyPI.

Set `AXOM_INSTALL` to the absolute Axom install prefix and pass it as `AXOM_DIR`:

```bash
AXOM_INSTALL=/absolute/path/to/axom/install
AXOM_PYTHON=$("$AXOM_INSTALL/bin/run_python_with_axom.sh" \
  -c 'import sys; print(sys.executable)')
uv build --wheel --python "$AXOM_PYTHON" \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL" src/python
```

Use the interpreter from the Axom install. Conduit's Python package contains a
CPython extension, so another Python minor version cannot load it.

The build resolves `AXOM_DIR` to `$AXOM_DIR/lib/cmake`, and also accepts a directory
that holds `axom-config.cmake` directly. CMake's own package variable, `axom_DIR`,
takes precedence when set. Do not use `CMAKE_PREFIX_PATH` for `uv build` or `uv pip install`
since scikit-build-core uses it internally for the isolated build environment.

CMake reads the Conduit install and Python package paths from `axom-config.cmake`.
Add `Conduit_DIR` only if Axom's recorded Conduit package path no longer resolves:

```bash
uv build --wheel \
  --python "$AXOM_PYTHON" \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL" \
  -C cmake.define.Conduit_DIR="$CONDUIT_INSTALL/lib/cmake/conduit" \
  src/python
```

Add `AXOM_PYTHON_CONDUIT_MODULE_DIR` only if Axom's recorded Conduit Python package path is
missing or stale:

```bash
uv build --wheel \
  --python "$AXOM_PYTHON" \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL" \
  -C cmake.define.AXOM_PYTHON_CONDUIT_MODULE_DIR="$CONDUIT_INSTALL/python-modules" \
  src/python
```

Use `AXOM_PYTHON_CONDUIT_MODULE_DIR` for this override. `find_package(axom)` loads
`ConduitConfig.cmake`, whose `set()` shadows caller-supplied cache values for
`CONDUIT_PYTHON_MODULE_DIR`.

To match the Axom install's compiler and MPI settings, pass the host-config
used to build it:

```bash
uv build --wheel \
  --python "$AXOM_PYTHON" \
  -C cmake.args=-C \
  -C cmake.args=/absolute/path/to/host-config.cmake \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL" \
  src/python
```

Build from the source tree that produced the install. The build compares the
wheel metadata version with the installed Axom version and fails if they differ.

The wheel also installs development helpers, which report their own paths:

```bash
axom-python-config --host-config  # path to axom/share/axom-python-host-config.cmake
axom-python-config --env-script   # path to axom/share/axom-python-env.sh
```

The host-config sets Axom and Conduit paths, compilers, MPI options, and the Python
interpreter for a downstream CMake project. The environment script exports
package paths and compilers that CMake reads, plus variables describing the build.
See the "pip / uv wheel" section of the Sidre user guide for examples.

### Editable installs

With `editable.rebuild=true`, scikit-build-core rebuilds the extension when a
Python process imports it after a source change. This feature is experimental.
If the rebuild fails, rerun the editable install:

```bash
uv pip install nanobind 'scikit-build-core[pyproject]' numpy
uv pip install -e 'src/python[test]' --no-build-isolation \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL" \
  -C build-dir=build/py -C editable.rebuild=true
source .venv/bin/activate
(cd "$(mktemp -d)" && python -m pytest -o python_files='*_Py.py' "$OLDPWD/src/axom/sidre/tests/")
```

The `python_files` option lets pytest discover Axom's `*_Py.py` tests. Run them
from a scratch directory because several tests write to the current directory.

### Stable ABI (abi3)

By default the wheel targets the CPython version that built it.
With CMake >= 3.26 and CPython >= 3.12, pass both flags below to build an abi3 wheel.
`AXOM_PYTHON_STABLE_ABI` enables nanobind's limited-API module, and `wheel.py-api`
sets the matching wheel tag:

```bash
uv build --wheel \
  -C cmake.define.AXOM_PYTHON_STABLE_ABI=ON \
  -C wheel.py-api=cp312 \
  -C cmake.define.AXOM_DIR="$AXOM_INSTALL" \
  src/python
```

The build requires CMake's `Development.SABIModule` component for this option
and reports an error if it is unavailable. Use a CPython 3.12+ interpreter.
The abi3 module can run on later compatible CPython versions, but still requires
the Axom install and toolchain it was built against, plus a compatible Conduit
Python module. This project does not build free-threaded `abi3t` wheels.

`wheel.py-api` is unset in `pyproject.toml` so ordinary builds retain their
CPython version tag.

### Source distributions (sdist)

Build from a full repository checkout with `pip install ./src/python` or
`uv build --wheel src/python`. An sdist omits the binding source at
`src/axom/sidre/nanobind_sidre.cpp` and the version file at
`src/cmake/AxomVersion.cmake`, both outside this project's root.
It cannot build on its own.

### Package metadata and extras

The `mpi` extra installs `mpi4py`, and the `test` extra installs `pytest`.
The package metadata does not select extras based on Axom's build options,
so users of an MPI build must request `mpi` when needed. NumPy is a required
dependency. The generated `axom-conduit.pth` exposes Conduit's Python module.

The GitHub wheel job builds against a prebuilt Axom install and passes its
host-config through `cmake.args=-C`. It selects the `mpi` extra when the installed
`axom-config.cmake` reports `AXOM_USE_MPI=ON`.
