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
   an already-installed Axom. See "Building and installing the wheel" below.


## Two ways to get the Python interface

The two paths named above differ only in *who* compiles the `_sidre` extension and
*how you import it*. They compile the **same** binding translation unit
(`src/axom/sidre/nanobind_sidre.cpp`) under the **same** nanobind domain
(`NB_DOMAIN axom`) and ship the **same** pure-Python tree from this directory,
so `import axom.sidre` behaves identically either way. Both also stand on top of a
fully built, installed Axom plus a matching Conduit -- neither path builds Axom's
C++ libraries or its third-party libraries (see "What `uv` builds" below).

### Path A -- in-tree CMake build, imported via `PYTHONPATH`

Enable the bindings in the same CMake build that compiles Axom
(a Python interpreter must be found; currently only Sidre is bound).
The `_sidre` extension is built next to `libaxom`/`libsidre`, staged into `<build>/python/axom/sidre/`,
and installed under `AXOM_PYTHON_MODULE_INSTALL_PREFIX` (default `lib/python<X.Y>/site-packages`).
To use it, put that directory -- plus Conduit's Python-module dir and numpy -- on `PYTHONPATH`, and `import axom.sidre` works.

The build configures a convenience wrapper that assembles that environment from the spack prefixes,
so an ad hoc script "just works" without a venv:

```bash
# runs the build's Python with axom + conduit + numpy (+ mpi4py) already on PYTHONPATH
<build>/bin/run_python_with_axom.sh my_script.py
```

This is the natural path during Axom development since it doesn't require a separate packaging step,
and rebuilding Axom rebuilds the bindings in place. 
The wrapper is bash-only and does not compose with Jupyter kernels, IDE runners, or debuggers.
For those, use the uv environment in **Quick start** below.

### Path B -- thin pip/uv wheel, imported into a venv

Build and install Axom first (the normal CMake/spack path, bindings enabled),
then build a **binding-only** wheel against that install and install it into a virtual environment:

```bash
# Axom + Conduit already built and installed; compile just the bindings against them.
# Point find_package at the install with axom_DIR
uv build --wheel -C cmake.define.axom_DIR="$AXOM_INSTALL/lib/cmake" src/python
uv pip install dist/axom-*.whl          # then expose the same-build Conduit; see below
```

The wheel's `CMakeLists.txt` runs `find_package(axom CONFIG REQUIRED)` and
compiles only `nanobind_sidre.cpp` -- it consumes the install, but does not rebuild it.
This is the path for distributing or consuming the bindings in an ordinary Python environment.
The full uv workflows (wheelhouse install, editable rebuild loop, the same-build-Conduit `.pth`,
and the matching-interpreter constraint) are in "Building and installing the wheel" below.

### What `uv` builds -- and what it does not

`uv build` (through scikit-build-core) runs **only** the wheel's `CMakeLists.txt`,
which compiles the single binding TU and links it against an already-installed
`axom::sidre` and `conduit::conduit_python`. 
It does **not** build:

- **Axom's C++ libraries** -- supplied by the `find_package(axom)` install.
- **Third-party libraries** (Conduit, HDF5, RAJA, Umpire, MPI, ...) -- provisioned
  by spack and not pip-installable. Conduit especially is deliberately *not* a wheel dependency;
  its Python module reaches the venv through a `.pth` pointing at the same-build Conduit 
  (a `pip install conduit` is an unrelated, ABI-incompatible package).

By design, `uv` does not generate the Axom libraries directly: 
Axom and its TPLs come from the CMake/spack world, and `uv` adds only the thin binding layer on top.
Keeping the wheel thin makes it fast and reproducible against a known install.

## Quick start: an Axom environment with uv (optionally for Jupyter)

The happy path for *using* the bindings from an ordinary Python environment (Path B).
It assumes Axom was already built and installed with its Python bindings enabled at `$AXOM_INSTALL`,
against a Conduit at `$CONDUIT_INSTALL`.
For the per-host-config wheelhouse, the editable developer loop, and the reasoning behind each step,
see "Building and installing the wheel" below.

```bash
# 1. A venv on the SAME interpreter your Axom/Conduit were built against.
uv venv --python $(which python3)

# 2. Install the axom wheel -- either from a prebuilt per-host-config wheelhouse:
uv pip install axom --find-links /path/to/wheelhouse/<hostconfig>
#    ...or build it from the source checkout against your install:
uv pip install /path/to/axom/src/python \
  -C cmake.define.axom_DIR="$AXOM_INSTALL/lib/cmake"

# 3. Expose the *same-build* Conduit Python module (never a PyPI 'conduit') 
#    with a one-line .pth dropped into the venv's site-packages:
echo "$CONDUIT_INSTALL/python-modules" > \
  "$(uv run python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')/conduit.pth"

# 4. Verify -- no PYTHONPATH needed; the installed wheel and the .pth do the work.
uv run python -c "import axom.sidre, conduit, numpy; print(axom.__version__)"
```

That environment works anywhere the venv's interpreter runs,
e.g. plain scripts, an IDE, a debugger, and Jupyter.

### Using it in Jupyter

Because the wheel and the Conduit `.pth` live in the venv's `site-packages`,
a Jupyter kernel *running in that venv* imports `axom.sidre` natively. 
There is nothing extra to wire up. Add Jupyter to the same venv and register it as a kernel:

```bash
uv pip install jupyterlab ipykernel
uv run python -m ipykernel install --user --name axom --display-name "Axom (uv)"
uv run jupyter lab
```

Then select the **Axom (uv)** kernel. See the Sidre user guide's "Python interface page"
for a worked notebook example and the two one-line checks to run when a kernel cannot import `axom.sidre`
(`src/axom/sidre/docs/sphinx/python_interface.rst`).


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

## Building and installing the wheel

The wheel is thin: it compiles only the binding code against an already-installed Axom and Conduit.
It never builds Axom or its third-party libraries, so a wheel is specific to the
Axom install (host-config / toolchain / glibc) it was built against.

These wheels are not portable and not intended for PyPI.
They carry absolute rpaths to the install's shared libraries
and skip the `auditwheel` / `delocate` repair a redistributable `manylinux` wheel needs;
they target controlled environments -- an LC host-config, a spack view, a CI image.
Producing portable, many-platform wheels would additionally need a tool such as `cibuildwheel` 
plus a bundling/repair step, which is out of scope here.

**Pointing the build at the install.** Pass one flag,
`-C cmake.define.axom_DIR="$AXOM_INSTALL/lib/cmake"`. Conduit needs no flag of its
own in the common case: `find_package(axom CONFIG)` pulls in Conduit via
`find_dependency`, using the Conduit prefix recorded in `axom-config.cmake` when
Axom was installed. 
Pass `-C cmake.define.Conduit_DIR="$CONDUIT_INSTALL/lib/cmake/conduit"` only if
that recorded path no longer resolves (e.g. a relocated install, a different mount,
or a container path).

Do *not* use `CMAKE_PREFIX_PATH`. Under scikit-build-core (which drives `uv build` and `uv pip install`)
it is force-set to the isolated build environment -- that is how the build locates its own bundled `nanobind` -- 
so a user-supplied value would be overwritten and ignored, and `find_package(axom)` 
would fail with a "Could not find axom" message.

**One more precondition.** Build the wheel from the source tree that produced the install:
the wheel's version comes from this checkout's `src/cmake/AxomVersion.cmake`
while the extension links the installed Axom, so the build fails with an explicit message
if the two versions disagree rather than shipping a wheel whose `axom.__version__` misreports its own binary.

Two constraints apply to every workflow below:

- **Same-build Conduit.** The bindings exchange `conduit::Node`s with the `conduit` Python module
  through Conduit's C capsule API, so that module must wrap the *same* `libconduit` the install links. 
  A `pip install conduit` is unsafe (unrelated PyPI package; a pip-built Conduit yields a second,
  ABI-incompatible `libconduit`). Expose the install's own Conduit instead,
  via a one-line `.pth` file (shown below).
- **Matching interpreter.** Build with the interpreter family whose toolchain/glibc matches the host-config. 
   On LC, pin it explicitly: `uv venv --python $(which python3)`.

### 1. User install from a per-host-config wheelhouse

Prebuilt wheels live in a per-host-config wheelhouse, e.g. `/usr/workspace/<...>/wheelhouse/<hostconfig>/` 
(one directory per host-config, since a wheel is not portable across toolchains).
Consume it with `--find-links` (or a `[tool.uv.sources]` entry):

```bash
uv venv --python $(which python3)
uv pip install axom --find-links /path/to/wheelhouse/dane-gcc13
```

Then add the Conduit `.pth` and verify exactly as in **Quick start** steps 3--4.

### 2. Build from source against an existing Axom install

```bash
uv pip install /path/to/axom/src/python \
  -C cmake.define.axom_DIR="$AXOM_INSTALL/lib/cmake"
```

`$AXOM_INSTALL` is a spack/uberenv-built Axom prefix (or any `cmake --install` tree)
whose `+python` bindings were enabled. Add the Conduit `.pth` afterward, as in **Quick start** step 3,
and add `Conduit_DIR` only if Conduit has moved since Axom was installed.

### 3. Developer loop (editable, rebuild-on-import)

nanobind's recommended editable flow rebuilds the extension automatically when
you re-import it after editing the binding source. Rebuild-on-import is a
scikit-build-core *experimental* feature (`editable.rebuild=true`) and may change;
if it misbehaves, reinstall the editable wheel to force a rebuild:

```bash
uv pip install nanobind 'scikit-build-core[pyproject]'
uv pip install -e src/python --no-build-isolation \
  -C cmake.define.axom_DIR="$AXOM_INSTALL/lib/cmake" \
  -C build-dir=build/py -C editable.rebuild=true
uv run pytest src/axom/sidre/tests -k _Py
```

### Stable ABI (abi3)

By default the wheel is tagged for the exact CPython that built it.
Opt into a single abi3 wheel that serves every CPython >= 3.12 on the machine
by passing both flags together (the CMake option makes nanobind build the limited-API module;
the scikit-build-core setting sets the wheel tag, and the two must agree):

```bash
uv build --wheel \
  -C cmake.define.AXOM_PYTHON_STABLE_ABI=ON \
  -C wheel.py-api=cp312 \
  -C cmake.define.axom_DIR="$AXOM_INSTALL/lib/cmake" \
  src/python
```

Below Python 3.12 nanobind silently builds a non-stable module, so only enable this on a 3.12+ interpreter;
the build fails with an explicit message if the interpreter does not provide `Development.SABIModule`,
rather than quietly producing a mislabelled wheel. Stable ABI relaxes the Python-version coupling,
not the toolchain coupling: an abi3 wheel is still specific to the host-config it was built against.

### Free-threaded Python (abi3t)

Not built today. scikit-build-core 1.0+ can emit free-threaded stable-ABI wheels
(`abi3t`, e.g. a `cp313t` / `cp315t` tag) once the bindings and Conduit run under
a free-threaded interpreter. Revisit if/when Axom targets free-threaded Python;
no action is needed for the GIL-enabled builds above.

### MPI

Wheel metadata is static, but whether the underlying Axom is an MPI build is a build-time choice,
so the wheel cannot force the MPI dependency at install time.
When you need mpi4py (to pass a communicator to `IOManager`, or to initialize MPI), install the extra explicitly:

```bash
uv pip install 'axom[mpi]'
```

pytest lives in the `test` extra (`uv pip install 'axom[test]'`),
never in the runtime dependencies.

## Notes

- These files are installed verbatim (no template substitution). They contain no CMake-configured values.
- This directory doubles as the root of the pip/uv wheel project (`pyproject.toml` + `CMakeLists.txt`),
  which reuses the very files above, so the two delivery modes cannot drift.
- End-user instructions for installing and importing the bindings also live in the Sidre user guide's
  "Python interface" page.
