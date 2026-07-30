.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Python interface
******************************************************

Sidre ships a Python interface, ``axom.sidre``, that mirrors much of the C++ API,
e.g. to create a ``DataStore``, navigate ``Group`` and ``View`` objects, allocate and describe data,
and exchange data with `Conduit <https://llnl-conduit.readthedocs.io>`_ ``Node`` objects and NumPy arrays without copying.
The interface is a compiled extension generated with `nanobind <https://nanobind.readthedocs.io>`_,
which is built when Axom is configured with the Sidre component and Python bindings enabled.

.. code-block:: python

   import axom.sidre as sidre

   ds = sidre.DataStore()
   root = ds.getRoot()

   grp = root.createGroup("fields")
   view = grp.createViewAndAllocate("density", sidre.TypeID.FLOAT64_ID, 10)

   # Zero-copy NumPy view onto the buffer Sidre owns
   arr = view.getDataArray()
   arr[:] = 1.0

   print(ds.getRoot().getView("fields/density").getNumElements())   # 10

The module carries a ``__version__`` matching the Axom release, and exposes
feature flags (e.g., ``AXOM_USE_HDF5``, ``AXOM_ENABLE_MPI``) so Python code can
branch on how Axom was built.

=======================================
Getting a working ``import axom.sidre``
=======================================

How to make the interface importable depends on whether you are using an
installed Axom package or a build tree from an Axom development environment.
The two workflows are intentionally different.

Development build tree
----------------------

Axom's uberenv-generated TPL environments intentionally use ``view: false``.
Those environments are for configuring and building Axom from a worktree.
They do not make the build-tree package importable by a plain interpreter when activated.

For development builds, use CTest or the generated ``run_python_with_axom.sh`` helper.
These paths set the build-tree ``PYTHONPATH`` entries needed for Axom's staged package
and its Python runtime dependencies.

.. code-block:: bash

   $ cd build-axom
   $ ctest -R sidre_smoke_Py --output-on-failure
   $ ./bin/run_python_with_axom.sh -c "import axom.sidre as sidre; print(sidre.__version__)"

Spack environment view
----------------------

Axom declares itself a Python extension (``extends("python")``), so a spack
environment view can expose the bindings in the view's ``site-packages`` alongside their dependencies.
This is useful for testing or using an installed Axom package with a plain interpreter,
but it is not the normal Axom development-build workflow.

To use this, install Axom with the ``+python`` variant in a dedicated
environment whose ``spack.yaml`` enables a view:

.. code-block:: yaml

   spack:
     specs:
       - axom+python
     view: true
     ...

After ``spack install``, the environment's interpreter should have a working Axom Python installation:

.. code-block:: bash

   $ spack env activate .
   $ python -c "import axom.sidre, conduit, numpy; print(axom.sidre.__version__)"


pip / uv wheel (thin, external Axom)
------------------------------------

The wheel compiles only the Sidre binding against an already-installed Axom.
It is tied to that Axom install, its Conduit install, and its host-config;
it is not a portable PyPI-style wheel.

.. note::
   **Do not install Conduit from PyPI.** ``axom.sidre`` must use the same
   ``libconduit`` that Axom was built against. The PyPI packages named
   ``conduit`` and ``llnl-conduit`` do not provide that same build.

   Use the Conduit Python package from the Conduit install recorded by Axom.
   Wheels built by Axom's Python project record that path in ``conduit.pth``.

Quick start
^^^^^^^^^^^

Use an absolute ``AXOM_DIR`` pointing at the directory containing
``axom-config.cmake``, usually ``$AXOM_INSTALL/lib/cmake``.

.. code-block:: bash

   $ uv venv --python $(which python3)

   $ uv pip install /path/to/axom/src/python \
       -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake"

   $ uv run python -c "import axom.sidre, conduit, numpy; print(axom.__version__)"

Optional dependencies use the normal Python extras syntax on the local source
path. Keep the same CMake ``-C`` options used for the Axom install:

.. code-block:: bash

   $ uv pip install '/path/to/axom/src/python[mpi]' \
       -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake"

   $ uv pip install '/path/to/axom/src/python[test]' \
       -C cmake.define.AXOM_DIR="$AXOM_INSTALL/lib/cmake"

Use ``[mpi]`` for ``mpi4py`` support, ``[test]`` for ``pytest``,
or combine extras as ``'/path/to/axom/src/python[mpi,test]'``.
If the Axom wheel is already installed and you only need the optional dependency package,
installing ``mpi4py`` or ``pytest`` directly is also fine.

If ``axom.sidre`` is already installed in a venv but ``import conduit`` fails,
add the same-build Conduit Python package with one ``.pth`` file. On current LC
installs this path is usually ``$CONDUIT_INSTALL/lib/pythonX.Y/site-packages``;
the value recorded by an Axom install is ``AXOM_CONDUIT_PYTHON_MODULE_DIR`` in ``axom-config.cmake``.

.. code-block:: bash

   $ CONDUIT_PYTHON_MODULE_DIR=/path/to/conduit/install/lib/pythonX.Y/site-packages
   $ printf '%s\n' "$CONDUIT_PYTHON_MODULE_DIR" > \
       "$(uv run python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')/conduit.pth"
   $ uv run python -c "import axom.sidre, conduit; print(conduit.__file__)"

If your site publishes a host-config-specific wheelhouse, install from the path
they provide with ``uv pip install axom --find-links <wheelhouse>``.
Axom does not assume a central wheelhouse.

The installed wheel also carries a CMake host-config for downstream projects:

.. code-block:: bash

   $ cmake -C "$(uv run axom-python-config --host-config)" -S /path/to/project -B build

For build details, including MPI compiler wrappers, editable installs,
and stable ABI wheels, see ``src/python/README.md``.

Using Axom in Jupyter
^^^^^^^^^^^^^^^^^^^^^

Because the wheel and the Conduit ``.pth`` live in the venv's ``site-packages``,
a Jupyter kernel running in that venv imports ``axom.sidre`` natively -- there is
nothing extra to configure, and no need to modify ``PYTHONPATH``.
Add Jupyter to the same venv and register it as a kernel:

.. code-block:: bash

   $ uv pip install jupyterlab ipykernel
   $ uv run python -m ipykernel install --user --name axom --display-name "Axom (uv)"
   $ uv run jupyter lab

For more IDE-like completions, signature help, and hover documentation in JupyterLab,
install the language-server packages in the same venv:

.. code-block:: bash

   $ uv pip install jupyterlab-lsp 'python-lsp-server[all]'

The Axom wheel installs PEP 561 type information and generated ``.pyi`` stubs for ``axom.sidre``.
JupyterLab's LSP extension can use those stubs for richer completion and overload help
than the classic notebook frontend usually shows.

Select the **Axom (uv)** kernel, then for example:

.. code-block:: python

   import axom.sidre as sidre
   import numpy as np

   ds = sidre.DataStore()
   grp = ds.getRoot().createGroup("fields")
   view = grp.createViewAndAllocate("velocity", sidre.TypeID.FLOAT64_ID, 4)
   np.asarray(view.getDataArray())[:] = [1.0, 2.0, 3.0, 4.0]   # zero-copy view
   print(np.asarray(grp.getView("velocity").getDataArray()))

.. warning::

   Sidre currently preserves the C++ API's no-op semantics for some invalid operations.
   For example, ``grp.createGroup("foo")`` followed by another ``grp.createGroup("foo")``
   returns ``None`` for the second call unless ``accept_existing=True`` is passed.
   The related SLIC diagnostic may be written to the process stderr/log stream
   instead of appearing as a notebook cell error, so notebook code should either check for ``None``
   or use the explicit ``accept_existing`` option when reusing a group is intended.

If the kernel cannot import ``axom.sidre``, it is nearly always either the wrong kernel
(one outside the venv) or a missing Conduit ``.pth``. Check both from inside the notebook:

.. code-block:: python

   import sys; print(sys.executable)          # expect <venv>/bin/python
   import conduit; print(conduit.__file__)    # expect $CONDUIT_INSTALL/lib/pythonX.Y/site-packages/...

If the underlying Axom is an MPI build and you need to pass a communicator to
``IOManager`` (or to initialize MPI), install the ``mpi`` extra.
For a local source install, use ``uv pip install '/path/to/axom/src/python[mpi]' -C ...``
as shown above; for a prebuilt wheel from a wheelhouse,
use ``uv pip install 'axom[mpi]' --find-links <wheelhouse>``.

====================================
Working with Conduit and NumPy
====================================

Arrays returned by ``View.getDataArray`` and ``Buffer.getDataArray`` are
zero-copy NumPy views onto memory Sidre owns. The array keeps the owning Sidre
object alive for as long as the array is reachable.

.. warning::

   One sharp edge remains, and the binding cannot defend against it:
   reallocating a buffer (for example growing a view) can move the underlying storage,
   leaving any previously obtained NumPy array pointing at freed memory.
   Re-acquire arrays after any operation that may reallocate,
   exactly as you would re-slice a NumPy array after resizing its base.

The ``conduit`` Python module is a hard runtime dependency of the bindings and
must wrap the same Conduit build Axom links. It imports alongside ``axom.sidre``:

.. code-block:: python

   import axom.sidre as sidre
   from conduit import Node

   n = Node()
   n["field"] = 100
   assert n["field"] == 100

For how Sidre's on-disk layout and its in-memory hierarchy relate to the
Conduit Blueprint data model, see :doc:`sidre_conduit`.

==================================================================
Running standalone scripts: the ``run_python_with_axom.sh`` helper
==================================================================

The methods above make ``import axom.sidre`` work in a plain interpreter.
If you are not in a spack environment view and just want to run a one-off
Python script that uses Axom's Python modules, the build generates a helper script,
``run_python_with_axom.sh``, that prepends directories for the required runtime dependencies
to ``PYTHONPATH`` and then runs the interpreter:

.. code-block:: bash

   $ ./bin/run_python_with_axom.sh my_script.py
   $ ./bin/run_python_with_axom.sh -c "import axom.sidre, conduit"

The helper is a ``PYTHONPATH`` prepend and is bash-only, so it does not compose
with Jupyter kernels, IDE runners, or debuggers
