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

The wheel compiles only the binding code against an already-installed Axom and Conduit;
it does not build Axom or its third-party libraries.
A wheel is therefore specific to the toolchain/glibc of the Axom install it was built against.
These wheels are not portable and not intended for PyPI; 
they target controlled environments (an LC host-config, a spack view, or a CI image).

.. note::
   **Do not install Conduit from PyPI.** Axom's bindings pass ``conduit::Node`` objects across the
   C++ boundary to the ``conduit`` Python module, so that module must wrap the *same* ``libconduit``
   that Axom was compiled and linked against.
   Two separate PyPI packages are easy to reach for by mistake, and neither works:

   * ``conduit`` is an **unrelated project** (a stream-transformation library for power-engineering
     analytics).
   * ``llnl-conduit`` **is** LLNL's Conduit, but installing it produces a *separate build* of the library,
     compiled by pip with its own compiler, flags and third-party configuration.
     It is unlikely to be ABI-compatible with the spack/CMake Conduit inside your Axom install,
     and using it would put two ``libconduit`` libraries in one process.

   Instead, expose the Conduit that your Axom was built against, using the one-line ``.pth`` file in step 3 below.

Quick start
^^^^^^^^^^^

Point the build at the Axom install with ``axom_DIR``.
Conduit is located transitively through Axom's own CMake config, so it normally needs no flag of its own:

.. code-block:: bash

   # 1. A venv on the same interpreter Axom and Conduit were built against.
   $ uv venv --python $(which python3)

   # 2a. Install a prebuilt wheel from a per-host-config wheelhouse...
   $ uv pip install axom --find-links /path/to/wheelhouse/<hostconfig>

   # 2b. ...or build it from a source checkout against your Axom install.
   $ uv pip install /path/to/axom/src/python \
       -C cmake.define.axom_DIR="$AXOM_INSTALL/lib/cmake"

   # 3. Expose the *same-build* Conduit Python module (not a PyPI package; see the note above).
   $ echo "$CONDUIT_INSTALL/python-modules" > \
       "$(uv run python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')/conduit.pth"

   # 4. Verify -- no PYTHONPATH and no wrapper script.
   $ uv run python -c "import axom.sidre, conduit, numpy; print(axom.__version__)"

Three details matter when building the wheel yourself:

* Use ``axom_DIR``, not ``CMAKE_PREFIX_PATH``. scikit-build-core (which backs
  ``uv build`` and ``uv pip install``) force-sets ``CMAKE_PREFIX_PATH`` to its own
  isolated build environment, so a user-supplied value would be overwritten and
  ``find_package(axom)`` would fail. A standalone ``cmake -S src/python`` has no such
  layer and can use ``CMAKE_PREFIX_PATH`` directly.
* Add ``-C cmake.define.Conduit_DIR="$CONDUIT_INSTALL/lib/cmake/conduit"`` only if
  Conduit has moved since Axom was installed: Axom's config records the Conduit
  prefix it was built against, and that recorded path is what the transitive lookup uses.
* Build from the source tree that produced the install.
  The wheel takes its version from ``src/cmake/AxomVersion.cmake`` in the checkout, 
  and the build fails with an explicit message if that disagrees with the installed Axom, 
  so a wheel can never misreport the version of the binary inside it.

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

Select the **Axom (uv)** kernel, then for example:

.. code-block:: python

   import axom.sidre as sidre
   import numpy as np

   ds = sidre.DataStore()
   grp = ds.getRoot().createGroup("fields")
   view = grp.createViewAndAllocate("velocity", sidre.TypeID.FLOAT64_ID, 4)
   np.asarray(view.getDataArray())[:] = [1.0, 2.0, 3.0, 4.0]   # zero-copy view
   print(np.asarray(grp.getView("velocity").getDataArray()))

If the kernel cannot import ``axom.sidre``, it is nearly always either the wrong kernel
(one outside the venv) or a missing Conduit ``.pth``. Check both from inside the notebook:

.. code-block:: python

   import sys; print(sys.executable)          # expect <venv>/bin/python
   import conduit; print(conduit.__file__)    # expect $CONDUIT_INSTALL/python-modules/...

If the underlying Axom is an MPI build and you need to pass a communicator to
``IOManager`` (or to initialize MPI), install the ``mpi`` extra with ``uv pip install 'axom[mpi]'``.

For the per-host-config wheelhouse convention, the editable developer loop,
and the optional stable-ABI build, see ``src/python/README.md``.

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
