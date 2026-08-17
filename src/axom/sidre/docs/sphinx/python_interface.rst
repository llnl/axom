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

.. note::

   The pip/uv-installable wheel is planned and not yet available.
   This section is a placeholder for the workflow it will enable.
   Until it lands, use the build-tree helper for development builds
   or a dedicated Spack environment view for installed-package testing.

The wheel will compile only the binding code against an already-installed Axom
(located via ``CMAKE_PREFIX_PATH``); it will not build Axom or its third-party
libraries. Because a pip-built Conduit would produce a second, ABI-incompatible
``libconduit`` in the same process, the wheel will rely on the Conduit Python
module from the same Axom/Conduit build, exposed via a ``.pth`` file rather
than a PyPI install.

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
