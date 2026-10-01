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
and exchange data with `Conduit <https://llnl-conduit.readthedocs.io>`_ ``Node`` objects and zero-copy NumPy arrays.

Axom builds this `nanobind <https://nanobind.readthedocs.io>`_ extension when
the Sidre component and Python bindings are enabled.

.. code-block:: python

   import axom.sidre as sidre

   ds = sidre.DataStore()
   root = ds.getRoot()

   grp = root.createGroup("fields")
   view = grp.createViewAndAllocate("density", sidre.TypeID.FLOAT64_ID, 10)

   # Modify the Sidre buffer through a NumPy array.
   arr = view.getDataArray()
   arr[:] = 1.0

   print(ds.getRoot().getView("fields/density").getNumElements())   # 10

The module's ``__version__`` matches the Axom release. Feature flags such as
``AXOM_USE_HDF5`` and ``AXOM_ENABLE_MPI`` let Python code check which options
the Axom build enabled.

========================
Importing ``axom.sidre``
========================

Use the build-tree helper during development, or install the bindings through
a Spack environment view or a pip/uv wheel.

Development build tree
----------------------

For development builds, use CTest or the generated ``run_python_with_axom.sh`` helper.
Both set ``PYTHONPATH`` to include Axom's build-tree package and its runtime dependencies.

.. code-block:: bash

   $ cd build-axom
   $ ctest -R sidre_smoke_Py --output-on-failure
   $ ./bin/run_python_with_axom.sh -c "import axom.sidre as sidre; print(sidre.__version__)"

Spack environment view
----------------------

A Spack environment view puts the installed bindings and their dependencies in
the view's ``site-packages``. Install Axom with the ``+python`` variant in an
environment whose ``spack.yaml`` enables a view:

.. code-block:: yaml

   spack:
     specs:
       - axom+python
     view: true

After ``spack install``, activate the environment and check the import:

.. code-block:: bash

   $ spack env activate .
   $ python -c "import axom.sidre, conduit, numpy; print(axom.sidre.__version__)"

pip / uv wheel
--------------

The wheel compiles Sidre's Python bindings against an installed Axom. It is tied
to that Axom install, its Conduit install, and the host-config used to build them.

.. note::
   Use the Conduit Python package from the Conduit install recorded by Axom.
   ``axom.sidre`` must use the ``libconduit`` that Axom links. The PyPI packages
   named ``conduit`` and ``llnl-conduit`` do not provide that build.
   The wheel records the package path in ``axom-conduit.pth``.

Quick start
^^^^^^^^^^^

Set ``AXOM_INSTALL`` to the absolute Axom install prefix and pass it as
``AXOM_DIR``. The prefix's ``lib/cmake`` directory must contain ``axom-config.cmake``:

.. code-block:: bash

   $ export AXOM_INSTALL=/absolute/path/to/axom/install
   $ export AXOM_PYTHON=$("$AXOM_INSTALL/bin/run_python_with_axom.sh" \
       -c 'import sys; print(sys.executable)')
   $ uv venv --python "$AXOM_PYTHON"

   $ uv pip install /path/to/axom/src/python \
       -C cmake.define.AXOM_DIR="$AXOM_INSTALL"

   $ uv run --no-project python -c "import axom.sidre, conduit, numpy; print(axom.__version__)"

Use the interpreter from the Axom install. Conduit's Python package contains a
CPython extension and must match the venv's Python minor version.

.. note::
   These examples use ``--no-project`` to keep ``uv run`` from installing a
   project it finds in the current directory or its parents.
   Inside ``src/python``, that could rebuild the Axom wheel without ``AXOM_DIR``.
   ``--no-sync`` still selects the project's environment, which may differ
   from the venv above. You can also activate the venv with
   ``source .venv/bin/activate`` and run ``python`` directly.

To install optional dependencies, add extras to the source path and keep
the CMake ``-C`` options used to build the wheel:

.. code-block:: bash

   $ uv pip install '/path/to/axom/src/python[mpi]' \
       -C cmake.define.AXOM_DIR="$AXOM_INSTALL"

   $ uv pip install '/path/to/axom/src/python[test]' \
       -C cmake.define.AXOM_DIR="$AXOM_INSTALL"

Use ``[mpi]`` to install ``mpi4py``, ``[test]`` to install ``pytest``,
or combine extras as ``'/path/to/axom/src/python[mpi,test]'``.
If the wheel is already installed, you can install ``mpi4py`` or ``pytest`` directly.

If ``axom.sidre`` is already installed in a venv but ``import conduit`` fails,
check ``axom-conduit.pth`` in the venv. If it is missing, create it with the
``AXOM_CONDUIT_PYTHON_MODULE_DIR`` path recorded in ``axom-config.cmake``:

.. code-block:: bash

   $ CONDUIT_PY_DIR=/path/to/conduit/install/python-modules
   $ printf '%s\n' "$CONDUIT_PY_DIR" > \
       "$(uv run --no-project python -c 'import sysconfig; print(sysconfig.get_paths()["platlib"])')/axom-conduit.pth"
   $ uv run --no-project python -c "import axom.sidre, conduit; print(conduit.__file__)"

If your site provides a wheelhouse for your host-config, install from it with
``uv pip install axom --find-links <wheelhouse>``.

The installed wheel also contains a CMake host-config for downstream projects:

.. code-block:: bash

   $ cmake -C "$(uv run --no-project axom-python-config --host-config)" -S /path/to/project -B build

For build details, including MPI compiler wrappers, editable installs,
and stable ABI wheels, see ``src/python/README.md``.

Using Axom in Jupyter
^^^^^^^^^^^^^^^^^^^^^

Install Jupyter in the wheel's venv and register that venv as a kernel.
The kernel can then find the wheel and ``axom-conduit.pth`` in ``site-packages``:

.. code-block:: bash

   $ uv pip install jupyterlab ipykernel
   $ uv run --no-project python -m ipykernel install --user --name axom --display-name "Axom (uv)"
   $ uv run --no-project jupyter lab

For code completion, signature help, and hover documentation in JupyterLab,
install the language-server packages in the same venv:

.. code-block:: bash

   $ uv pip install jupyterlab-lsp 'python-lsp-server[all]'

The wheel build generates ``.pyi`` stubs for ``axom.sidre`` by default and packages
them with a PEP 561 marker. JupyterLab's LSP extension reads those stubs for type
and overload information.

Select the Axom (uv) kernel and run:

.. code-block:: python

   import axom.sidre as sidre
   import numpy as np

   ds = sidre.DataStore()
   grp = ds.getRoot().createGroup("fields")
   view = grp.createViewAndAllocate("velocity", sidre.TypeID.FLOAT64_ID, 4)
   np.asarray(view.getDataArray())[:] = [1.0, 2.0, 3.0, 4.0]   # zero-copy view
   print(np.asarray(grp.getView("velocity").getDataArray()))

.. warning::

   Calling ``grp.createGroup("foo")`` twice returns ``None`` on the second call
   unless you pass ``accept_existing=True``. The SLIC diagnostic may go to stderr
   or a log instead of appearing as a notebook cell error.
   Check for ``None``, or pass ``accept_existing=True`` to reuse the group.

If the kernel cannot import ``axom.sidre``, check that it uses
the venv's Python and can find the Conduit package:

.. code-block:: python

   import sys
   print(sys.executable)    # Expect <venv>/bin/python
   import conduit
   print(conduit.__file__)  # Expect the Conduit install recorded by Axom

For an MPI build, install the ``mpi`` extra to initialize MPI or pass a
communicator to ``IOManager``. Use the local source command above,
or ``uv pip install 'axom[mpi]' --find-links <wheelhouse>`` for a prebuilt wheel.

====================================
Working with Conduit and NumPy
====================================

``View.getDataArray`` and ``Buffer.getDataArray`` return NumPy arrays without
copying the data. Changes through the array affect the underlying storage,
which may belong to Sidre or, for an external View, another owner.
The array keeps the Sidre View or Buffer and its DataStore alive while Python
holds a reference to it.

.. warning::

   Reallocating a buffer can move its storage and leave existing NumPy arrays
   pointing at freed memory. After an operation that may reallocate,
   call ``getDataArray`` again and use the new array.

Use the ``conduit`` Python module from the Conduit build that Axom links:

.. code-block:: python

   import axom.sidre as sidre
   from conduit import Node

   n = Node()
   n["field"] = 100
   assert n["field"] == 100

See :doc:`sidre_conduit` for the relationship between Sidre's file layout,
its in-memory hierarchy, and the Conduit Blueprint data model.

==========================
Running standalone scripts
==========================

To run a script against a development build, use ``run_python_with_axom.sh``
from the build directory. It adds the Axom package and dependency directories
to ``PYTHONPATH``, then runs Python with your arguments:

.. code-block:: bash

   $ ./bin/run_python_with_axom.sh my_script.py
   $ ./bin/run_python_with_axom.sh -c "import axom.sidre, conduit"

The helper requires Bash and sets ``PYTHONPATH`` only for the process it launches.
For a notebook, IDE, or debugger that launches Python directly, select an
interpreter from the Spack view or venv where you installed the bindings.
