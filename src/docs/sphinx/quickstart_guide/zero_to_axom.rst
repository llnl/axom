.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)


Zero to Axom: Quick install of Axom and Third Party Dependencies
================================================================

The easiest path to build Axom's dependencies is to use `Spack <https://github.com/spack/spack>`_.
This has been encapsulated using `Uberenv <https://github.com/LLNL/uberenv>`_. Uberenv helps by
doing the following:

* Pulls a blessed version of Spack locally
* If you are on a known operating system, such as the clusters on Livermore Computing (LC), or our provided docker images, we have defined Spack configuration files
  to keep Spack from building the world
* Installs our Spack packages into the local Spack
* Simplifies whole dependency build into one command

Uberenv will create a directory containing a Spack instance with the required Axom
dependencies installed.

It also generates a host-config file (``<config_dependent_name>.cmake``)
at the root of Axom repository. This host-config defines all the required information for building
Axom.

.. code-block:: bash

   $ python3 scripts/uberenv/uberenv.py

.. note::
  On LC machines, it is good practice to do the build step in parallel on a compute node
  (rather than using a login node). Here is an example command:
  ``salloc -ppdebug -N1-1 python3 scripts/uberenv/uberenv.py``

Unless otherwise specified, Spack will default to a compiler.  This is generally not a good idea when
developing large codes. To specify which compiler to use, add the compiler specification to the ``--spec`` Uberenv
command line option. Supported compiler specs can be found in the Spack compiler files in our repository:
``scripts/spack/configs/<platform>/spack.yaml``.

We currently regularly test the following Spack configuration files:

* Linux Ubuntu 22.04 (via Windows WSL 2)
* TOSS 4 (On Dane at LC)
* TOSS 4 Cray (On Tioga at LC)
* BlueOS (On Lassen at LC)

To install Axom on a new platform, it is a good idea to start with a known Spack configuration directory
(located in the Axom repo at ``scripts/spack/configs/<platform>``). The ``spack.yaml`` file
describes the compilers, associated flags required for the platform, and the
the low-level libraries on the system to prevent Spack from building the world. Documentation on
these configuration files is located in the `Spack docs <https://spack.readthedocs.io/en/latest/configuration.html>`_.

Some helpful uberenv options include :

* ``--spec=+cuda`` (build Axom with CUDA support)
* ``--spec=+devtools`` (also build the devtools with one command)
* ``--spec=%clang@10.0.0`` (build with a specific compiler as defined in the ``compiler.yaml`` file)
* ``--spack-env-file=<Path to Spack Environment file>`` (use specific Spack configuration file)
* ``--prefix=<Path>`` (required, build and install the dependencies in a particular location)

There is more thorough `uberenv documentation here <https://uberenv.readthedocs.io/en/latest/>`_.

The modifiers to the Spack specification ``spec`` can be chained together, e.g. ``--spec=%clang@10.0.0+debug+devtools``.

If you already have a Spack instance from another project that you would like to reuse,
you can do so by changing the uberenv command as follows:

.. code-block:: bash

   $ python3 scripts/uberenv/uberenv.py --upstream=</path/to/my/spack>/opt/spack


.. _using_in_your_project:

Preparing Windows WSL/Ubuntu for Axom installation
--------------------------------------------------

For faster installation of the Axom dependencies via Spack on Windows WSL/Ubuntu systems,
install CMake, MPICH, openblas, OpenGL, and the various developer tools using the following commands:

**Ubuntu 22.04**

.. code-block:: bash

   $ sudo apt-get update
   $ sudo apt-get upgrade
   $ sudo apt-get install cmake libopenblas-dev libopenblas-base mpich cppcheck doxygen libreadline-dev python3-sphinx python3-pip clang-format-10 m4
   $ sudo ln -s /usr/lib/x86_64-linux-gnu/* /usr/lib


Note that the last line is required since Spack expects the system libraries to exist in a directory named ``lib``. 

During the third party library build phase, the appropriate Spack config directory must be specified, e.g.:

.. code-block:: bash

   $ python3 scripts/uberenv/uberenv.py --spack-env-file=scripts/spack/configs/linux_ubuntu_22/spack.yaml --prefix=path/to/install/libraries



Using Axom in Your Project
--------------------------

The install includes examples that demonstrate how to use Axom
in CMake-based, BLT-based and Makefile-based build systems.

CMake-based build system example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. literalinclude:: ../../../examples/using-with-cmake/CMakeLists.txt
   :language: cmake
   :start-after: _zerotoaxom_docs_start
   :end-before: _zerotoaxom_docs_end

See:  ``examples/axom/using-with-cmake``

BLT-based build system example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. literalinclude:: ../../../examples/using-with-blt/CMakeLists.txt
   :language: cmake
   :start-after: _zerotoaxom_docs_start
   :end-before: _zerotoaxom_docs_end

See:  ``examples/axom/using-with-blt``


Makefile-based build system example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../examples/using-with-make/Makefile.in
   :language: make
   :start-after: _zerotoaxom_docs_start
   :end-before: _zerotoaxom_docs_end

See: ``examples/axom/using-with-make``
