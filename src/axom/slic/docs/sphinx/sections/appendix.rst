.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _sections/appendix:

Appendix
---------

 .. _SlicApplicationCodeExample:

Slic Application Code Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is the complete :ref:`SlicApplicationCodeExample` presented in
the :ref:`sections/getting_started` section. The code can be found in the Axom
source code under ``src/axom/slic/examples/basic/logging.cpp``.

 .. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_BASIC_EXAMPLE_BEGIN
   :end-before: SPHINX_SLIC_BASIC_EXAMPLE_END
   :language: C++
   :linenos:


.. _axomProcessAbort:

axom::utilities::processAbort()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`axomProcessAbort` function gracefully aborts the application by:

#. Calling ``abort()`` if it is a serial application.

#. Calls ``MPI_Abort()`` if the `Axom Toolkit`_ is compiled with MPI and the
   application has initialized MPI, i.e., it's a distributed MPI application.


.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst