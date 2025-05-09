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

.. _SlicMacros:

Slic Macros Used in Axom
^^^^^^^^^^^^^^^^^^^^^^^^^
Slic provides a set of convenience macros that can be used to streamline
logging within an application.

.. note::

  The :ref:`SlicMacros` are not the only interface
  to log messages with Slic. They are used within the `Axom Toolkit`_ for
  convenience. Applications or libraries that adopt Slic, typically, use the
  C++ API directly, e.g., call ``slic::logMessage()`` and  wrap the
  functionality in application specific macros to better suit the requirements
  of the application.

.. _CollectiveSlicMacros:

Collective Slic Macros
^^^^^^^^^^^^^^^^^^^^^^^^^
A subset of SLIC macros are collective operations when used with
MPI-aware :ref:`LogStream` instances such as :ref:`SynchronizedStream`
or :ref:`LumberjackStream`.

Additionally, macros such as ``SLIC_WARNING`` and ``SLIC_CHECK`` become collective
operations when certain flags are toggled on or functions are called. Other macros
such as ``SLIC_ERROR`` and ``SLIC_ASSERT`` can be made not collective when certain
functions are called.

The table below details which SLIC macros are collective:

+----------------------------+----------------------------------------------------------------------------+
| Macro                      | Collective                                                                 |
+============================+============================================================================+
| | ``SLIC_ASSERT``          | | Collective by default.                                                   |
| | ``SLIC_ASSERT_MSG``      | | Collective after calling ``slic::enableAbortOnError()``.                 |
| |                          | | No longer collective after calling ``slic::disableAbortOnError()``.      |
+----------------------------+----------------------------------------------------------------------------+
| | ``SLIC_CHECK``           | | Not collective by default.                                               |
| | ``SLIC_CHECK_MSG``       | | Collective after ``slic::debug::checksAreErrors`` is set to ``true``,    |
| |                          | |   defaults to ``false``.                                                 |
+----------------------------+----------------------------------------------------------------------------+
| | ``SLIC_DEBUG``           | | Never                                                                    |
| | ``SLIC_DEBUG_IF``        | |                                                                          |
| | ``SLIC_DEBUG_ROOT``      | |                                                                          |
| | ``SLIC_DEBUG_ROOT_IF``   | |                                                                          |
+----------------------------+----------------------------------------------------------------------------+
| | ``SLIC_INFO``            | | Never                                                                    |
| | ``SLIC_INFO_IF``         | |                                                                          |
| | ``SLIC_INFO_ROOT``       | |                                                                          |
| | ``SLIC_INFO_ROOT_IF``    | |                                                                          |
| | ``SLIC_INFO_TAGGED``     | |                                                                          |
+----------------------------+----------------------------------------------------------------------------+
| | ``SLIC_ERROR``           | | Collective by default.                                                   |
| | ``SLIC_ERROR_IF``        | | Collective after calling ``slic::enableAbortOnError()``.                 |
| | ``SLIC_ERROR_ROOT``      | | No longer collective after calling ``slic::disableAbortOnError()``       |
| | ``SLIC_ERROR_ROOT_IF``   | |                                                                          |
+----------------------------+----------------------------------------------------------------------------+
| | ``SLIC_WARNING``         | | Not collective by default.                                               |
| | ``SLIC_WARNING_IF``      | | Collective after calling ``slic::enableAbortOnWarning()``.               |
| | ``SLIC_WARNING_ROOT``    | | No longer collective after calling ``slic::disableAbortOnWarning()``     |
| | ``SLIC_WARNING_ROOT_IF`` | |                                                                          |
+----------------------------+----------------------------------------------------------------------------+

Doxygen generated API documentation on Macros can be found here: `SLIC Macros <../../../../doxygen/html/slic__macros_8hpp.html>`_

.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst