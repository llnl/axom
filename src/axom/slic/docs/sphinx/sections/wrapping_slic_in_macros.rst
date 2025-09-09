.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _sections/wrapping_slic_in_macros:

Wrapping Slic in Macros
------------------------

The recommended way of integrating Slic into an application is to wrap the
Slic API for logging messages into a set of convenience application macros that
are used throughout the application code.

This allows the application code to:

* Centralize all use of Slic behind a thin macro layer,
* Insulate the application from API changes in Slic,
* Customize and augment the behavior of logging messages if needed, e.g.,
  provide macros that are only active when the code is compiled with debug
  symbols etc.

The primary function used to log messages is ``slic::logMessage()``, which in
its most basic form takes the following arguments:

#. The :ref:`logMessageLevel` associated with the message

#. A string corresponding to the user-supplied message

#. The name of the file where the message was emitted

#. The corresponding line number within the file where the message was emitted

There are additional variants of the ``slic::logMessage()`` function that
allow an application to specify a ``TAG`` for different types of messages, etc.
Consult the `Slic Doxygen API Documentation`_ for more details.

For example, an application, ``MYAPP``,  may want to define macros to log
``DEBUG``, ``INFO``, ``WARNING`` and ``ERROR`` messages as illustrated below

.. code-block:: c++
   :linenos:

   #define MYAPP_LOGMSG( LEVEL, msg )                                         \
   {                                                                          \
     std::ostringstream oss;                                                  \
     oss << msg;                                                              \
     slic::logMessage( LEVEL, oss.str(), __FILE__, __LINE__ );                \
   }

   #define MYAPP_ERROR( msg ) MYAPP_LOGMSG( slic::message::Error, msg )
   #define MYAPP_WARNING( msg ) MYAPP_LOGMSG( slic::message::Warning, msg )
   #define MYAPP_INFO( msg ) MYAPP_LOGMSG( slic::message::Info, msg )
   #define MYAPP_DEBUG( msg ) MYAPP_LOGMSG( slic::message::Debug, msg )

These macros can then be used in the application code as follows:

.. code-block:: c++

   MYAPP_INFO( "this is an info message")
   MYAPP_ERROR( "this is an error message" );
   ...

.. note::

   Another advantage of encapsulating the Slic API calls in macros is that this
   approach alleviates the burden from application developers to have to
   pass the ``__FILE__`` and ``__LINE__`` to the ``logMessage()`` function
   each time.

   Macros that use ``slic::logMessage()`` with a :ref:`logMessageLevel` of
   ``WARNING`` or ``ERROR`` are collective operations when used with
   MPI-aware :ref:`LogStream` instances. Consult :ref:`CollectiveSlicMacros`
   for a list of collective Axom macros.

The :ref:`SlicMacros` provide a good resource for the type of macros that an
application may want to adopt and extend. Although these macros are tailored
for use within the `Axom Toolkit`_, these are also callable by application code.


.. _SlicMacros:

Slic Macros Used in Axom
^^^^^^^^^^^^^^^^^^^^^^^^^
Slic provides a set of convenience macros that can be used to streamline
logging within an application, as summarized in the table below.

.. note::

  The :ref:`SlicMacros` are not the only interface
  to log messages with Slic. They are used within `Axom Toolkit`_ for
  convenience. Applications or libraries that adopt Slic, typically, use the
  C++ API directly, e.g., call ``slic::logMessage()`` and  wrap the
  functionality in application specific macros to better suit the requirements
  of the application.

.. _CollectiveSlicMacros:

Collective Slic Macros
""""""""""""""""""""""

A subset of SLIC macros are collective operations when used with
MPI-aware :ref:`LogStream` instances such as :ref:`SynchronizedStream`
or :ref:`LumberjackStream`.

Additionally, macros such as ``SLIC_WARNING`` and ``SLIC_CHECK`` become collective
operations when certain flags are toggled on or functions are called. Other macros
such as ``SLIC_ERROR`` and ``SLIC_ASSERT`` can be made not collective when certain
functions are called.

The table below details the built-in SLIC macros as well as some notes about when they are collective calls:

+----------------------------+------------------------------------------------+----------------------------------------------------------------------------+
| Macro                      |    Availability                                |   Collective status                                                        |
+============================+================================================+============================================================================+
| | ``SLIC_ASSERT``          | |   Only available in debug configurations     | | Collective by default.                                                   |
| | ``SLIC_ASSERT_MSG``      | |   (i.e. when `AXOM_DEBUG` is defined).       | | Collective after calling ``slic::enableAbortOnError()``.                 |
| |                          | |   Not available in device code.              | | No longer collective after calling ``slic::disableAbortOnError()``.      |
+----------------------------+------------------------------------------------+----------------------------------------------------------------------------+
| | ``SLIC_CHECK``           | |   Only available in debug configurations     | | Not collective by default.                                               |
| | ``SLIC_CHECK_MSG``       | |   (i.e. when `AXOM_DEBUG` is defined).       | | Collective after ``slic::debug::checksAreErrors`` is set to ``true``,    |
| |                          | |   Not available in device code.              | |   defaults to ``false``.                                                 |
+----------------------------+------------------------------------------------+----------------------------------------------------------------------------+
| | ``SLIC_DEBUG``           | |   Only available in debug configurations     | | Never                                                                    |
| | ``SLIC_DEBUG_IF``        | |   (i.e. when `AXOM_DEBUG` is defined)        | |                                                                          |
| | ``SLIC_DEBUG_ROOT``      | |                                              | |                                                                          |
| | ``SLIC_DEBUG_ROOT_IF``   | |                                              | |                                                                          |
+----------------------------+------------------------------------------------+----------------------------------------------------------------------------+
| | ``SLIC_INFO``            | |   Always                                     | | Never                                                                    |
| | ``SLIC_INFO_IF``         | |                                              | |                                                                          |
| | ``SLIC_INFO_ROOT``       | |                                              | |                                                                          |
| | ``SLIC_INFO_ROOT_IF``    | |                                              | |                                                                          |
| | ``SLIC_INFO_TAGGED``     | |                                              | |                                                                          |
+----------------------------+------------------------------------------------+----------------------------------------------------------------------------+
| | ``SLIC_ERROR``           | |   Always                                     | | Collective by default.                                                   |
| | ``SLIC_ERROR_IF``        | |                                              | | Collective after calling ``slic::enableAbortOnError()``.                 |
| | ``SLIC_ERROR_ROOT``      | |                                              | | No longer collective after calling ``slic::disableAbortOnError()``       |
| | ``SLIC_ERROR_ROOT_IF``   | |                                              | |                                                                          |
+----------------------------+------------------------------------------------+----------------------------------------------------------------------------+
| | ``SLIC_WARNING``         | |   Always                                     | | Not collective by default.                                               |
| | ``SLIC_WARNING_IF``      | |                                              | | Collective after calling ``slic::enableAbortOnWarning()``.               |
| | ``SLIC_WARNING_ROOT``    | |                                              | | No longer collective after calling ``slic::disableAbortOnWarning()``     |
| | ``SLIC_WARNING_ROOT_IF`` | |                                              | |                                                                          |
+----------------------------+------------------------------------------------+----------------------------------------------------------------------------+

Doxygen generated API documentation on Macros can be found here: `SLIC Macros <../../../../../doxygen/html/slic__macros_8hpp.html>`_

Consider the following rules of thumb when choosing from the above logging macros:

* The `SLIC_ABORT` and `SLIC_CHECK` macros are typically used to check preconditions/postconditions of functions
  and help catch developer errors. They are only available in debug configurations (i.e. when `AXOM_DEBUG` is available).
* `SLIC_WARNING` and `SLIC_ERROR` are available in all configurations and can be used to check for conditions that might affect the results.
  They are also useful for validating user inputs.
* `SLIC_INFO` and `SLIC_DEBUG` macros are typically used to provide information about the state of an application. 
  The `SLIC_*_IF` variants can be used to conditionally log messages. `SLIC_DEBUG` macros are compiled out in non-debug configurations 
  (i.e. their messages will not get logged), while `SLIC_INFO` macros are always available.
* The `SLIC_*_ROOT` variants can help reduce logging verbosity when called in an MPI application, especially if all
  MPI ranks are expected to have the same data (for example, if a value was broadcast from one rank to all the other ranks).


.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst
