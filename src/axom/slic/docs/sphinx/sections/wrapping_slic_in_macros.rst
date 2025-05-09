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

.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst
