.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _sections/getting_started:

Getting Started with Slic
--------------------------

This section illustrates some of the key concepts and capabilities of Slic by
presenting a short walk-through of a C++ application. The complete
:ref:`SlicApplicationCodeExample` is included in the
:ref:`sections/appendix` section and is also available within the Slic source
code, under ``src/axom/slic/examples/basic/logging.cpp``.

This example illustrates the following concepts:

* Initializing the Slic Logging Environment,
* Prescribing the :ref:`logMessageFormat`,
* Basic logging to the console using the :ref:`GenericOutputStream`, and,
* Using some of the various :ref:`SlicMacros` to log messages.

Step 1: Add Header Includes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, the Slic header must be included to make all the Slic
functions and classes accessible to an application:

.. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_INCLUDES_BEGIN
   :end-before: SPHINX_SLIC_INCLUDES_END
   :language: C++
   :linenos:

.. note::

   All the classes and functions in Slic are encapsulated within the
   ``axom::slic`` namespace.

Step 2: Initialize Slic
^^^^^^^^^^^^^^^^^^^^^^^

Prior to logging any messages, the Slic Logging Environment is initialized
by the following:

.. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_INIT_BEGIN
   :end-before: SPHINX_SLIC_INIT_END
   :language: C++
   :linenos:

This creates the root logger instance. However, in order to log messages,
an application must first specify an output destination and optionally,
prescribe the format of the log messages. These steps are demonstrated in
the following sections.

.. warning::

   If you do not initialize Slic, Slic will call ``slic::initialize()``,
   setup a minimal configuration (perform Steps 2 through 5), and issue a warning
   message. It is recommended that you call ``slic::initialize()`` to get rid
   of the warning and perform your own configuration.


.. _slicExampleStep3:

Step 3: Set the Message Format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`logMessageFormat` is specified as a string consisting of keywords,
enclosed in ``< ... >``, that Slic knows how to interpret when assembling
the log message.

.. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_FORMAT_MSG_BEGIN
   :end-before: SPHINX_SLIC_FORMAT_MSG_END
   :language: C++
   :linenos:

For example, the ``format`` string in the code snippet above indicates that
the resulting log messages will have the following format:

* A line with the message time stamp
* A line consisting of the :ref:`logMessageLevel`, enclosed in
  brackets ``[ ]``, followed by the user-supplied message,
* A third line with the name of the file where the message was emitted and
* The corresponding line number location within the file, in the fourth line.

The format string is used in :ref:`slicExampleStep5`. Specifically, it is passed
as an argument to the :ref:`GenericOutputStream` object constructor to
prescribe the format of the messages.

See the :ref:`logMessageFormat` section for the complete list of keyword options
available that may be used to customize the format of the messsages.

.. note::

   This step is optional. If the format is not specified, a
   :ref:`defaultMessageFormat` will be used to assemble the message.

Step 4: Set Severity Level
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The severity of log messages to be captured may also be adjusted at runtime to
the desired :ref:`logMessageLevel` by calling ``slic::setLoggingMsgLevel()``.
This provides another knob that the application can use to filter the type and
level of messages to be captured.

All log messages with the specified severity level or higher are captured.
For example, the following code snippet sets the severity level to *debug*.

.. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_SET_SEVERITY_BEGIN
   :end-before: SPHINX_SLIC_SET_SEVERITY_END
   :language: C++
   :linenos:

This indicates that all log messages that are *debug* or higher
are captured otherwise, the messages are ignored. Since *debug* is the lowest
severity level, all messages will be captured in this case.

.. warning::

   All messages will be ignored until the first call to ``slic::setLoggingMsgLevel()``.

.. _slicExampleStep5:

Step 5: Register a Log Stream
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Log messages can have one or more output destination. The output destination
is specified by registering a corresponding :ref:`LogStream` object to each
:ref:`logMessageLevel`.

The following code snippet uses the :ref:`GenericOutputStream` object,
one of the :ref:`BuiltInLogStreams` provided by Slic, to specify ``std::cout``
as the output destination for messages at each :ref:`logMessageLevel`.

.. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_SET_STREAM_BEGIN
   :end-before: SPHINX_SLIC_SET_STREAM_END
   :language: C++
   :linenos:

.. note::

   Instead of calling ``slic::addStreamToAllMsgLevels()`` an application
   may use ``slic::addStreamToMsgLevel()`` that allows more fine grain
   control of how to bind :ref:`LogStream` objects to each
   :ref:`logMessageLevel`. Consult the `Slic Doxygen API Documentation`_
   for more information.

The :ref:`GenericOutputStream`,  takes two arguments in its constructor:

* A C++ ``std::ostream`` object that specifies the destination
  of messages. Consequently, output of messages can be directed to the console,
  by passing ``std::cout`` or ``std::cerr``, or to a file by passing a C++
  ``std::ofstream`` object. In this case, ``std::cout`` is specified as the
  output destination.

* A string corresponding to the :ref:`logMessageFormat`, discussed in
  :ref:`slicExampleStep3`.

.. note::

   Slic maintains ownership of all registered :ref:`LogStream` instances and
   will deallocate them when ``slic::finalize()`` is called.

Step 5.1: Tagged Log Streams
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Slic has limited support for tags, where users can bind streams
to user-defined tags. The bound streams only output messages with the
given tag, disregarding the message's :ref:`logMessageLevel`. The tagged
:ref:`LogStream` can be created by calling ``slic::addStreamToTag()``.

The following code snippet uses the :ref:`GenericOutputStream` object to
specify ``std::cout`` as the output destination for messages tagged with
the custom tag ``myTag``.

.. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_SET_TAGGED_STREAM_BEGIN
   :end-before: SPHINX_SLIC_SET_TAGGED_STREAM_END
   :language: C++
   :linenos:

Step 6: Log Messages
^^^^^^^^^^^^^^^^^^^^^

Once the output destination of messages is specified, messages can be logged
using the :ref:`SlicMacros`, as demonstrated in the code snippet
below.

.. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_LOG_MESSAGES_BEGIN
   :end-before: SPHINX_SLIC_LOG_MESSAGES_END
   :language: C++
   :linenos:

.. note::

   By default, ``SLIC_ERROR()`` will print the specified message and a stacktrace
   to the corresponding output destination and call :ref:`axomProcessAbort` to
   gracefully abort the application. This behavior can be toggled by calling
   ``slic::disableAbortOnError()``. Additionally, a custom abort function can be
   registered with ``slic::setAbortFunction()``. See the `Slic Doxygen API Documentation`_
   for more details.

.. note::

   A subset of SLIC macros are collective operations when used with
   MPI-aware :ref:`LogStream` instances such as :ref:`SynchronizedStream`
   or :ref:`LumberjackStream`. Consult :ref:`CollectiveSlicMacros`
   for a list of collective Axom macros.

Step 7: Finalize Slic
^^^^^^^^^^^^^^^^^^^^^^

Before the application terminates, the Slic Logging Environment must be
finalized, as follows:

.. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_FINALIZE_BEGIN
   :end-before: SPHINX_SLIC_FINALIZE_END
   :language: C++
   :linenos:

Calling ``slic::finalize()`` will properly deallocate the registered
:ref:`LogStream` objects and terminate the Slic Logging Environment.

.. note::

   ``slic::finalize()`` is a collective operation when used with
   MPI-aware :ref:`LogStream` instances such as :ref:`SynchronizedStream`
   or :ref:`LumberjackStream`. See the `Slic Doxygen API Documentation`_
   for more details.

Step 8: Run the Example
^^^^^^^^^^^^^^^^^^^^^^^^

After building the `Axom Toolkit`_ the :ref:`SlicApplicationCodeExample` may be
run from within the build space directory as follows:

.. code-block:: bash

   > ./example/slic_logging_ex

The resulting output should look similar to the following:

.. code-block:: bash

   Fri Apr 26 14:29:53 2019
   [DEBUG]: Here is a debug message!
   FILE=/Users/zagaris2/dev/AXOM/source/axom/src/axom/slic/examples/basic/logging.cpp
   LINE=44

   Fri Apr 26 14:29:53 2019
   [INFO]: Here is an info mesage!
   FILE=/Users/zagaris2/dev/AXOM/source/axom/src/axom/slic/examples/basic/logging.cpp
   LINE=45

   Fri Apr 26 14:29:53 2019
   [WARNING]: Here is a warning!
   FILE=/Users/zagaris2/dev/AXOM/source/axom/src/axom/slic/examples/basic/logging.cpp
   LINE=46

   Fri Apr 26 14:29:53 2019
   [ERROR]: Here is an error message!
   ** StackTrace of 3 frames **
   Frame 1: axom::slic::logErrorMessage(std::__1::basic_string<char, std::__1::char_traits<char>,
   std::__1::allocator<char> > const&, std::__1::basic_string<char, std::__1::char_traits<char>,
   std::__1::allocator<char> > const&, int)
   Frame 2: main
   Frame 3: start
   =====


   FILE=/Users/zagaris2/dev/AXOM/source/axom/src/axom/slic/examples/basic/logging.cpp
   LINE=47

   Abort trap: 6
.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst
