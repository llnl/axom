.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core utilities
******************************************************

The ``axom::utilities`` namespace contains basic commonly-used functions.
Often these started out in another Axom component and were "promoted" to Axom 
Core when they were useful in other Axom libraries. In some cases, 
``axom::utilities`` brings functionality from recent C++ standards to 
platforms restricted to older compilers.

Axom can print a self-explanatory message and provide a version string.

.. literalinclude:: ../../examples/core_utilities.cpp
   :start-after: _about_start
   :end-before: _about_end
   :language: C++

Here a function showing the usage of string and filesystem facilities available
in Axom Core.

.. literalinclude:: ../../examples/core_utilities.cpp
   :start-after: _fs_string_start
   :end-before: _fs_string_end
   :language: C++

Axom Core also includes a ``Timer`` class.  Here, we time the preceding
filesystem example snippet.

.. literalinclude:: ../../examples/core_utilities.cpp
   :start-after: _timer_start
   :end-before: _timer_end
   :language: C++

Axom provides the NumericArray class, which implements arithmetic
operations on numerical tuples.

There are several other utility functions.  Some are numerical functions such as
variations on ``clamp`` (ensure a variable is restricted to a given range) and
``swap`` (exchange the values of two variables).  There are also functions for
testing values with tolerances, such as ``isNearlyEqual`` and
``isNearlyEqualRelative``. These are provided since they work in device (GPU) 
kernel code, wheras C++ standard library routines do not. There is also 
``processAbort``, to gracefully exit an application.  For details, please see 
the `Axom Core API Documentation <../../../../doxygen/html/coretop.html>`_.
