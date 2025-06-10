.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _core-acceleration:

******************************************************
Core acceleration
******************************************************

Axom's `core` component provides several utilities for user applications that 
support execution on hardware accelerators. Axom lets users control the execution 
space (i.e. *where* code runs) using RAJA and movement between memory spaces using Umpire. 
As noted in the `RAJA documentation`_,
developers of high-performance computing applications have many options for 
running code: on the CPU, using OpenMP, or on GPU hardware accelerators, and 
the options are constantly evolving.

.. note::
   * Axom's memory management and execution space APIs have default 
     implementations when Axom is not configured with Umpire and RAJA, 
     and can be used in all Axom configurations.
   * GPU execution support for Axom features is continually exapanding. 
     Axom Core offers an interface that uses RAJA and Umpire internally and 
     provides easy access to for-loop level acceleration via the parallel-for 
     idiom, which invokes a given lambda function for every index in range.
     Additional constructs for atomics, reductions, scans, and sorts are also
     available. These are templated off an execution space and typically wrap
     RAJA functionality, though serial-only implementations are also provided.

The memory management API allows users to leverage either C++ memory functions 
or Umpire, depending on the availability of Umpire at compilation. It supports 
a simple operation set: allocate, deallocate, reallocate, and copy. If Umpire 
is used, an allocator can be specified -- 
see `Umpire documentation`_ 
for details. Note that the fall-back to C++ memory functions is automatic, 
so the same piece of code can handle standard C++ or C++ with Umpire. However, 
to use advanced features, such as accessing unified memory, Umpire must be 
enabled, otherwise errors will occur at compilation.

.. note::
   We refrain from including the output of the example in this section in the 
   documentation since the output is verbose. If you are curious, please 
   run the examples yourself.

Here is an example of using Axom's memory management tools:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _membasic_start
   :end-before: _membasic_end
   :language: C++

Here is an Axom example showing sequential execution:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _exebasic_start
   :end-before: _exebasic_end
   :language: C++

Here's the same loop from the above snippet, this time with CUDA or HIP:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _deviceexebasic_start
   :end-before: _deviceexebasic_end
   :language: C++

For more advanced functionality, users can directly call RAJA and Umpire.
See the `RAJA documentation`_
and the `Umpire documentation`_.

.. _RAJA documentation: https://raja.readthedocs.io/en/main/index.html
.. _Umpire documentation: https://umpire.readthedocs.io
