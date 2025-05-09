.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _gpu-porting-label:

*******************************
GPU Porting in Axom
*******************************

Axom uses the following two libraries as the main workhorses for GPU porting:

  * `RAJA <https://github.com/LLNL/RAJA>`_:
    Handles software abstractions that enable architecture and programming
    model portability for HPC applications
  * `Umpire <https://github.com/LLNL/Umpire>`_:
    Resource management library that allows the discovery, provision, and
    management of memory on machines with multiple memory devices like NUMA
    and GPUs.

From RAJA and Umpire, Axom derives a set of convenience macros and function
wrappers in the ``axom`` namespace encapsulating commonly-used RAJA and Umpire
functions, and preset execution spaces for host/device execution.

For the user's guide on using GPU utilities, see also
:ref:`Core Acceleration<core-acceleration>`.


.. _gpu-macros-label:

===============
Macros 
===============

Axom's macros can be found in the file
`axom/core/Macros.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/Macros.hpp>`_.

Most of the GPU-related macros are used to guard device code for compilation.

For guarding device code:

.. literalinclude:: ../../../axom/core/Macros.hpp
   :start-after:  _guarding_macros_start
   :end-before:  _guarding_macros_end
   :language: C++

.. note::

  * Functions called in CUDA or HIP GPU device code require the ``__device__``
    annotation.
  * Functions that will be called in device code *and* CPU host code require
    the ``__host__ __device__`` annotation.

The following code shows the macros used in Axom to apply these annotations:

.. literalinclude:: ../../../axom/core/Macros.hpp
   :start-after:  _decorating_macros_start
   :end-before:  _decorating_macros_end
   :language: C++

Below is a function that uses Axom macros to apply a ``__host__ __device__``
annotation and guard the use of a CUDA intrinsic to inside a kernel:

.. literalinclude:: ../../../axom/core/utilities/BitUtilities.hpp
   :start-after:  gpu_macros_example_start
   :end-before:  gpu_macros_example_end
   :language: C++

.. _gpu-memory-label:

===============
Memory 
===============

Axom's memory management routines can be found in the file
`axom/core/memory_management.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/memory_management.hpp>`_.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Memory Management Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Umpire has the concept of "allocators" associated with each
`memory resource type`_ (``umpire::resource::MemoryResourceType``).

To allocate memory on a particular resource, you use the ID for the allocator
associated with the ``umpire::resource::MemoryResourceType``.

You are able to set a default allocator, whereby all your memory allocations
will go on the resource associated with the allocator unless otherwise
specified:

.. literalinclude:: ../../../axom/core/memory_management.hpp
   :start-after:  _memory_management_routines_start
   :end-before:  _memory_management_routines_end
   :language: C++

.. note::

  When Axom is built without Umpire, the getters and setters shown above
  become no-ops or are undefined, while the memory allocation functions
  default to C++ standard library functions with only allocation on the
  host (CPU):

  * ``axom::allocate`` calls ``std::malloc``
  * ``axom::deallocate`` calls ``std::free``
  * ``axom::reallocate`` calls ``std::realloc``
  * ``axom::copy`` calls ``std::memcpy``


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MemorySpace
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../axom/core/memory_management.hpp
   :start-after:  _memory_space_start
   :end-before:  _memory_space_end
   :language: C++

Axom provides the ``axom::MemorySpace`` enum type to define values indicating
the memory space where data in
``axom::Array`` and ``axom::ArrayView`` lives.

``Dynamic`` allows you to define the location at run time, with some caveats
(see :ref:`Core Containers<core-containers>` for more details and examples).

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Useful Links
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`Umpire Tutorial`_ - First two sections cover ``Allocators`` and ``Resources``.


.. _gpu-kernels-label:

===============
Kernels
===============

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
axom::for_all
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``axom::for_all`` can be found in the file
`axom/core/execution/for_all.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/execution/for_all.hpp>`_.

``axom::for_all`` is a wrapper around `RAJA forall`_, which is used to execute
simple for-loop kernels.

This is used in Axom to execute for-loop style kernels that will be run on a
GPU device, or on both a GPU device and a CPU host. For example::

  template <typename ExecSpace, typename KernelType>
  void axom::for_all(const IndexType& N, KernelType&& kernel)

  template <typename ExecSpace, typename KernelType>
  void axom::for_all(const IndexType& begin, const IndexType& end, KernelType&& kernel)

.. note::

  When Axom is built without RAJA, ``axom::for_all`` becomes a ``for``-loop on
  host (CPU).

%%%%%%%%%%%%%%%%
Portability
%%%%%%%%%%%%%%%%

Adherence to the GPU porting guidelines generally results in code that will compile and run on
multiple backends. However, backends such as CUDA require additional guidelines.

**Do not use ``auto`` lambda parameters** with ``axom::for_all`` or the code will not
compile under nvcc.

Do this:

  .. code-block:: cpp

      axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType index) { /* body */});


Do NOT do this:

  .. code-block:: cpp

      axom::for_all<ExecSpace>(n, AXOM_LAMBDA(auto index) { /* body */});


**Pass ArrayView by value**. Data in Axom is often contained in useful containers such as ``axom::Array``.
Use ``axom::ArrayView`` to access array data from within kernels. Views contain
a pointer to the data and memory shape information and can be constructed within
device code to access data arrays. Views are captured by kernels and passed to
the device code; this copies the pointer and shape information into an object that
the device code can access. If the view was passed into a method where the ``axom::for_all``
function is calling a kernel, be sure to pass the view by value so the compiler
does not capture a host reference to the view, causing the kernel to fail.

Do this:

  .. code-block:: cpp

      template <typename ExecSpace>
      void doSomething(axom::ArrayView<int> dataView)
      {
        axom::for_all<ExecSpace>(dataView.size(), AXOM_LAMBDA(axom::IndexType index)
        {
          /* body uses dataView[index] */
        });
      }

Do NOT do this:

  .. code-block:: cpp

      template <typename ExecSpace>
      void doSomething(axom::ArrayView<int> &dataView)
      {
        axom::for_all<ExecSpace>(dataView.size(), AXOM_LAMBDA(axom::IndexType index)
        {
          /* body uses dataView[index] */
          /* It will crash on GPU devices because the host reference was
             captured rather than the object.
           */
        });
      }

If pass by reference is used, create a new view inside the method and then use that
"device" view in the kernel so the device code uses an object captured by value.


**Do not call for_all from protected or private methods**.
The nvcc compiler requires the method containing the kernel to be publicly accessible.
Conditional compilation can be used to make methods public when they should otherwise
be private or protected.

Do this:

  .. code-block:: cpp

      class Algorithm
      {
      public:
        // stuff
      #if !defined(__CUDACC__)
      private:
      #endif
        void helperMethod()
        {
          axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType index) { /* body */});
        }
      };

Do NOT do this:

  .. code-block:: cpp

      class Algorithm
      {
      public:
        // stuff
      private:
        void helperMethod()
        {
          axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType index) { /* body */});
        }
      };

**Avoid for_all within a lambda**. When calling a kernel via ``axom::for_all`` from a surrounding lambda function, 
consider calling ``axom::for_all`` from a class member method instead. The nvcc compiler will
not compile kernel invocations inside lambda functions. This pattern comes up when an intermediate
function is supplied a lambda that uses ``axom::for_all`` such as when handling many data types.

Do this:

  .. code-block:: cpp

      template <typename ExecSpace>
      class Algorithm
      {
      public:
        void execute(conduit::Node &data)
        {
          // Handle any data type
          Node_to_ArrayView(data, [&](auto dataView)
          {
            // Delegate operation to a template member method
            handleData(dataView);
          });
        }
        template <typename DataView>
        void handleData(DataView dataView)
        {
          // Call the kernel here in the member method
          axom::for_all<ExecSpace>(AXOM_LAMBDA(axom::IndexType) { /* body */ });
        }
      };

Do NOT do this:

  .. code-block:: cpp

      template <typename ExecSpace>
      class Algorithm
      {
      public:
        void execute(conduit::Node &data)
        {
          // Handle any data type
          Node_to_ArrayView(data, [&](auto dataView)
          {
            // nvcc will not compile this
            axom::for_all<ExecSpace>(AXOM_LAMBDA(axom::IndexType) { /* body */ });
          });
        }
      };

**Avoid calling lambdas from kernels.** This can work on some systems and not on others.
For best odds at a portable algorithm, design your kernel so it is "one level deep",
and does not result in calling other functions that are also marked ``AXOM_LAMBDA``.

**Specialize templates outside other classes.** It is necessary to extract
a nested class/struct from the containing class before specializing it.

Do this:

  .. code-block:: cpp

     namespace internal
     {
       template <int NDIMS>
       struct B { static void method() { } };

       template <>
       struct B<2> { static void method() { /* 2D-specific method*/ } };
     }

     template <typename ViewType>
     struct A
     {
       void method() { internal::B<ViewType::dimension()>::method(); }
     };

Do NOT do this:

  .. code-block:: cpp

     template <typename ViewType>
     struct A
     {
       template <int NDIMS>
       struct B { static void method() { } };

       template <>
       struct B<2> { static void method() { /* 2D-specific method*/ } };

       void method() { B<ViewType::dimension()>::method(); }
     };

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RAJA::kernel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``RAJA::kernel`` is used to execute kernels implemented using nested loops.

This is used infrequently, mainly seen only in a `few unit tests`_.

Your general go-to will be ``axom::for_all``.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Useful Links
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`RAJA Loops`_ - Covers ``RAJA::forall``, ``RAJA::kernel``, ``RAJA::launch``
kernel execution methods.


.. _gpu-execution-spaces-label:

================================
Execution Spaces & Policies
================================

Axom's execution spaces can be found in the file
`axom/core/execution/execution_space.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/execution/execution_space.hpp>`_.

Axom's execution spaces are derived from an ``axom::execution_space<ExecSpace>``
traits class containing RAJA execution policies and default Umpire memory
allocators associated with each space.

Axom currently supports `four execution spaces`_, each one a type with the
following specialization of the ``execution_space`` class:

* ``SEQ_EXEC`` - Sequential execution policies on host
* ``OMP_EXEC`` -  OpenMP execution policies on host
* ``CUDA_EXEC`` - CUDA execution policies in Unified Memory (host + device)
* ``HIP_EXEC`` - HIP execution policies in Unified Memory (host + device)

Additionally, ``HIP_EXEC`` and ``CUDA_EXEC`` types are templated by the
number of threads and SYNCHRONOUS or ASYNC execution:

.. literalinclude:: ../../../axom/core/execution/internal/cuda_exec.hpp
   :start-after:  _cuda_exec_start
   :end-before:  _cuda_exec_end
   :language: C++

Each execution space provides:

* Axom policies that are type aliases of RAJA policies to be used with kernels,
  RAJA types, and RAJA operations

  * ``loop_policy`` - For `RAJA scans`_ and other operations; ``axom::for_all``
    uses the loop_policy from the templated execution space.
  * ``reduce_policy`` - For  `RAJA reduction types`_ that perform reduction
    operations:

  .. literalinclude:: ../../../axom/core/examples/core_acceleration.cpp
     :start-after:  _gpu_reduce_start
     :end-before:  _gpu_reduce_end
     :language: C++

  * ``atomic_policy`` - For `RAJA atomic operations`_ that avoid race
    conditions when updating data values:

  .. literalinclude:: ../../../axom/core/examples/core_acceleration.cpp
     :start-after:  _gpu_atomic_start
     :end-before:  _gpu_atomic_end
     :language: C++

  * ``sync_policy`` - For Axom's `synchronize`_ function, which is a wrapper
    around ``RAJA::synchronize()``. Synchronizes execution threads when using
    an asynchronous ``loop_policy``:

  .. literalinclude:: ../../../axom/core/execution/synchronize.hpp
     :start-after:  _gpu_synchronize_start
     :end-before:  _gpu_synchronize_end
     :language: C++

* Umpire allocator defaults

  * ``memory_space`` - The memory space abstraction for use by
    :ref:`Core Containers<core-containers>` like ``axom::Array``.
  * ``allocatorID()`` - Gets the allocator ID for the Umpire resource to use in
    this execution space.

* General information on the execution space

  * ``name()`` - Name of the execution space
  * ``onDevice()`` - Is the execution space on device? (True/False)
  * ``valid()`` - Is the execution space valid? (True)
  * ``async()`` - Is the execution space asynchronous? (True/False)

The :doc:`Core <../../../axom/core/docs/sphinx/index>` component provides a set of nested execution policies
located at
`axom/axom/execution/nested_for_exec.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/execution/nested_for_exec.hpp>`_
to be used with
``RAJA::kernel`` e.g. for iterating over mint meshes.  (These generic policies formerly resided in
:doc:`Mint <../../../axom/mint/docs/sphinx/index>` and have been moved to
:doc:`Core <../../../axom/core/docs/sphinx/index>`.)

.. note::

  When Axom is built without RAJA, only ``SEQ_EXEC`` is available
  for host (CPU) execution. When Axom is built with RAJA but without
  Umpire for memory management on device, only
  ``SEQ_EXEC`` and ``OMP_EXEC`` is available for host (CPU) execution.


.. _gpu-tips-label:

================================
General, Rough Porting Tips
================================

* Start with figuring out what memory you need on device, and use
  ``axom::Array``, ``axom::ArrayView``, and
  :ref:`memory_managment routines<gpu-memory-label>`
  to do the allocations:

  .. code-block:: c

    // Allocate 100 2D Triangles in unified memory
    using cuda_exec = axom::CUDA_EXEC<256>;
    using TriangleType = axom::primal::Triangle<double, 2>;
    axom::Array<Triangle> tris (100, axom::execution_space<cuda_exec>::allocatorID()));
    axom::ArrayView<Triangle> tris_view(tris);

    // Allocate the sum of Triangle areas
    using reduce_pol = typename axom::execution_space<cuda_exec>::reduce_policy;
    RAJA::ReduceSum<reduce_pol, double> totalArea(0);

* Using an ``axom::for_all`` kernel with a device policy, attempt to
  access and/or manipulate the memory on device:

  .. code-block:: c

      axom::for_all<cuda_exec>(
      100,
      AXOM_LAMBDA(int idx) {
        // Set values on device
        tris_view[idx] = Triangle();
        totalArea = 0;
      });

* Add the functions you want to call on device to the ``axom::for_all`` kernel:

  .. code-block:: c

      axom::for_all<cuda_exec>(
      100,
      AXOM_LAMBDA(int idx) {
        tris_view[idx] = Triangle();
        totalArea = 0;

        // Call area() method on device
        double area = tris_view[idx].area();
      });

* Apply a ``__host__ __device__`` annotation to your functions if you see the
  following error or similar::

    error: reference to __host__ function 'area' in __host__ __device__ function

  * Recompiling will likely introduce complaints about more
    functions (the functions being the non-decorated functions your
    newly-decorated functions are calling)::

      error: reference to __host__ function 'abs<double>' in __host__ __device__ function

      error: reference to __host__ function 'signedArea<2>' in __host__ __device__ function

    Keep decorating until all the complaints are gone.

  * Most of the C++ standard library is not available on device. Your options
    are Axom's equivalent functions/classes if it exists, or to add your
    own or rewrite the code to not use standard library.

  * It may not be possible to remove all such warnings on some platforms that support
    both CPU/GPU backends since AXOM_LAMBDA will expand to ``__host__ __device__`` and then
    the compiler will issue warnings about host functions such as RAJA::ReduceSum::~ReduceSum
    for the OpenMP backend being called from ``__host__ __device__`` code. This warning
    can be ignored.

* With no more decorating complaints from the compiler, write the logically
  correct kernel:

  .. code-block:: c

      // Computes the total area of a 100 triangles
      axom::for_all<cuda_exec>(
        100,
        AXOM_LAMBDA(int idx) {
          totalArea += tris_view[idx].area();
      });

  * If at this point your kernel is not working/segfaulting, it is hopefully a
    logical error, and you can debug the kernel without diving into debugging
    tools.
  * Utilize ``printf()`` for debugging output
  * Try using the ``SEQ_EXEC`` execution space

* If your kernel compiles/links but fails at runtime, the cause is often:

  * The data arrays referenced in the kernel are not in the correct memory space for the kernel.
  * A reference to host memory has been captured for use in the kernel.
  * There is nothing wrong with your code and the tooling has failed. In this case, try
    separating the ``axom::for_all`` and your kernel into a separate method or function.
    This can limit the available objects that the compile will attempt to capture and
    increase the likelihood that the kernel will run correctly.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Useful Links
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* List of debugging tools:

  * `Totalview <https://help.totalview.io/>`_ (CUDA)
  * `Nvidia Nsight Developer Tools <https://developer.nvidia.com/tools-overview>`_ (CUDA)

    * `LLNL Guide for Sierra <https://lc.llnl.gov/confluence/display/SIERRA/NVIDIA+Nsight+Developer+Tools>`_

    * ORNL Guides for `Nsight Compute <https://www.olcf.ornl.gov/wp-content/uploads/2020/02/OLCF-Webinar-Nsight-Compute.pdf>`_
      and `Nsight Systems <https://www.olcf.ornl.gov/wp-content/uploads/2020/02/Summit-Nsight-Systems-Introduction.pdf>`_

  * `HPCToolkit <https://github.com/HPCToolkit/hpctoolkit>`_ (CUDA, HIP)
  * `ROCprof <https://rocmdocs.amd.com/en/latest/ROCm_Tools/ROCm-Tools.html>`_ (HIP)
  * `ROCgdb <https://rocmdocs.amd.com/en/latest/ROCm_Tools/ROCgdb.html>`_ (HIP)
* `ORNL Training Archives <https://docs.olcf.ornl.gov/training/training_archive.html>`_
* `LLNL HPC Tutorials <https://hpc.llnl.gov/documentation/tutorials>`_



.. _memory resource type: https://github.com/LLNL/Umpire/blob/develop/src/umpire/resource/MemoryResourceTypes.hpp#L63
.. _Umpire Tutorial: https://umpire.readthedocs.io/en/develop/sphinx/tutorial.html
.. _RAJA forall: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/loop_basic.html#simple-loops-raja-forall
.. _few unit tests: https://github.com/LLNL/axom/search?q=RAJA%3A%3Akernel
.. _RAJA Loops: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/loop_basic.html#elements-of-loop-execution
.. _four execution spaces: https://github.com/LLNL/axom/tree/develop/src/axom/core/execution/internal
.. _RAJA scans: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/scan.html#scans
.. _RAJA reduction types: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/reduction.html#reduction-operations
.. _RAJA atomic operations: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/atomic.html#atomics
.. _synchronize: https://github.com/LLNL/axom/blob/482749b322057ca7ef9a6d786e948692d9f72246/src/axom/core/execution/synchronize.hpp
