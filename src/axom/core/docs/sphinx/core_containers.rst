.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _core-containers:

******************************************************
Core containers
******************************************************

.. warning::
   Axom's containers are currently the target of a refactoring effort. The 
   following information is subject to change, as are the interfaces for the 
   containers themselves. When changes are made to the Axom containers, the
   changes will be reflected here.

Axom Core contains the ``Array``, ``ArrayView``, and ``StackArray`` classes.  
Among other things, these data containers facilitate porting code that uses 
``std::vector`` to GPUs.

#####
Array
#####

``Array`` is a multidimensional contiguous container template. In the 
1-dimensional case, this class behaves similar to ``std::vector``. In higher 
dimensions, some vector-like functionality, such as ``push_back``, are not
available and multidimensional-specific operations will mirror the 
``ndarray`` provided by ``numpy`` when possible. 

The ``Array`` object manages all memory. Typically, the ``Array`` object will
allocate extra space to facilitate the insertion of new elements and minimize the number of reallocations.
The actual capacity of the array (i.e., total number of elements that
the ``Array`` can hold) can be queried via the ``capacity()`` method.
When allocated memory is used up, inserting a new element triggers a
reallocation.  At each reallocation, extra space is allocated
according to the **resize_ratio** parameter, which is set to 2.0
by default. To return all extra memory, an application can call
`shrink()`.

.. warning:: Reallocations tend to be costly operations in terms of performance.
  Use ``reserve()`` when the number of nodes is known a priori, or
  use a constructor that takes an actual size and capacity when possible.
  
.. note:: The ``Array`` destructor deallocates and returns all memory associated
  with it to the system.

Here's an example showing how to use ``Array`` instead of ``std::vector``.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraybasic_start
   :end-before: _arraybasic_end
   :language: C++

The output of this example is::

  Length of a = 3
  After appending a value, a's length = 4
  Array a = [2, 5, 11, 4]
  After inserting two values, Array a = [2, 5, 6, 11, 1, 4]

Applications commonly store tuples of data in a flat array or a ``std::vector``.
In this sense, it can be thought of as a two-dimensional array.  ``Array``
supports arbitrary dimensionalities but an alias, ``MCArray``
(or Multi-Component Array), is provided for this two-dimensional case.

The ``MCArray`` (i.e., ``Array<T, 2>``) class formalizes tuple storage, as shown in the next example.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraytuple_start
   :end-before: _arraytuple_end
   :language: C++

The output of this example is::

  MCArray b with 2 3-tuples = [
    [1, 4, 2]
    [8, 0, -1]
  ]
  MCArray b with 4 3-tuples = [
    [1, 4, 2]
    [0, -1, 1]
    [1, -1, 0]
    [8, 0, -1]
  ]

#########
ArrayView
#########

It is also often useful to wrap an external, user-supplied buffer without taking ownership
of the data.  For this purpose Axom provides the ``ArrayView`` class, which is a lightweight
wrapper over a buffer that provides one- or multi-dimensional indexing/reshaping semantics.
For example, it might be useful to reinterpret a flat (one-dimensional) array as a two-dimensional array.
This is accomplished via ``MCArrayView`` which, similar to the ``MCArray`` alias, is an alias for ``ArrayView<T, 2>``.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _extbuffer_start
   :end-before: _extbuffer_end
   :language: C++

The output of this example is::

  Array a = [2, 5, 6, 11, 1, 4]
  MCArrayView c with 3 2-tuples = [
    [2, 5]
    [6, 11]
    [1, 4]
  ]
  Array a and MCArrayView c use the same memory, a's internal buffer.
  Array a = [1, 5, 6, 9, 1, 4]
  MCArrayView c with 3 2-tuples = [
    [1, 5]
    [6, 9]
    [1, 4]
  ]

.. note:: The set of permissible operations on an ``ArrayView`` is somewhat limited,
  as operations that would cause the buffer to resize are not permitted.

``ArrayView`` can view memory spaces at regular intervals while
ignoring spaces in between.  By default, the spacing is one, meaning
adjacent elements are one space apart.  A spacing of 2 sees every
other elements in memory.  A spacing of ``N`` sees every ``N``-th
element.  Spacing is set in the ``ArrayView`` constructor.

The following example creates a 2D array of 4-tuples and uses an
``ArrayView`` to access the 3rd item in each 4-tuple.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _spacing_start
   :end-before: _spacing_end
   :language: C++

The output of this example is::

  Third components of 2D array of 4-tuples:
  a(0,0) = (0,0).2
  a(0,1) = (0,1).2
  a(0,2) = (0,2).2
  a(1,0) = (1,0).2
  a(1,1) = (1,1).2
  a(1,2) = (1,2).2

In the future, it will also be possible to restride an ``ArrayView``.

Iteration is also identical between the ``Array`` and ``ArrayView`` classes.
In particular:

  * ``operator()`` indexes into multidimensional data. Currently, the number of indexes 
    passed must match the dimensionality of the array.
  * ``operator[]`` indexes into the full buffer, i.e., ``arr[i]`` is equivalent to ``arr.data()[i]``.
  * ``begin()`` and ``end()`` refer to the full buffer.

Consider the following example:

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _iteration_start
   :end-before: _iteration_end
   :language: C++

The output of this example is::

  In ArrayView c, index (0, 0) yields 1
  In ArrayView c, index (0, 1) yields 5
  In ArrayView c, index (1, 0) yields 6
  In ArrayView c, index (1, 1) yields 9
  In ArrayView c, index (2, 0) yields 1
  In ArrayView c, index (2, 1) yields 4
  Range-based for loop over ArrayView c yields: 1 5 6 9 1 4
  Standard for loop over ArrayView c yields: 1 5 6 9 1 4

########################
Using Arrays in GPU Code
########################

Instead of writing kernels and device functions that operate on raw pointers, we can use ``ArrayView``
in device code. The basic "workflow" for this process is as follows:

  1. Create an ``Array`` allocated in device-accessible memory via either specifying an allocator ID
     or using a class template parameter for the desired memory space.
  2. Write a kernel that accepts an ``ArrayView`` parameter **by value**, not by reference or pointer.
  3. Create an ``ArrayView`` from the ``Array`` to call the function.  For non-templated kernels
     an implicit conversion is provided.


The full template signature for ``Array`` (``ArrayView`` has an analogous signature) is
``Array<typename T, int DIM = 1, MemorySpace SPACE = MemorySpace::Dynamic>``.  Of particular interest
is the last parameter, which specifies the memory space in which the array's data are allocated.
The default, ``Dynamic``, means that the memory space is set via an allocator ID at runtime.

.. note:: Allocating ``Array`` s in different memory spaces is only possible when Umpire is available.
   To learn more about Umpire, see the `Umpire documentation <https://readthedocs.org/projects/umpire>`_

Setting the ``MemorySpace`` to an option other than ``Dynamic`` (for example, ``MemorySpace::Device``) provides
a compile-time guarantee that data can always be accessed from a GPU.  "Locking down" the memory space at
compile time can help to prevent illegal memory accesses and segmentation faults when pointers are dereferenced
from the wrong execution space.

To summarize, there are a couple different options for creating an ``ArrayView``.
Consider a function that takes as an argument an ``ArrayView`` on the device:

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _basic_array_function_start
   :end-before: _basic_array_function_end
   :language: C++

To create an argument to this function we can select the space either at runtime or at compile-time as follows:

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _basic_array_device_create_start
   :end-before: _basic_array_device_create_end
   :language: C++

The first way we can create the required ``ArrayView`` is by implicit conversion, which also simplifies
the process of "locking down" a ``MemorySpace::Dynamic`` array to an explicit memory space - ``MemorySpace:Device`` in this case.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _basic_array_device_implicit_start
   :end-before: _basic_array_device_implicit_end
   :language: C++

.. warning:: If we had attempted to convert from a ``MemorySpace::Dynamic`` array that had been allocated in host memory,
  for example, an error would be produced at runtime.

We can also explicitly construct the ``ArrayView`` before calling the function.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _basic_array_device_explicit_start
   :end-before: _basic_array_device_explicit_end
   :language: C++

A more realistic example of this functionality involves a GPU kernel requiring
that its argument arrays be allocated in a specific memory space.
To illustrate how different memory spaces can be required, the following kernel requires that its
input arrays ``A`` and ``B`` are in unified memory and its output array ``C`` is in device memory.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _device_kernel_start
   :end-before: _device_kernel_end
   :language: C++

The following snippet illustrates how one would create and initialize the inputs/outputs to this kernel.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _device_array_create_start
   :end-before: _device_array_create_end
   :language: C++

.. note:: Unless the Dynamic memory space is in use, the ``Array`` constructor will
   ignore an allocator ID that doesn't match its memory space, and in debug
   builds will print a warning at runtime.

We can now launch the kernel and display the results via a transfer back to host-accessible memory:

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _device_array_call_start
   :end-before: _device_array_call_end
   :language: C++

If RAJA is available, we can also use Axom's acceleration utilities to perform an operation on the GPU
via a lambda:

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _array_w_raja_start
   :end-before: _array_w_raja_end
   :language: C++

By default, ``Array`` copies and moves will propagate the allocator ID; this ensures that objects
with ``Array`` members do not accidentally move their data to the host when copied or moved:

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _device_array_propagate_start
   :end-before: _device_array_propagate_end
   :language: C++

##########
StackArray
##########

The ``StackArray`` class is a work-around for a limitation in older versions
of the nvcc compiler, which do not capture arrays on the stack in device 
lambdas.  More details are in the API documentation and in the tests.
