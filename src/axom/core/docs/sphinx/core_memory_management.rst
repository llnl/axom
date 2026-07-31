.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core memory management
******************************************************

The Axom Core component provides mechanisms to control which memory spaces
are used for allocations to support code execution on CPUs and GPUs and integration
with tools like `Umpire`_. Memory spaces are specified in the Core interface
using enum values:

.. literalinclude:: ../../memory_management.hpp
   :start-after: _memory_space_start
   :end-before: _memory_space_end
   :language: C++

For CPU-accessible memory, Axom uses the concepts of a **global allocator** and
a **host allocator**. Referring to the enum class values above, the Axom
global default allocator is the default for ``MemorySpace::Dynamic`` and for
interface routines that do not specify an allocator; for example, many
``axom::Array`` and ``axom::ArrayView`` APIs.

The Axom host allocator is a process-wide *default* used by legacy convenience
paths that resolve ``MemorySpace::Host`` through global state. New and updated
APIs prefer that host allocation intent be expressed explicitly via the
``axom::HostAllocator`` wrapper type, rather than by relying on the process
default.

Axom must be configured with Umpire enabled to have access to GPU memory
resources. Whether or not Axom is configured with Umpire also controls
default behavior for CPU memory allocations. Specifically, when Umpire is
enabled:

  * the default host allocator is ``malloc`` (i.e., ``MemorySpace::Malloc``)
  * the default global allocator is Umpire's default allocator

When Umpire is disabled:

  * the default host allocator and the default global allocator both refer to
    ``MemorySpace::Malloc``

.. note:: The Axom default host allocator is ``axom::MemorySpace::Malloc`` regardless of whether Axom is configured with Umpire.

The separation of host and global allocation in Axom is intentional because,
in Umpire-enabled builds, the global default may be device, unified, or some
other Umpire allocator, and the host allocator can still be Axom's malloc or
Umpire's Host allocator.

Explicit host allocation
^^^^^^^^^^^^^^^^^^^^^^^^

When an API argument semantically means "host-resident storage" (including host
staging/scratch for device operations), prefer passing an ``axom::HostAllocator``
object. This makes host allocation behavior independent of process-global defaults and
reduces sensitivity to initialization order.

The following illustrates the general pattern::

  axom::HostAllocator hostAlloc {axom::MALLOC_ALLOCATOR_ID};
  // Pass hostAlloc into APIs that allocate host memory or host staging.

For APIs that need both an execution-space allocator and host staging, keep the
two choices separate::

  int dataAllocId = axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID();
  axom::HostAllocator hostAlloc {axom::MALLOC_ALLOCATOR_ID};
  // Pass dataAllocId for primary storage and hostAlloc for host scratch/staging.

Several Core APIs now expose explicit host allocator overloads. For example,
``axom::Array`` constructors can receive both the primary allocator ID and an
``axom::HostAllocator``, and ``axom::fill()`` can receive a host allocator for
temporary host scratch used when filling non-host allocations::

  axom::fill(devicePtr, numValues, value, hostAlloc);

Overloads that do not receive ``axom::HostAllocator`` remain available for
compatibility. They forward through Axom's current default host allocator and
should be treated as legacy convenience paths in new production code.

Changing the default host and global allocators (legacy)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The contents of this section refer to Axom *legacy* convenience routines,
meaning those that do not take an explicit host allocator argument.
The default host allocator controls what Axom uses when code asks for
``MemorySpace::Host`` through these legacy convenience paths. This is separate from
the global default allocator used for ``MemorySpace::Dynamic``.

In practice, that means:

* ``axom::setDefaultHostAllocator(...)`` changes where legacy
  ``MemorySpace::Host`` allocations go.
* ``axom::setDefaultAllocator(...)`` changes the global default allocator for
  ``MemorySpace::Dynamic``.
* Changing one does not automatically change the other.

Axom host allocator choice can still be selected at run time for compatibility
paths. For example, to switch between Axom's malloc-backed host allocator and
the platform host allocator::

  // set Axom host allocator to Axom malloc
  axom::setDefaultHostAllocator(axom::MemorySpace::Malloc);

  // set Axom host allocator to Umpire Host allocator
  axom::setDefaultHostAllocator(axom::MemorySpace::Host);

You can also inspect the current legacy selection::

  int hostAllocId = axom::getDefaultHostAllocatorID();
  int hostSpaceId = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host);

  // hostAllocId and hostSpaceId refer to the same allocator.

If you need to preserve legacy ``MemorySpace::Host`` behavior in an
Umpire-enabled build, you can resolve the host resource allocator ID yourself
and install it as Axom's default host allocator::

  // set Axom host allocator to an explicitly chosen Umpire allocator
  int hostId =
    axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
  axom::setDefaultHostAllocator(hostId);

The following example shows that the host allocator and global allocator are
independent. Here, ``MemorySpace::Dynamic`` is set to Umpire Host while
``MemorySpace::Host`` still uses Axom malloc::

  axom::setDefaultHostAllocator(axom::MemorySpace::Malloc);
  axom::setDefaultAllocator(axom::MemorySpace::Host);

  int dynamicId = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Dynamic);
  int hostId = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host);

  // dynamicId is the Umpire Host allocator
  // hostId is axom::MALLOC_ALLOCATOR_ID

Once Axom has made an allocation from the current host allocator, that host
allocator selection is fixed for the remainder of the process. Configure the
desired host allocator before creating ``MemorySpace::Host`` allocations or
other allocations that use the current Axom host allocator. For example::

  axom::setDefaultHostAllocator(axom::MemorySpace::Host);

  int* values =
    axom::allocate<int>(16, axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host));

  // After this allocation, changing the default host allocator will abort.
  // For example, axom::setDefaultHostAllocator(axom::MemorySpace::Malloc);
  // will abort when Axom is configured with Umpire since, in that case,
  // axom::MemorySpace::Host is the Umpire Host allocator, which is different than
  // Axom's malloc allocator.

Similarly, the Axom global allocator can be changed. For example::

  // set Axom global allocator to Umpire Host
  axom::setDefaultAllocator(axom::MemorySpace::Host);

  // set Axom global allocator to Umpire Unified memory allocator
  axom::setDefaultAllocator(axom::MemorySpace::Unified);

  // set Axom global allocator to an explicitly chosen Umpire allocator
  int allocId =
    axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Pinned);
  axom::setDefaultAllocator(allocId);

One important thing to note is that, when Umpire is enabled, you cannot set
the global allocator to Axom malloc::

  // currently, can't do this -- code will abort
  axom::setDefaultAllocator(axom::MALLOC_ALLOCATOR_ID);

.. _Umpire: https://umpire.readthedocs.io
