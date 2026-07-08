.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core memory management
******************************************************

The Axom Core component provides mechanisms to control which memory spaces
allocations occur to support code execution on CPUs and GPUs and integration
with tools like `Umpire`_. Memory spaces are specified in the Core interface
using enum values:

.. literalinclude:: ../../memory_management.hpp
   :start-after: _memory_space_start
   :end-before: _memory_space_end
   :language: C++

For CPU-accessible memory, Axom uses the concepts of a **global allocator** and
a **host allocator**. Referring to the enum class values above, the Axom
global default allocator is the default for ``MemorySpace::Dynamic`` and for
interface routines that do not specify an allocator; for example, the
``axom::Array`` and ``axom::ArrayView`` APIs. The Axom host allocator is the
default for ``MemorySpace::Host``, which indicates the allocator Axom should
use when someone explicitly asks for host-allocated CPU memory.

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

Changing the default host and global allocators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default host allocator controls what Axom uses when code asks for
``MemorySpace::Host``. This is separate from the global default allocator used
for ``MemorySpace::Dynamic``.

In practice, that means:

* ``axom::setDefaultHostAllocator(...)`` changes where ``MemorySpace::Host``
  allocations go.
* ``axom::setDefaultAllocator(...)`` changes the global default allocator for
  ``MemorySpace::Dynamic``.
* Changing one does not automatically change the other.

Axom host allocator choice can be selected at run time. For example, to switch
between Axom's malloc-backed host allocator and the platform host allocator::

  // set Axom host allocator to Axom malloc
  axom::setDefaultHostAllocator(axom::MemorySpace::Malloc);

  // set Axom host allocator to Umpire Host allocator
  axom::setDefaultHostAllocator(axom::MemorySpace::Host);

You can also inspect the current selection::

  int hostAllocId = axom::getDefaultHostAllocatorID();
  int hostSpaceId = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host);

  // hostAllocId and hostSpaceId refer to the same allocator.

If you need to be explicit in an Umpire-enabled build, you can resolve the host
resource allocator ID yourself and install it as Axom's default host
allocator::

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
