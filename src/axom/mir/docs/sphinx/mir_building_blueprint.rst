.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

**************************
Building Blueprint Output
**************************

While views are provided to help write algorithms that process existing Blueprint data,
the views are primarily meant to be read-only. Blueprint output from algorithms is built
in the usual Conduit way by adding new key:value data into paths within the Conduit node
hierarchy. Axom provides a ``ConduitAllocateThroughAxom`` object to help allocate Conduit's
bulk data through Axom's memory allocation routines. Since data processing algorithms in
Axom are often templated on an execution space, the ConduitAllocateThroughAxom object is
also templated on execution space, which enables it to install various allocators for
Conduit. To make Conduit allocate data through Axom, set a Conduit node's allocator
to the allocator returned by ``ConduitAllocateThroughAxom::getConduitAllocatorID()``.
After setting the allocator, the node's memory can be allocated by calling the ``Node::set()``
method and passing a ``conduit::DataType`` object that encodes the data type and size.
After allocating data in the Conduit node, wrap the data in an ``axom::ArrayView`` using
the ``make_array_view`` function. Once Conduit nodes have a corresponding ArrayView, the
ArrayViews can be passed to kernels to fill in the memory of the Conduit node. Remember
to pass ArrayViews and not the Conduit nodes to device kernels.

.. code-block:: cpp
    template <typename ExecSpace>
    void makeNewField(conduit::Node &n_mesh)
    {
      namespace bputils = axom::mir::utilities::blueprint;

      // This object registers Axom's allocation functions with Conduit.
      bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

      // Make the new field normally by adding members to the n_mesh Conduit node.
      conduit::Node &n_field = n_mesh["fields/newField"];
      n_field["topology"] = "mesh";
      n_field["association"] = "element";     
      conduit::Node &n_values = n_field["values"];

      // Set the allocator so Axom will be used to allocate the node's memory.
      // This is key when working on GPU platforms.
      n_values.set_allocator(c2a.getConduitAllocatorID());

      // Allocate memory in the right memory space for ExecSpace.
      // The cpp2conduit template gets the Conduit data type id for supported C++ types.
      const int nzones = 100;
      n_values.set(conduit::DataType(bputils::cpp2conduit<double>::id, nzones));

      // Wrap the node in an ArrayView.
      auto values = bputils::make_array_view<double>(n_values);

      // Fill in the values in a kernel through the "values" ArrayView.
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          values[index] = sin(M_PI * (index / 100.));
        });
    }

*********************
Common Memory Issues
*********************

If your program terminates while executing a device kernel like the one above, and the problem
is not accessing memory out of bounds, then the most likely problem is the memory was not
allocated in the right memory space. The most common issue is attempting to access host
memory in the device kernel. Check whether the Conduit node whose memory backs the
``ArrayView`` set the right allocator.

If your program terminates while printing the Conduit node or saving it to disk, the node's
memory was likely allocated on the device and the memory needs to be moved to the host.

.. code-block:: cpp

    // Assume we have some algorithm output in "n_device".
    // Move that memory to the host.
    conduit::Node n_host;
    axom::mir::utilities::blueprint::copy<axom::SEQ_EXEC>(n_host, n_device);
    // Now, access the memory on the host through n_host.
    n_host.print();
