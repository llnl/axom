.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _buffer-label:

==========
Buffer
==========

A Sidre ``Buffer`` object holds a contiguous array of data *described* by a
data type and length. The data owned by a buffer is unique to that buffer;
i.e., buffers do not share data.

A buffer can be created without a data description and then described
later in a separate operation, or it can be described when it is created.
In either case, data description and allocation are distinct operations. This
allows an application to create buffers it needs, then assess the types and
amount of data they will hold before deciding how and when to allocate data.

.. note:: * ``Buffer`` objects can only be created and destroyed using
            ``DataStore`` class methods. Buffer constructor and destructor
            methods are private (see :ref:`datastore-label`).
          * Each ``Buffer`` object has a unique integer identifier generated
            when the buffer is created. If you want to retrieve a buffer or
            interact with it directly, you must keep a pointer to it or note
            its id so that you can retrieve it from a ``DataStore`` when needed.

``Buffer`` objects hold data for Sidre ``View`` objects in most cases. Each
buffer maintains references to the views that refer to its data. These
references are created when a buffer object is attached to a view, or when
data is allocated with a view. Data stored in a buffer may be accessed through
a view or through the buffer. See :ref:`view-label` for more information about
``View`` objects.

The ``Buffer`` class provides the following operations:

 * Retrieve the unique id of the buffer
 * Query whether a buffer is *described* or *allocated*
 * Describe buffer data (type and number of elements)
 * Allocate, reallocate, deallocate buffer data
 * Copy a given number of bytes of data from a given pointer to a buffer
   allocation
 * Get data held by a buffer as a pointer or ``conduit::Node::Value`` type
 * Query the data owned by a buffer: type, number of elements, total number
   of bytes, number of bytes per element, etc.
 * Retrieve the number of views the buffer is attached to
 * Copy buffer description and its data to/from a ``conduit::Node``

One can iterate through the buffers in the datastore using either "range-for" or iterator syntax:

.. code-block:: C++

    // 'range-for' syntax
    for(auto& buff: datastore->buffers()) { /* ... */ }

    // 'iterator' syntax:
    for(auto it = datastore->buffers().begin(),
          itEnd = datastore->buffers().end(); it != itEnd; ++it)
    {
      auto& buff = *it;
      /* ... */
    }
