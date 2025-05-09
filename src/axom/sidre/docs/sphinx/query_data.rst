.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _datastore-queries-label:

===========================
Datastore Queries
===========================

Sidre provides a couple of methods that can be called to query certain
information about the contents of a datastore object. Each method takes
a Conduit Node and inserts information related to the query into it. The
caller can then inspect the information by retrieving values of specifically 
named fields from the Node. This is often helpful for debugging.

-------------
Buffer Query
-------------

A datastore object can be queried for information about the buffer objects it
contains. For example::

  Datastore* ds = ...;

  conduit::Node n;
  ds->getBufferInfo(n);

  // Print Conduit node contents to stdout, if desired.
  // Note that when running with multiple MPI ranks, each rank will 
  // independently print its node contents to stdout, which may not be wanted.
  n.print();

This method call inserts four fields into the Node that have numeric values
accessible as type ``axom::IndexType``. For example::

  using IndexType = axom::IndexType;

  IndexType num_buffers = n["num_buffers"].value();
  IndexType num_buffers_referenced = n["num_buffers_referenced"].value();
  IndexType num_buffers_detached = n["num_buffers_detached"].value();
  IndexType num_bytes_allocated = n["num_bytes_allocated"].value();

The variables have the following values:

  * ``num_buffers`` : Total number of buffer objects in the datastore.
  * ``num_buffers_referenced`` : Number of buffers referenced (i.e. attached to) a view.
  * ``num_buffers_detached`` : Number of buffers not referenced by a view. Note: ``num_buffers_detached`` = ``num_buffers`` - ``num_buffers_referenced``.
  * ``num_bytes_allocated`` : Total number of bytes allocated in all buffers.

--------------------
Group Subtree Query
--------------------

A group can be queried for information about the data associated with it, or
associated with the entire subtree rooted at the group. For example::

  Group* group = ...;

  bool recursive = false; 

  // get information about a single group
  conduit::Node n;
  group->getDataInfo(n, recursive);

  // Print Conduit node contents to stdout, if desired.
  // Note that when running with multiple MPI ranks, each rank will 
  // independently print its node contents to stdout, which may not be wanted.
  n.print();

  // get information about entire subtree rooted at group
  recursive = true;
  conduit::Node n1;
  group->getDataInfo(n1, recursive);
  n1.print();   // print Conduit node contents to std out, if desired

Similar to the ``Datastore::getBufferInfo`` method described above, the 
``Group::getDataInfo`` method inserts fields into the given Conduit Node
that have numeric values accessible as type ``axom::IndexType``.  For example::

  using IndexType = axom::IndexType;

  IndexType num_groups = n["num_groups"].value();
  IndexType num_views = n["num_views"].value();
  IndexType num_views_empty = n["num_views_empty"].value();
  IndexType num_views_buffer = n["num_views_buffer"].value();
  IndexType num_views_external = n["num_views_external"].value();
  IndexType num_views_scalar = n["num_views_scalar"].value();
  IndexType num_views_string = n["num_views_string"].value();
  IndexType num_bytes_assoc_with_views = n["num_bytes_assoc_with_views"].value();
  IndexType num_bytes_external = n["num_bytes_external"].value();
  IndexType num_bytes_in_buffers = n["num_bytes_in_buffers"].value();

The variables have the following values describing the single group or
entire group subtree:

  * ``num_groups`` : Total number of groups.
  * ``num_views`` : Total number of views.
  * ``num_views_empty`` : Number of views with no associated buffer or data (may or may not be described).
  * ``num_views_buffer`` : Number of views associated with a buffer.
  * ``num_views_external`` : Number of views associated with external data.
  * ``num_views_scalar`` : Number of views associated with a single scalar data item.
  * ``num_views_string`` : Number of views associated with string data.
  * ``num_bytes_assoc_with_views`` : Total number of bytes associated with views (buffer, string, and scalar). Note that this may be an over-count if two or more views share a buffer and their data overlap, for example.
  * ``num_bytes_external`` : Total number of bytes described by external views. Note that this may be an over-count. For example, if two or more views share an external allocation and their data overlaps more bytes may be reported than exist.
  * ``num_bytes_in_buffers`` : Total number of bytes allocated in buffers that are attached to views. Each buffer is counted exactly once if attached to multiple views.

.. important:: ``num_bytes_assoc_with_views`` and ``num_bytes_external`` may over-count the actual data that exists. For example, if two or more views share a buffer or if two or more views share an external allocation and the view data overlaps in either case.
