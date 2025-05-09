.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _datastore-label:

==========
DataStore
==========

A Sidre ``DataStore`` object provides the main access point for Sidre contents,
including the data managed by Sidre. A datastore maintains the ``Group`` object
at the root of a Sidre data hierarchy, a collection of ``Buffer`` objects, and 
a collection of ``Attribute`` objects. Generally, the first thing a Sidre user 
does is create the datastore object; this operation also creates the root group.
Apart from providing access to the root group, a datastore provides methods to 
retrieve and interact with buffers and attributes.

The ``DataStore`` class provides methods to retrieve error state and messages
arising from I/O errors reading or writing ``Group`` data:

 * Query or set Conduit error flag
 * Query, append, or clear exception messages from a Conduit I/O error

.. note:: ``Buffer`` and ``Attribute`` objects can only be created and 
          destroyed using ``DataStore`` methods noted below. The ``Buffer`` and
          ``Attribute`` class constructors and destructors are private.

The ``DataStore`` class provides the following methods to manage and interact 
with ``Buffer`` objects:

 * Create, destroy, and allocate data in buffers
 * Query the number of buffers that exist
 * Query whether a buffer exists with a given id
 * Retrieve the buffer with a given id
 * Iterate over the set of buffers in a datastore

Please see :ref:`buffer-label` for more information about using ``Buffer`` 
objects.

The ``DataStore`` class provides the following methods to manage and interact 
with ``Attribute`` objects:

 * Create and destroy attributes
 * Query the number of attributes that exist
 * Query whether an attribute exists with a given name or id
 * Retrieve attribute with a given name or id
 * Iterate over the set of attributes in a datastore

Please see :ref:`attribute-label` for more information about using ``Attribute`` objects.
