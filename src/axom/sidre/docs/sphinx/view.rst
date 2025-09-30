.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _view-label:

==========
View
==========

A Sidre ``View`` object describes data and provides access to it. A view can 
describe a (portion of a) data allocation in **any way that is compatible** 
with the allocation. Specifically, the allocation must contain enough bytes 
to support the description. In particular, the data type of a view description 
need not match the types associated with the allocation or types of other views
that refer to that data.

.. note:: ``View`` objects can only be created and destroyed using ``Group``
          class methods. The ``View`` class constructor and destructor are 
          private.

Each view has a name and is owned by one group in a Sidre group hierarchy; i.e.,
its *owning* group. A view maintains a pointer to the group that owns it.

.. note:: * The name (string) of a view **must be unique** within its
            owning group.
          * A view has a unique integer identifier within its owning group, 
            which is generated when the view is created.
          * Views in a group can be accessed by name or integer id.

A view can describe and provide access to data referenced by a 
pointer in one of four ways described below. A view data description includes: 
data type, length (number of elements), offset (from base data pointer) and 
stride (based on the pointer address and data type). 

  * A view can describe (a subset of) data owned by a pre-existing buffer.
    In this case, the buffer is manually *attached* to the view and the
    view's data description is *applied* to the buffer data. Buffer data can be 
    (re)allocated or deallocated by the view if and only if it is the only 
    view attached to the buffer. **In general, a buffer can be attached
    to more than one view.**
  * A view description can be used to allocate data for the view using a 
    ``View`` class ``allocate()`` method similar to buffer data description 
    and allocation (see :ref:`buffer-label`). In this case, the view is 
    usually exclusively associated with a buffer and no other view is allowed 
    to (re)allocate or deallocate the data held by the buffer.
  * A view can **describe** data associated with a pointer to an *external* 
    data allocation. In this case, the view cannot (re)allocate or deallocate 
    the data. However, all other view operations can be applied to the data
    in essentially the same ways as the previous two cases.
  * A view can hold a pointer to an undescribed (*opaque*) data pointer. In 
    this case, the view knows nothing about the type or structure of the data; 
    it can only provide access to it. A user is entirely responsible for 
    casting the pointer to a proper type when accessed through the view, 
    knowing the size of the data, etc.

A view may also refer to a scalar quantity, a tuple (which is like a
scalar but may have multiple values like an array) or a string. Such
views hold their data differently than the pointer cases described
above.

Before we describe the Sidre ``View`` class interface, we present some view 
concepts that describe various *states* a view can be in at any given time. 
Hopefully, this will provide some useful context for the method descriptions 
that follow.

The key view concepts that users should be aware of are: 

  * View data description (data type, number of elements, stride, offset, etc.)
  * View data association (data lives in an attached Sidre buffer,
    accessed via an external pointer, or is a scalar, tuple or string
    owned by the view)
  * Whether the view data description has been applied to the data

The table below summarizes View data associations (rows) and view states with 
respect to the data (columns).

.. figure:: figs/sidre-view-states.png

   This table summarizes Sidre view *data associations* and *data states*. 
   Each row is a data association and each column refers to a data state.
   The True/False entries in the cells indicate return values of the
   ``View`` class methods at the tops of the columns. The circumstances under 
   which those values are returned are noted as well.

The three View data state methods at the tops of columns and their return 
values are:

  * **isDescribed()** returns true if a view has a data description, and
    false otherwise.
  * **isAllocated()** returns true if a view is associated with data, such as
    a non-null pointer, and false otherwise.
  * **isApplied()** returns true if a view has a data description and is
    associated with data that is compatible with that description, and the 
    description has been applied to the data; otherwise false is returned.

The rows indicate data associations; the ``View`` class has methods to query
these as well; e.g., ``isEmpty()``, ``hasBuffer()``, etc. The associations are:

  * **EMPTY.** A view has no associated data; the view may or may not have
    a data description.
  * **BUFFER.** A view has an attached buffer; the view may or may not have 
    a data description and the buffer may or may not be allocated and the
    description (if view has one) may or may not be applied to the buffer data
    (if allocated).
  * **EXTERNAL.** A view has a non-null pointer to external data; the view
    may or may not have a data description and the description, if the view 
    has one, may or may not be applied to the external data.
  * **TUPLE.** A view was created to hold a tuple value; such a view always
    has a valid data description, is allocated, and the description is applied.
    A scalar value is equivalent to a 1-tuple.
  * **STRING.** A view was created to hold a string; such a view always
    has a valid data description, is allocated, and the description is applied.

Note that there are specific consequences that follow from each particular
association/state that a view is in. For example, an EMPTY view cannot have an
attached buffer. Neither can an EXTERNAL, TUPLE or STRING view. A view that
is EMPTY, BUFFER, TUPLE, or STRING cannot be EXTERNAL, etc.

The following lists summarize the ``View`` class interface:

.. note:: Most ``View`` class methods return a pointer to the view object on 
          which the method is called. This allows operations to be chained; 
          e.g., ::

             View* view = ...;
             view->describe(...)->allocate(...)->apply(); 

.. _view-interface-label:

View Property Operations
-----------------------------

 * Retrieve the name or id of the view object.
 * Retrieve the view path name from the root of the tree or the path to the
   group that owns it.
 * Get a pointer to the group that owns the view.
 * Is the view equivalent to another view; i.e., are names and data descriptions
   the same?
 * Rename a view.
 * Clear the view by removing description and data.

Data Association Queries
--------------------------

 * Is view empty?
 * Does view have a buffer attached?
 * Is view associated with external data?
 * Is it a scalar view?
 * Is it a tuple view?
 * Is it a string view?

Data State Queries
-------------------

 * Does view have a data description?
 * Is view data allocated?
 * Is view data description applied to data?
 * Is view opaque; i.e., it has an external pointer and no description?

Data Description Queries
--------------------------

 * Get the type of the data described by a view.
 * Get total number of bytes of data.
 * Get number of elements (total bytes / size of type).
 * Get number of bytes per data element (for type).
 * Get data offet.
 * Get data stride.
 * Get number of dimensions and shape of multi-dimensional data.
 * Get a ``conduit::Schema`` object that contains the view data description.

Data Management Operations
---------------------------

 * Allocate, reallocate, and deallocate view data.
 * Attach buffer to view (with or without data description), 
   and detach buffer from view.
 * Apply current view description to data or apply a new description.
 * Set view scalar value.
 * Set view tuple value.
 * Set view string. 
 * Set external data pointer, with or without a data description. 

Data Access Methods
-----------------------

 * Get a pointer to the view data, actual type or void*.
 * Get scalar value for a scalar view.
 * Retrieve pointer to buffer attached to view.
 * Get a ``conduit::Node`` object that holds the view data.

Attribute Methods
-------------------

 * Query whether a view has an attribute with given id or name.
 * Get attribute associated with a view by id or name.
 * Query whether aAttribute has been set explicitly for view.
 * Reset attribute with given id or name to its default value.
 * Set attribute with given id or name to a given scalar value, tuple or string.
 * Retrieve scalar value, tuple pointer or string of an attribute.
 * Iterate over attributes of a view.

I/O Operations
--------------
 
 * Copy view data description to a ``conduit::Node``.
