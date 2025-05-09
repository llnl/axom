.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
An introductory example
******************************************************

As an introduction to the core concepts in Sidre and how they work, here is an
example where we construct the Sidre data hierarchy structure shown in the 
following figure:

.. figure:: figs/sidre_datastore_example.png
   :figwidth: 650px
   :alt: diagram of an example datastore

   ..

   +------------------------------------+-------------+
   | Symbol                             | Sidre class |
   +====================================+=============+
   | .. image:: figs/roundrectangle.png | Group       |
   |    :width: 120px                   |             |
   +------------------------------------+-------------+
   | .. image:: figs/rectangle.png      | View        |
   |    :width: 61px                    |             |
   +------------------------------------+-------------+
   | .. image:: figs/hexagon.png        | Attribute   |
   |    :width: 72px                    |             |
   +------------------------------------+-------------+

The diagram represents a datastore container, holding all Sidre objects and
provides the main interface to access those objects.
Rounded rectangles represent Sidre ``Group`` objects. Each ``Group`` has a 
name and one parent group, except for the root group (i.e. "/"), which 
has no parent. A datastore provides one root group (i.e. "/") that is 
created when a ``DataStore`` object is constructed; thus, an application does 
not create the root. A ``Group`` may have zero or more child groups
(indicated by an arrow from the parent to each child). Each group also owns
zero or more ``View`` objects, which are shown as rectangles. An arrow points 
from a group to each view it owns.

A Sidre ``View`` object has a name and, typically, some data associated with it.
This example shows various types of data that can be described by a view, 
including scalars, strings, and data in arrays (both externally allocated 
and owned by Sidre ``Buffer`` objects). Each array view has a data pointer 
and describes data in terms of data type, number of elements, offset, and 
stride. Data pointers held by array view are shown as dashed arrows.

A ``DataStore`` also contains a collection of ``Buffer`` objects, shown as 
segmented rectangles.

A ``DataStore`` contains a list of ``Attribute`` objects. Each ``Attribute`` 
is outlined with a hexagon and defines a metadata label and a default value 
associated with that label. In the example, the datastore has "vis" (with 
default value 0) and "restart" (with default value 1) attributes.
Default attributes apply to all views unless explicitly set for individual 
views. In the example, the "temp" and "rho" views have the "vis" attribute 
set to 1.

Various aspects of Sidre usage are illustrated in the C++ code shown next. 
Sidre provides full C and Fortran APIs that can also be used to generate 
the same result.

First, we create a ``DataStore`` object, define some ``Attribute`` objects 
along with their default values, and add some child ``Group`` objects to the
root ``Group``.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_creategroups_start
   :end-before: _first_example_creategroups_end
   :language: C++

The ``createViewScalar()`` method lets an application store scalar 
values in a view owned by the group object the method is called on.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_state_start
   :end-before: _first_example_state_end
   :language: C++

This example stores (x, y, z) node position data in one array.
The array is managed by a ``Buffer`` object and three ``View`` objects 
point into it. C++ Sidre operations that create a ``Buffer``, a ``Group``, and 
a ``View``, as shown in the following code, return a pointer to the object 
that is created. This allows chaining operations. (Chaining is supported in 
the C++ API but not in C or Fortran.)

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_nodes_start
   :end-before: _first_example_nodes_end
   :language: C++

The last two integral arguments to the ``createView()`` method specify the
offset from the beginning of the array and the stride of the data. Thus,
the x, y, z values for each position are stored contiguously with the
x values, y values, and z values each offset from each other by a stride of
three in the array. 

The next snippet creates two views ("temp" and "rho") and allocates 
each of their data as an array of type double with length 'eltcount'. Then, 
it sets an attribute ("vis") on each view with a value of 1. Lastly, it 
creates a group ("ext") that has a view holding an external pointer 
("region"). The ``apply()`` method describes the view data as an array of
integer type and length 'eltcount'. Note that it is the responsibility of
the caller to ensure that the allocation to which the "region" pointer 
references is adequate to contain that data description.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_fields_start
   :end-before: _first_example_fields_end
   :language: C++

The next code example shows various methods to retrieve a group object and data 
from a view in the group hierarchy.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_access_start
   :end-before: _first_example_access_end
   :language: C++

In the last section, the code accesses the arrays associated with the views 
"y", "temp", and "region". While "temp" and "region" have the default offset 
(0) and stride (1), "y" has offset 1 and stride 3 (as described earlier). 
The pointer returned by ``View::getPointer()`` always points to the first 
data element described by the view (the view takes care of the offset),
but use of a stride other than 1 must be done by the code itself.

Unix-style path syntax using the slash ("/") delimiter is supported for
traversing Sidre group and view hierarchies and accessing their 
contents.  However, '..' and '.' syntax (up-directory and current directory) 
is not supported. This usage is shown in the last call to ``getView()`` in 
the code example above. The method call retrieves the view named "region" 
in the group named "ext" that is a child of the "fields" group. Character 
sequences before the first slash and between two consecutive slashes are 
assumed to be group names (describing parent-child relationships). For 
this method, and others dealing with view objects, the sequence following 
the last slash is assumed to be the name of a view. Similar path syntax 
can be used to retrieve a group, create a group and a view, a
nd so forth.
