.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

***********************
Delaunay Triangulation
***********************

Quest provides a ``quest::Delaunay`` class for incremental construction of a
2D or 3D Delaunay complex from a point set. The class currently builds the
triangulation by inserting points one at a time into an initial bounding mesh,
retriangulating the local cavity after each insertion.

At a high level, the workflow is:

#. define a bounding box that contains the points to be inserted
#. initialize the Delaunay object with that boundary
#. insert points incrementally
#. remove the artificial boundary elements
#. optionally validate the resulting complex and write it to VTK

The current implementation is useful for applications that need a simplicial
mesh over an unstructured point set. Quest also uses this capability in higher
level algorithms such as scattered interpolation.

The example application
``<axom>/src/axom/quest/examples/delaunay_triangulation.cpp`` demonstrates the
basic usage pattern.

The Delaunay class is templated on the dimension:

.. literalinclude:: ../../examples/delaunay_triangulation.cpp
   :start-after: _quest_delaunay_include_start
   :end-before: _quest_delaunay_include_end
   :language: C++

Creating the triangulation
--------------------------

The example creates a bounding box, initializes the Delaunay object, inserts
points, and then removes the artificial boundary:

.. literalinclude:: ../../examples/delaunay_triangulation.cpp
   :start-after: _quest_delaunay_basic_start
   :end-before: _quest_delaunay_basic_end
   :language: C++

The call to ``initializeBoundary()`` is required before inserting points. The
inserted points must lie inside that bounding box.

Validation
----------

Quest provides validation helpers for both the underlying mesh structure and
the Delaunay property itself:

.. literalinclude:: ../../examples/delaunay_triangulation.cpp
   :start-after: _quest_delaunay_validate_start
   :end-before: _quest_delaunay_validate_end
   :language: C++

Output
------

The resulting triangulation can be written to a VTK file for inspection:

.. literalinclude:: ../../examples/delaunay_triangulation.cpp
   :start-after: _quest_delaunay_output_start
   :end-before: _quest_delaunay_output_end
   :language: C++

Current scope
-------------

This page is intentionally conservative. It documents the current user-visible
workflow of Quest's Delaunay triangulation example without trying to freeze the
interface or fully characterize future use cases while the implementation is
still evolving.
