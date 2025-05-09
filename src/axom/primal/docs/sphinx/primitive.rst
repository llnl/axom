Primitives
----------

Primal includes the following primitives:

- Point
- Segment, Ray, Vector
- Plane, Triangle, Polygon
- Quadrilateral
- Sphere
- Tetrahedron
- Hexahedron
- BoundingBox, OrientedBoundingBox
- Polyhedron

.. note:: Primitives in Axom use a right-handed coordinate system.

Classes in Primal are templated on coordinate type (double, float, etc.) and
dimension.  The primitives do not inherit from a common base class.  This was a
design choice in favor of simplicity and performance.  Geometric primitives can
be tested for equality and can be printed to strings.

Primal also includes functions to merge a pair of BoundingBox or a pair of 
OrientedBoundingBox objects and to create new OrientedBoundingBox objects 
from a list of points.

The following includes header files for primal's primitives
as well as some ``using`` directives and ``typedef`` statements that will be used in
the examples. Header files for operations will be shown next to code examples.
Although the examples ``#include`` separate class header files, it is easier and
less error-prone to write ``#include axom/primal.hpp``.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _prims_header_start
   :end-before: _prims_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _using_start
   :end-before: _using_end
   :language: C++

