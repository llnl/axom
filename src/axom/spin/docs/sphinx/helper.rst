RectangularLattice
^^^^^^^^^^^^^^^^^^

The ``RectangularLattice`` is a helper class that maps all of N-D space into
a regular, rectangular grid of cells identified by integer coordinates.  The
grid is defined by an origin point and a vector indicating spacing in each dimension.

.. figure:: figs/showRectangularLattice.png
   :figwidth: 300px
   :alt: Diagram showing points and their RectangularLattice bins

The figure shows an example ``RectangularLattice`` in 2D, with its origin
(circled) at (-0.6, -0.2) and spacing (1.2, 0.8).  Given a query point, the
``RectangularLattice`` will return the coordinates of the cell that contains the
point.  It will also return the bounding box of a cell, or the coordinates of a
cell's lower-left corner.

The following example shows the use of the ``RectangularLattice``.  First, include
the header and (if desired) declare type aliases.  Using ``const int in2d = 2``
makes a 2D lattice.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _rectlattice_header_start
   :end-before: _rectlattice_header_end
   :language: C++

Use the ``RectangularLattice`` to find grid cells.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _rectlattice_use_start
   :end-before: _rectlattice_use_end
   :language: C++

Mortonizer
^^^^^^^^^^

The ``Mortonizer`` (along with its associated class ``MortonBase``) implements
the Morton index, an operation that associates each point in N-D space with a
point on a space-filling curve [#f1]_.  The ``PointHash`` class adapts the
``Mortonizer`` to provide a hashing functionality for use with
``std::unordered_map``, ``axom::FlatMap`` or similar container classes.

The math of the Morton index works with integers.  Thus, the ``Mortonizer`` and
its dependent class ``PointHash`` will not work with floating point coordinates.
The following code example shows how the cells of a ``RectangularLattice``,
which have integer coordinates, can be used with a hash table.

The ``Mortonizer`` works by interleaving bits from each coordinate of the input
point and finite computer hardware limits its resolution.  If the
``MortonIndexType`` is 64-bits, then in 2D it can uniquely index the least
significant 32 bits of the x- and y-coordinates.  In 3D, it can uniquely index
the least significant 21 bits of the x-, y-, and z-coordinates.

To use the ``PointHash``, include the header and (as desired) declare type aliases.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _morton_header_start
   :end-before: _morton_header_end
   :language: C++

The ``RectangularLattice`` grid cell associated with a query point can be
stored, using a ``PointHash``, in a ``std::unordered_map``.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _morton_use_start
   :end-before: _morton_use_end
   :language: C++

.. rubric:: Footnotes

.. [#f1] The Morton index is also known, among other things, as the Z-order curve: 
         see its `Wikipedia page <https://wikipedia.org/wiki/Z-order_curve>`_.
