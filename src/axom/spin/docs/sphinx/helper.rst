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

``Mortonizer`` maps a 2D, 3D, or 4D integer grid point to a single integer
called a Morton index. Sorting points by this index puts them in the order
of the Z-order curve [#f1]_. The ``PointHash`` class uses ``Mortonizer`` to
hash points for ``std::unordered_map``, ``axom::FlatMap`` or similar containers.

The Morton index works with integers, so the ``Mortonizer`` and
its dependent class ``PointHash`` use the  ``RectangularLattice``
to quantize floating points coordinates in the following example.

``Mortonizer`` builds the index by interleaving coordinate bits.
For example, if the binary coordinates are ``x = x1 x0`` and ``y = y1 y0``,
the index is ``y1 x1 y0 x0``, written from the most significant bit 
to the least significant bit. Each pair contains one bit from each coordinate.
The x-coordinate, axis zero, supplies the less significant bit of each pair.

The index type limits how many coordinate bits fit. A 64-bit ``MortonIndexType``
holds 32 bits per coordinate in 2D, 21 in 3D, or 16 in 4D.
Only the least significant bits within this limit contribute to the index.

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

Hilbertizer
^^^^^^^^^^^

``Hilbertizer<CoordType, IndexType, DIM>`` assigns an integer index to each
point on a bounded integer grid in 2D, 3D, or 4D. The index is the point's
position along a Hilbert curve, a path that visits every grid point once.
Consecutive indices correspond to points one grid step apart along one axis.
``toIndex()`` finds a point's index, and ``toPoint()`` recovers its coordinates.

.. code-block:: cpp

   #include "axom/spin/HilbertIndex.hpp"
   #include <cstdint>

   using H = axom::spin::Hilbertizer<std::int64_t, std::uint64_t, 3>;
   const H::PointType point {5, 10, 20};
   const int bits = 5;
   const auto index = H::toIndex(point, bits);  // 7865
   const auto decoded = H::toPoint(index, bits);  // {5, 10, 20}

The depth ``bits`` is the number of bits per coordinate. For example,
``bits = 5`` gives 32 positions along each axis, numbered 0 through 31.
At a valid depth, each integer point in ``[0, 2^bits - 1]^DIM`` maps
to exactly one index in ``[0, 2^(DIM * bits) - 1]``, and vice versa. Use the
same depth to encode and decode. Store the dimension and depth with the key.
Every conversion and check takes ``bits`` explicitly, since a point's index
changes with the depth.

At depth one, indices 0 through 3 visit the 2D points ``(0,0)``, ``(0,1)``,
``(1,1)``, and ``(1,0)``, in that order.

.. rubric:: Depth and valid inputs

Coordinate and index types must be integers other than ``bool``, at most 64
bits wide. Signed types use only nonnegative values. The maximum depth is
``MAX_UNIQUE_BITS = min(COORD_BITS, INDEX_BITS / DIM)``, where the bit counts
exclude sign bits. ``maxBitsPerCoord()`` returns this limit.
Narrow coordinate types can reduce the limits shown below.

.. list-table:: Maximum bits per coordinate
   :header-rows: 1

   * - Index type
     - 2D
     - 3D
     - 4D
   * - ``uint32_t``
     - 16
     - 10
     - 8
   * - ``int32_t``
     - 15
     - 10
     - 7
   * - ``uint64_t``
     - 32
     - 21
     - 16
   * - ``int64_t``
     - 31
     - 21
     - 15

Use ``isValidDepth(bits)``, ``isValidPoint(point, bits)``, and
``isValidIndex(index, bits)`` to check inputs. Depth zero, negative values, and
values beyond the selected depth are invalid. Conversion functions require
valid inputs. Debug assertions check these preconditions. Release builds leave
these checks to the caller. Neither conversion clamps inputs.

The free functions ``convertPointToHilbert<IndexType>(point, bits)`` and
``convertHilbertToPoint<CoordType, DIM>(index, bits)`` call ``toIndex`` and ``toPoint``.
They take the same depth argument and have the same preconditions.

Quantize floating-point positions before encoding, for example with
``RectangularLattice::gridCell()``. The lattice is unbounded. Check that each
cell coordinate is nonnegative and fits the selected depth before encoding.
Decoding returns the integer cell coordinates.

``Hilbertizer`` uses a finite state machine with precomputed tables on the
host and the scalar Skilling algorithm on CUDA and HIP devices. Both produce
the same indices at every supported depth.


.. rubric:: Footnotes

.. [#f1] The Morton index is also known, among other things, as the Z-order curve: 
         see its `Wikipedia page <https://wikipedia.org/wiki/Z-order_curve>`_.
