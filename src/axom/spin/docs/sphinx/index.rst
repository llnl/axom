.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Spin User Guide
===============

Axom's Spin component provides several data structures to accelerate spatial queries.

- Mortonizer and Hilbertizer map integer points to positions on one-dimensional space-filling curves.
- RectangularLattice computes integer bin coordinates for spatial points.
- UniformGrid and ImplicitGrid build one-level indexes of non-intersecting bins.
- BVHTree and SpatialOctree build hierarchies of nested bounding boxes.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/spintop.html>`_


.. toctree::
   :maxdepth: 2
   :caption: Helper classes, single-level indexes

   helper
   spatial_ordering
   uniformgrid
   implicitgrid

.. toctree::
   :maxdepth: 2
   :caption: Tree-structure indexes

   bvh
   octree

