.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Spin User Guide
===============

The Spin component of Axom provides several index data structures to accelerate
spatial queries.  The Morton code classes relate each point in a region of
interest to a point on a one-dimensional space filling curve, and the
RectangularLattice helps in the computation of bin coordinates.  The UniformGrid
and ImplicitGrid classes build one-level indexes of non-intersecting bins, while
the BVHTree and SpatialOctree classes build nesting hierarchies of bounding
boxes indexing a region of interest.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/spintop.html>`_


.. toctree::
   :maxdepth: 2
   :caption: Helper classes, single-level indexes

   helper
   uniformgrid
   implicitgrid

.. toctree::
   :maxdepth: 2
   :caption: Tree-structure indexes

   bvh
   octree

