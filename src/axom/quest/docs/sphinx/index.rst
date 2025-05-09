.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Quest User Guide
================

The Quest component of Axom provides several spatial operations and queries
on a ``mint::Mesh``.

  - Operations

    - :ref:`Read a surface mesh<reading-mesh>` from an STL file
    - :ref:`Check for some common mesh errors; deduplicate vertices<check-and-repair>`

      - vertex welding: merge vertices closer than a specified distance
        "epsilon"
      - find self-intersections and degenerate triangles in a surface mesh
      - watertightness test: is a surface mesh a watertight manifold?

  - Point queries

    - Surface mesh point queries :ref:`in C<surface-query-c>` or
      :ref:`in C++<surface-query-cpp>`

      - in/out query: is a point inside or outside a surface mesh?
      - signed distance query: find the minimum distance from a query point
        to a surface mesh

    - :ref:`Point in cell query<point-in-cell>`: for a query point, find the
      cell of the mesh that holds the point and the point's isoparametric
      coordinates within that cell
    - :ref:`All nearest neighbors<all-nearest>`: given a list of point
      locations and regions, find all neighbors of each point in a different
      region
    - :ref:`Isosurface detection<isosurface-detection>`: generate an
      isosurface mesh from a nodal scalar field and an isovalue.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/questtop.html>`_


.. toctree::
   :caption: Contents
   :maxdepth: 2

   read_mesh
   check_and_repair
   point_mesh_query
   point_mesh_query_cpp
   point_in_cell
   all_nearest_neighbors
   isosurface_detection

