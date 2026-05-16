.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Quest User Guide
================

Axom's Quest component provides spatial queries, geometry readers, contouring,
and shaping algorithms for simulation workflows. Quest works with several Axom
mesh representations, including ``mint::Mesh`` objects, Conduit Blueprint
meshes, and MFEM meshes.

This guide focuses on the most common Quest workflows:

* :ref:`Read geometry and meshes <reading-mesh>` from STL, Pro/E, STEP,
  C2C, and MFEM-based inputs.
* :ref:`Check and repair surface meshes <check-and-repair>` before using
  algorithms that require watertight geometry.
* Run surface and mesh queries, including :ref:`surface containment and signed
  distance <surface-query-cpp>`, :ref:`point-in-cell <point-in-cell>`, and
  :ref:`all nearest neighbors <all-nearest>`.
* Generate :ref:`isocontours and isosurfaces <isosurface-detection>` from
  nodal scalar fields on Blueprint meshes.
* Build point-set simplicial meshes with :doc:`Quest's Delaunay triangulation <delaunay>`.
* Query curved and linearized shapes with :doc:`Winding Numbers <winding_number>`.
* Approximate curved contour geometry with :doc:`Linearize Curves <linearize_curves>`.
* Build :ref:`shaping pipelines <shaping-overview>` that convert Klee shape
  descriptions into material volume fractions on target meshes.

The Sphinx pages describe the user-facing workflows and show representative
examples from Quest's sources and tests. For a full API reference, use the
generated Doxygen documentation below.

API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/questtop.html>`_


.. toctree::
   :caption: Contents
   :maxdepth: 1

   read_mesh
   check_and_repair
   point_mesh_query
   point_mesh_query_cpp
   point_in_cell
   all_nearest_neighbors
   isosurface_detection
   delaunay
   winding_number
   linearize_curves
   shaping
