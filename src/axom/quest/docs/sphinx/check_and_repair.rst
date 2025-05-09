.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _check-and-repair:

***********************
Check and repair a mesh
***********************

The STL file format specifies triangles without de-duplicating their
vertices.  Vertex welding is needed for several mesh algorithms that
require a watertight manifold.  Additionally, mesh files often contain
errors and require some kind of cleanup.
The following code examples are excerpts from the file
``<axom>/src/tools/mesh_tester.cpp``.

Quest provides a function to weld vertices within a distance of some
specified epsilon.  This function takes arguments
``mint::UnstructuredMesh< mint::SINGLE_SHAPE > **surface_mesh`` and
``double epsilon``, and modifies ``surface_mesh``.  In addition to
the mint Mesh and UnstructuredMesh headers (see previous page), we
include the headers declaring the functions for checking and repairing
surface meshes.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_repair_include_start
   :end-before: _check_repair_include_end
   :language: C++

The function call itself:

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_repair_weld_start
   :end-before: _check_repair_weld_end
   :language: C++

One problem that can occur in a surface mesh is self-intersection.  A
well-formed mesh will have each triangle touching the edge of each of its
neighbors.  Intersecting or degenerate triangles can cause problems for some
spatial algorithms.  To detect such problems using Quest, we first make
containers to record defects that might be found.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_repair_intersections_containers_start
   :end-before: _check_repair_intersections_containers_end
   :language: C++

Then, we call the function to detect self-intersections and degenerate
triangles.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_repair_intersections_start
   :end-before: _check_repair_intersections_end
   :language: C++

After calling ``findTriMeshIntersections``, ``collisions`` will hold
the indexes of each pair of intersecting triangles and ``degenerate`` will
contain the index of each degenerate triangle.  The user code can then address or
report any triangles found.  Mesh repair beyond welding close vertices is
beyond the scope of the Quest component.

Check for watertightness
------------------------

Before using Quest's surface point queries, a mesh must be watertight, with
no cracks or holes.  Quest provides a function to test for watertightness,
declared in the same header file as the tests self-intersection and an enum
indicating watertightness of a mesh.
If the code is working with a mesh read in from an STL file,
*weld the vertices* (see above) before checking for watertightness!

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_watertight_start
   :end-before: _check_watertight_end
   :language: C++

This routine builds the face relation of the supplied triangle surface mesh.
The face of a triangle is a one-dimensional edge.  If the mesh is
big, building the face relation may take some time.  Once built, the routine
queries face relation: each edge
of every triangle must be incident in two triangles.  If the mesh has a
defect where more than two triangles share an edge, the routine returns
``CHECK_FAILED``.  If the mesh has a hole, at least one triangle edge
is incident in only one triangle and the routine returns ``NOT_WATERTIGHT``.
Otherwise, each edge is incident in two triangles, and the routine returns
``WATERTIGHT.``

After testing for watertightness, report the result.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _report_watertight_start
   :end-before: _report_watertight_end
   :language: C++

After an STL mesh has

  - been read in with ``STLReader``,
  - had vertices welded using ``weldTriMeshVertices()``,
  - contains no self-intersections as reported by ``findTriMeshIntersections()``,
  - and is watertight as reported by ``isSurfaceMeshWatertight()``,

the in-out and distance field queries will work as designed.
