.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _isosurface-detection:

********************
Isosurface detection
********************

Quest generates isocontours from node-centered scalar fields on Conduit Blueprint meshes.
The fixed-stride output contains line segments in 2D or triangles in 3D
and records the input cell and domain for each element.

.. Note::

   The legacy backend implements the original algorithm:

   William E. Lorensen,  and Harvey E. Cline (1 August 1987).
   "Marching cubes: A high resolution 3D surface construction algorithm".
   *ACM SIGGRAPH Computer Graphics*. 21 (**4**): 163-169

.. Note::

   An isosurface saddle point makes the contour topology in a cell ambiguous.
   Each backend uses a fixed lookup table to resolve these cases,
   and their choices may differ.

.. figure:: figs/planar_and_spherical_isosurfaces.png
   :width: 400px

   Planar and spherical isocontours generated from
   :math:`f(\mathbf{r}) = f_0 + \mathbf{r} \cdot \mathbf{n}` and
   :math:`g(\mathbf{r}) = |\textbf{r} - \textbf{r}_0|`, respectively.
   Colors denote domain indices in the multi-domain cubic mesh.

The algorithm is implemented in the class ``quest::MarchingCubes``.

The inputs are:

#. The mesh containing the scalar field, in `Conduit Blueprint format
   <https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html>`__.
#. The name of the Blueprint topology to contour.
#. The name of the scalar field data within the input mesh.
#. The contour value.

The following example shows usage of the ``MarchingCubes`` class.
A complete example is in ``src/axom/quest/examples/quest_marching_cubes_example.cpp``.

Relevant header files:

.. sourcecode:: C++

   #include "axom/core.hpp"
   #include "axom/quest/MarchingCubes.hpp"
   #include "axom/mint/mesh/UnstructuredMesh.hpp"
   #include "conduit_relay_io_blueprint.hpp"

Set up the Blueprint mesh and the ``MarchingCubes`` object.

``MarchingCubes`` accepts single-domain and multi-domain Blueprint meshes.
A domain is one local part of a mesh. A multi-domain mesh may contain any
number of local domains, including zero.

You can pass a single-domain mesh directly. ``MarchingCubes`` wraps it internally.

Blueprint meshes have named topologies and fields.
The example uses the topology ``mesh`` and the nodal field ``scalarFieldName``.

The ``axom::runtime_policy::Policy::seq`` argument runs the extraction on the host.
Builds configured with OpenMP, CUDA, or HIP can use those policies instead.

The ``MarchingCubesDataParallelism`` constructor argument selects the scan
strategy used by the legacy structured-mesh backend. The Bump backend manages
its own parallelism and ignores this argument.

The two backends accept different mesh types.

Legacy backend
  This is the default. It accepts only a ``structured`` topology with an
  ``explicit`` coordset, including ghost-padded structured input.

Bump backend
  In a build configured with Bump, call ``setUseBumpBackend(true)`` before
  ``setMesh``. This backend accepts the legacy formats, ``uniform`` and
  ``rectilinear`` topologies, and single-shape unstructured meshes made of
  quadrilaterals in 2D or hexahedra in 3D.

  Explicit and rectilinear coordinate arrays must use ``float64``.
  Function fields must also use ``float64``, and mask fields must use ``int32``.

The Bump backend welds contour vertices, so adjacent facets share vertex IDs.
``populateContourMeshBlueprint`` and ``relinquishContourDataBlueprint``
return this representation. In 3D, Bump's native ``CutField`` output may contain
triangles, quadrilaterals, or polygons with more than four vertices.
The fixed-stride array and Mint APIs require triangles, so the adaptor
fan-triangulates each polygonal face for those outputs.

The backends differ in these ways:

#. *Precision.* The Bump intersector converts the ``float64`` function values
   to ``float`` and computes edge-crossing positions in single precision.
   The legacy backend uses ``double``.
#. *Ambiguity.* Like the legacy 1987 tables, Bump's VisIt-derived cut
   tables resolve ambiguous saddle cell configurations with a single
   fixed triangulation per case. The result may differ from the bilinear or
   trilinear interpolant. Neither backend uses an asymptotic decider or a
   three-way negative, zero, and positive classification.
#. *Fan triangulation.* Each polygonal Bump face is fan-triangulated from its
   first corner for the array and Mint outputs. For a non-planar polygon, the
   resulting triangle areas depend on the first corner.
   ``populateContourMeshBlueprint`` returns the original polygons unless its
   ``triangulate`` argument is true.

Bump classifies a corner with a strict ``>``, while the legacy backend uses
``>=``. To match the legacy behavior at nodal values, ``MarchingCubes`` passes
Bump the next lower ``float`` value. A direct
``axom::bump::extraction::CutField`` call does not apply this adjustment.

``MarchingCubesRobustnessPolicy::robust`` currently behaves the same as ``standard``.

.. sourcecode:: C++

   conduit::Node blueprintMesh = blueprint_mesh_from_user();
   axom::quest::MarchingCubes mc(
     axom::runtime_policy::Policy::seq,
     axom::getDefaultAllocatorID(),
     axom::quest::MarchingCubesDataParallelism::byPolicy);
   mc.setUseBumpBackend(true);
   mc.setMesh(blueprintMesh, "mesh");
   mc.setFunctionField("scalarFieldName");

Run the algorithm:

.. sourcecode:: C++

   double contourValue = 0.5;
   mc.computeIsocontour(contourValue);

Place the isocontour in an output ``axom::mint::UnstructuredMesh`` object:

``MarchingCubes`` generates the isocontour mesh in an internal format.
Use ``populateContourMesh`` to copy it to an
``axom::mint::UnstructuredMesh``. In 3D this method always produces triangles.
When the Bump backend is enabled,
``populateContourMeshBlueprint`` and ``relinquishContourDataBlueprint``
provide the welded Blueprint output directly.

Repeated calls to ``computeIsocontour`` append to the array and Mint outputs.
The Blueprint methods return only the most recent extraction for each input
domain. Call ``clearOutput`` before computing a replacement contour.

``populateContourMesh`` provides two scalar fields for the generated
mesh:

#. the ID of the cell from the input mesh that generated the
   isocontour cell.
#. the ID of the domain from the input mesh that generated the
   isocontour cell.

The names of these fields are user-specified.  Use empty strings if
you don't need these fields.  This example puts cell IDs in
"cellIds" and domain IDs in "domainIds".

.. sourcecode:: C++

   axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>
     contourMesh(3, axom::mint::TRIANGLE);
   mc.populateContourMesh(contourMesh, "cellIds", "domainIds");

After putting the isosurface in the ``UnstructuredMesh`` object,
the ``MarchingCubes`` object is no longer needed.

MPI-parallel runs
-----------------

Each MPI rank passes its local domains to ``MarchingCubes``. Extraction does
not communicate between ranks, and output node and cell IDs are unique only
within a rank. Applications that need globally unique IDs must renumber them.
