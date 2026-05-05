.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _reading-mesh:

*******************
Reading in geometry
*******************

Quest contains several readers that translate geometry files into Axom data
structures. The most common cases are STL surface meshes, Pro/E tetrahedral
meshes, and geometry inputs used by shaping workflows such as STEP, C2C, and
MFEM contour files.

The STL and Pro/E readers produce ``mint::Mesh`` objects directly. Other
readers expose geometry in forms that are better suited to downstream Quest
algorithms, such as NURBS patches, NURBS curves, or curved polygons.

Quest currently provides the following reader families:

* ``STLReader`` and ``PSTLReader`` for ASCII or binary STL triangle meshes.
* ``ProEReader`` and ``PProEReader`` for ASCII Pro/E tetrahedral meshes.
* ``STEPReader`` and ``PSTEPReader`` for STEP B-Rep geometry represented as
  trimmed NURBS patches, with optional triangulated output.
* ``C2CReader`` and ``PC2CReader`` for C2C contour files represented as
  NURBS curves, when Axom is built with the C2C dependency.
* ``MFEMReader`` for MFEM contour files represented as curves or curved
  polygons, when Axom is built with MFEM support.

.. _STL: https://en.wikipedia.org/wiki/STL_(file_format)

Reading an STL file
-------------------

STL (stereolithography) is a common file format for triangle surface meshes.
Quest's STL readers load the file from disk and populate an
``mint::UnstructuredMesh`` containing triangles. STL stores triangles as a
"triangle soup", so downstream algorithms often need a cleanup pass to weld
duplicate vertices and check for defects. The next page describes that
workflow.

The following example is excerpted from ``<axom>/src/tools/mesh_tester.cpp``.

We include the STL reader header

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _read_stl_include1_start
   :end-before: _read_stl_include1_end
   :language: C++

and also the mint Mesh and UnstructuredMesh headers.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _read_stl_include2_start
   :end-before: _read_stl_include2_end
   :language: C++

For convenience, we use typedefs in the axom namespace.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _read_stl_typedefs_start
   :end-before: _read_stl_typedefs_end
   :language: C++

The following example shows usage of the STLReader class:

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _read_stl_file_start
   :end-before: _read_stl_file_end
   :language: C++

After reading the STL file, the ``STLReader::getMesh`` method gives access to the
underlying mesh data.  The reader may then be deleted.

Reading a Pro/E file
--------------------

Quest's ``ProEReader`` reads ASCII Pro/E tetrahedral meshes and can optionally
filter the tetrahedra during input.

As read by Axom, an ASCII Pro/E tet file contains:

- Zero or more comment lines starting with a ``#`` character
- One line with two integers: the number of nodes ``n`` and the number of
  tetrahedra ``t``
- ``n`` lines, one for each node; each line contains a contiguous integer ID
  starting at 1 and three floating-point numbers specifying the node location
- ``t`` lines, one for each tetrahedron; each line contains a contiguous
  integer ID starting at 1 and four integers specifying the tet's nodes

Reading an ASCII Pro/E tet file is similar to reading an STL file. The code
examples are excerpts from ``<axom>/src/axom/quest/examples/quest_proe_bbox.cpp``.
The example demonstrates one of the reader's useful features: selecting a
subset of the input tetrahedra using a predicate. In this case, the predicate
keeps only tetrahedra whose nodes fall inside a user-supplied bounding box.

We include the ProEReader header

.. literalinclude:: ../../examples/quest_proe_bbox.cpp
   :start-after: _read_proe_include1_start
   :end-before: _read_proe_include1_end
   :language: C++

and also the mint Mesh and UnstructuredMesh headers.

.. literalinclude:: ../../examples/quest_proe_bbox.cpp
   :start-after: _read_proe_include2_start
   :end-before: _read_proe_include2_end
   :language: C++

For convenience, we specify some type aliases.

.. literalinclude:: ../../examples/quest_proe_bbox.cpp
   :start-after: _read_proe_typealiases_start
   :end-before: _read_proe_typealiases_end
   :language: C++

The following example shows how to use the ProEReader class.
Calling ``reader.setTetPredFromBoundingBox(bbox, false)``, as shown in the
code, makes a tetrahedron predicate that accepts tets with all four nodes
falling in ``bbox`` and rejects others.  Alternately, the user can specify
an arbitrary predicate function with ``setTetPred()``.  If the user specifies
no tetrahedron predicate, the reader reads all tets in the file.

.. literalinclude:: ../../examples/quest_proe_bbox.cpp
   :start-after: _read_proe_file_start
   :end-before: _read_proe_file_end
   :language: C++

After reading the Pro/E file, the ``ProEReader::getMesh`` method gives access
to the underlying mesh data.  The reader may then be deleted.

Other Quest readers
-------------------

Quest includes several other readers that are commonly used with shaping and
CAD-oriented workflows.

``STEPReader``
^^^^^^^^^^^^^^

``STEPReader`` reads trimmed STEP surfaces using Open Cascade. It can expose
the model as NURBS patches and trimming curves, query metadata such as patch
IDs and bounding boxes, and generate a triangulated ``mint::UnstructuredMesh``
approximation of the model when a surface mesh is needed.

This reader is available when Axom is configured with Open Cascade support.
The following example, excerpted from
``<axom>/src/axom/quest/examples/quest_winding_number_3d.cpp``, reads a STEP
file, loads its trimmed NURBS patches, and queries the model bounding box.

.. literalinclude:: ../../examples/quest_winding_number_3d.cpp
   :start-after: _read_step_file_start
   :end-before: _read_step_file_end
   :language: C++

When a triangle mesh approximation is needed, ``STEPReader`` can also
triangulate the loaded B-Rep:

.. literalinclude:: ../../examples/quest_winding_number_3d.cpp
   :start-after: _read_step_triangulate_start
   :end-before: _read_step_triangulate_end
   :language: C++

``C2CReader``
^^^^^^^^^^^^^

``C2CReader`` reads contour files and stores the result as NURBS curves. This
reader is available when Axom is configured with the C2C dependency and is
primarily used by shaping workflows that revolve or sample contour geometry.

The following example, excerpted from
``<axom>/src/axom/quest/examples/containment_driver.cpp``, reads a contour file
and linearizes the resulting curves into a segment mesh:

.. literalinclude:: ../../examples/containment_driver.cpp
   :start-after: _read_c2c_file_start
   :end-before: _read_c2c_file_end
   :language: C++

``MFEMReader``
^^^^^^^^^^^^^^

``MFEMReader`` reads MFEM contour files and can return either individual curves
or grouped curved polygons. It is available when Axom is configured with MFEM
support. Quest uses this representation in workflows such as sampling-based
shaping with winding-number containment tests.

The following example, excerpted from
``<axom>/src/axom/quest/examples/quest_winding_number_2d.cpp``, reads an MFEM
contour file into an array of NURBS curves:

.. literalinclude:: ../../examples/quest_winding_number_2d.cpp
   :start-after: _read_mfem_file_start
   :end-before: _read_mfem_file_end
   :language: C++
