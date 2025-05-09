.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _reading-mesh:

*****************
Reading in a mesh
*****************

Applications commonly need to read a mesh file from disk.  Quest provides the
``STLReader`` class, which can read binary or ASCII `STL`_ files, as well as the
``PSTLReader`` class for use in parallel codes.  STL (stereolithography)
is a common file format for triangle surface meshes.  The STL reader classes
will read the file from disk and build a ``mint::Mesh`` object.  Quest also
provides the ``ProEReader`` class, for ASCII Pro/E files containing tetrahedra,
and the ``PProEReader`` class for use in parallel codes.  PTC Creo is a modeling
application formerly known as Pro/ENGINEER, and its file format is in use among
Axom's users.

.. _STL: https://en.wikipedia.org/wiki/STL_(file_format)

Reading an STL file
-------------------

The code examples are excerpts from the file ``<axom>/src/tools/mesh_tester.cpp``.

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

As read by Axom, an ASCII Pro/E tet file contains:

- Zero or more comment lines starting with a ``#`` character
- One line with two integers: the number of nodes ``n`` and the number of
  tetrahedra ``t``
- ``n`` lines, one for each node; each line contains a contiguous integer ID
  starting at 1 and three floating-point numbers specifying the node location
- ``t`` lines, one for each tetrahedron; each line contains a contiguous
  integer ID starting at 1 and four integers specifying the tet's nodes

Reading an ASCII Pro/E tet file is similar to reading an STL file.  The code
examples are excerpts from the file ``<axom>/src/axom/quest/examples/quest_proe_bbox.cpp``.
The Pro/E reader has the ability to read a subset of the mesh in the file,
defined by a user-supplied predicate function.  The example code shows how
to use a convenience function to specify a predicate that keeps only tets
fully included in a user-supplied bounding box.

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
