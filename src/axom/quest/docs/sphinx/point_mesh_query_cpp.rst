.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _surface-query-cpp:

***********************************
Surface mesh point queries: C++ API
***********************************

Codes written in C++ may use the object-oriented C++ APIs to perform in/out and
signed distance queries.  In addition to language choice, the C++ API lets a code work
with more than one mesh at the same time.  Unlike the C API, the C++ API for
in/out and signed distance queries has no initializer taking a file name:
readying the mesh is a separate, prior step.

In/out Octree
-------------

The C++ in/out query is provided by the ``quest::InOutOctree`` class, from the
following header.  See ``<axom>/src/axom/quest/tests/quest_inout_octree.cpp``.

.. literalinclude:: ../../tests/quest_inout_octree.cpp
   :start-after: _quest_inout_cpp_include_start
   :end-before: _quest_inout_cpp_include_end
   :language: C++

Some type aliases are useful for the sake of brevity.  The class is templated on
the dimensionality of the mesh.  Currently, only meshes in 3D are supported;
here ``DIM`` equals 3.

.. literalinclude:: ../../tests/quest_inout_octree.cpp
   :start-after: _quest_inout_cpp_typedef_start
   :end-before: _quest_inout_cpp_typedef_end
   :language: C++

Instantiate the object using ``GeometricBoundingBox bbox`` and a mesh, and
generate the index.

.. literalinclude:: ../../tests/quest_inout_octree.cpp
   :start-after: _quest_inout_cpp_init_start
   :end-before: _quest_inout_cpp_init_end
   :language: C++

Test a query point.
::

   SpacePt pt = SpacePt::make_point(2., 3., 1.);
   bool inside = octree.within(pt);

All cleanup happens when the index object's destructor is called 
(in this case, when the variable ``octree`` goes out of scope).

Signed Distance
---------------

The C++ signed distance query is provided by the ``quest::SignedDistance`` class,
which wraps an instance of ``primal::BVHTree``.
Examples from ``<axom>/src/axom/quest/tests/quest_signed_distance.cpp``.

Class header:

.. literalinclude:: ../../tests/quest_signed_distance.cpp
   :start-after: _quest_distance_cpp_include_start
   :end-before: _quest_distance_cpp_include_end
   :language: C++

The constructor takes several arguments:

- ``const mint::Mesh* surfaceMesh``: A pointer to a surface mesh with triangles
  and/or quadrilaterals.
- ``bool isWatertight``: Indicates the mesh is a watertight mesh, a manifold.
  The signed distance from a point to a manifold is mathematically well-defined.
  When the input is not a closed surface mesh, the mesh must span the entire
  computational mesh domain, dividing it into two regions.
- ``bool computeSign`` (default ``true``): Optional. Enables or disables the
  computation of signs in distance queries.
- ``int allocatorID``: Optional. Sets a custom Umpire allocator to use in
  constructing the underlying BVH; by default, this is set to a default allocator
  for the execution space the ``SignedDistance`` class is instantiated in
  (host-side memory for CPU and OpenMP, unified memory for GPUs).

Note that the second and subsequent arguments to the constructor correspond to
``quest::signed_distance_set`` functions in the C API.

As with the ``InOutOctree``, the class is templated on the dimensionality
of the mesh, with only 3D meshes being supported. The class also accepts a
template parameter for execution space, for running signed distance queries with
OpenMP or on a GPU.

.. literalinclude:: ../../tests/quest_signed_distance.cpp
   :start-after: _quest_distance_cpp_init_start
   :end-before: _quest_distance_cpp_init_end
   :language: C++

Test a query point.

.. code-block:: C++

   axom::primal::Point< double,3 > pt =
     axom::primal::Point< double,3 >::make_point(2., 3., 1.);
   double signedDistance = signed_distance.computeDistance(pt);

Test a batch of query points.

.. code-block:: C++

   const int numPoints = 20;
   axom::primal::Point<double, 3>* pts =
     axom::allocate<axom::primal::Point<double, 3>>(numPoints);
   for (int ipt = 0; ipt < numPoints; ipt++)
   {
     // fill pts array
     pts[ipt] = ...;
   }
   double signedDists = axom::allocate<double>(20);
   signed_distance.computeDistances(numPts, pts, signedDists);

The object destructor takes care of all cleanup.
