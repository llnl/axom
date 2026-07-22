BVH
^^^

The ``BVH`` class implements a
`bounding volume hierarchy <https://en.wikipedia.org/wiki/Bounding_volume_hierarchy>`_.
This data structure recursively subdivides a rectilinear region of interest
into a "tree" of subregions.

``BVH`` is well-suited for particle-mesh or ray-mesh intersection tests.
It is also well-suited to data sets where the contents are unevenly distributed,
since the bins are subdivided based on their contents.

The following code example shows how a ``BVH`` can be used to accelerate a
point-mesh intersection algorithm.  First, we generate bounding boxes for all
triangles in the mesh, and call ``BVH::initialize()`` with the bounding boxes.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvh_header_start
   :end-before: _bvh_header_end
   :language: C++

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvh_build_start
   :end-before: _bvh_build_end
   :language: C++

After the structure is built, we can use the BVH to generate a list of element IDs
that are candidate neighbors to the query points.  Call ``BVH::findPoints()``
to get the list of element IDs that the query point intersects.
The key idea of ``findPoints()`` is that testing for probe intersection with a
child node (bounding box) is cheap.  If a node intersection test fails (misses),
the child node, and all of its corresponding elements, can be skipped during the
BVH traversal.  If the probe does intersect a child node, the node's children are also
tested for probe intersection.  Without the acceleration data structure, each
probe point must be tested against each triangle.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvh_candidate_start
   :end-before: _bvh_candidate_end
   :language: C++

Note that the returned packed candidate intersection array (``candidatesPtr`` above)
needs to be deallocated by the caller.

When using candidate-query overloads that allocate ``axom::Array`` output, pass
an ``axom::HostAllocator`` when host staging or host fallback allocation should
be controlled explicitly. Query overloads without a host allocator remain
available as compatibility paths and use Axom's current default host allocator.

Finally, test the point against all candidate neighbor triangles.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvh_cand_int_start
   :end-before: _bvh_cand_int_end
   :language: C++

Device Traversal API
--------------------

The ``BVH`` class contains a ``getTraverser()`` method, which returns an object
that can be used to traverse a BVH with user-defined actions.

The returned traverser type has one function, ``traverse_tree()``, which takes the
following arguments:

- ``const QueryObject& p``: the object to traverse the BVH with. This is passed
  into each invocation of the traversal predicate.
- ``LeafAction&& lf``: a function or lambda which is executed on each leaf node
  of the BVH that is reached during traversal. It should take in two arguments,
  the index of the leaf node in the BVH, and a pointer to an array mapping leaf
  node indices to the original index of elements.
- ``Predicate&& predicate``: a function which determines whether to traverse
  down to a given internal node. It should take in two arguments: the query
  object, and the tentative node's bounding box.

This object may be used within a CUDA kernel, so long as the execution space
parameter of ``BVH`` is set correctly.

This method can be used to avoid the extra memory allocation needed for holding an
array of candidate intersections. For example, if we only wish to count the number
of intersections for a given query point, our leaf action could get the underlying
mesh element based on the element index, check whether the query point intersects it
and then increment a per-point counter.

This method also allows uses of the BVH beyond intersection testing.
``quest::SignedDistance`` uses the BVH traversal object to search for the closest
surface elements to a query point. The leaf action that is used checks each candidate
leaf against a current-minimum candidate; if closer, the current-minimum candidate
is set to the new surface element. The predicate used for traversal also utilizes the
current-minimum candidate data to avoid traversing internal nodes that are farther than
the current minimum squared distance.

Example: Broad-phase collision detection
========================================

The following example tests each element of a surface mesh for intersection with
each other element. It provides an example of how the device traversal object
might be used in a broad-phase collision detection problem. First, we initialize
the BVH with the bounding boxes of all the query objects (mesh elements), and
create a traverser object:

.. literalinclude:: ../../../quest/examples/quest_bvh_two_pass.cpp
   :start-after: _bvh_traverse_init_start
   :end-before: _bvh_traverse_init_end
   :language: C++

Next, we define a traversal predicate. In this case, since we are testing objects
in the mesh against each other, our query objects are bounding boxes. Thus, the
traversal predicate tests if the query bounding box intersects with the bounding
boxes of nodes in the BVH:

.. literalinclude:: ../../../quest/examples/quest_bvh_two_pass.cpp
   :start-after: _bvh_traverse_predicate_start
   :end-before: _bvh_traverse_predicate_end
   :language: C++

Since we do not know the total number of candidate intersections yet, we must
traverse the BVH twice; for the first traversal we count the number of candidate
intersections for each query object, allowing us to compute offset indices and
total storage requirements for the collision pairs:

.. literalinclude:: ../../../quest/examples/quest_bvh_two_pass.cpp
   :start-after: _bvh_traverse_first_pass_start
   :end-before: _bvh_traverse_first_pass_end
   :language: C++

After computing offset indices and allocating output arrays, we can then perform
a second traversal through the BVH. This time, we will store candidate collision
pairs when we reach a leaf node:

.. literalinclude:: ../../../quest/examples/quest_bvh_two_pass.cpp
   :start-after: _bvh_traverse_second_pass_start
   :end-before: _bvh_traverse_second_pass_end
   :language: C++

The result of the two-pass query is a list of candidate collision pairs. A code
could then do further operations on these pairs of elements, such as test them
for actual intersection.
