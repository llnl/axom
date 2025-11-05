# Lesson 04: Acceleration squared: Accelerating our device-based self-intersection query with spatial indexes

In this lesson, we add a GPU-enabled spatial index, the Bounding Volume Hierarchy (BVH tree) to further accelerate our triangle mesh self-intersection algorithm.

### BVH trees
A BVH tree is a nested spatial index over the objects in our mesh. Specifically, we build a linear BVH tree in a bottom-up fashion over the axis-aligned bounding boxes of the triangles in our mesh. We then use a built-in batch query to find all overlapping bounding boxes for a collection of input bounding boxes.

> :information_source: Our linear BVH tree implementation is adapted from Lauterbach et al's "Fast BVH Construction on GPUs" (Computer Graphics Forum 28:375--384 2009) and Karras's "Maximizing parallelism in the construction of BVHs, octrees, and k-d trees." (High-Performance Graphics, 2012).

> :bulb:  Axom provides a GPU-accelerated BVH tree as well as several other spatial indexes, such as ``Octree``s and ``ImplicitGrid``s in its ``spin`` component.

### Approach
Our BVH-based algorithm consists of three phases:
1. We use the BVH tree to find, for each triangle ``t``, the index of all candidate triangles whose bounding boxes overlap that of ``t``. This uses the relatively inexpensive ``primal::intersect(BoundingBox, BoundingBox)`` predicate
2. We then filter the list to remove degenerate triangles, and to remove duplicate pairs
3. Finally, from the above list of potential candidates, we retain the pairs of triangles that actually intersect using the relatively expensive ``primal::intersect(Triangle, Triangle)`` predicate

> :key:  This approach significantly reduces the computational work for reasonable meshes, such as those that Axom users would use for simulations, since each triangle's bounding box is expected to overlap with only a few others in the mesh. Thus, to find the candidates, we expect to only check ``O(lg(N))`` other triangles, rather than ``O(N)`` triangles for a mesh with ``N`` triangles.

### First step: Finding overlapping bounding boxes
In our first step, we first initialize a BVH tree over a collection of bounding boxes:
```cpp
  axom::spin::BVH<3, ExecSpace, double> bvh;
  bvh.setAllocatorID(kernel_allocator);
  bvh.initialize(bbox_v, bbox_v.size());
```
and then query the tree with a collection of bounding boxes using the built-in ``findBoundingBoxes()`` function.

```cpp
  IndexArray offsets_d(bbox_v.size(), bbox_v.size(), kernel_allocator);
  IndexArray counts_d(bbox_v.size(), bbox_v.size(), kernel_allocator);
  IndexArray candidates_d(0, 0, kernel_allocator);

  auto offsets_v = offsets_d.view();
  auto counts_v = counts_d.view();
  bvh.findBoundingBoxes(offsets_v, counts_v, candidates_d, bbox_v.size(), bbox_v);
```
The input to this query is a collection of bounding boxes, and the output is returned as a CSR-like collection of arrays indexed by the input bounding boxes (which correspond to the triangles from the mesh). After the call:
* ``counts[idx]`` contains the number of entities that intersect with input entity ``idx``
* ``offsets[idx]`` contains the offset within the returned ``candidates`` array of the first result of input entity ``idx``
* The ``candidates`` array contains all the potential intersection candidates.

> :memo: In this case, the input bounding boxes are identical to the query bounding boxes, but this is not generally the case

### Second step: Filtering and deduplicating the candidates
The results from the previous step will, in general, contain numerous duplicates as well as the indices of degenerate triangles. For example, since every bounding box intersects itself, the results from the previous query contain all pairs ``(idx, idx)``  for a triangle with index ``idx``.  For pairs of distinct triangles ``idx1`` and ``idx2`` whose bounding boxes overlap, the results will contain both ``(idx1, idx2)`` and ``(idx2, idx1)``

In this phase, we filter out these duplicates, as well as any pairs containing a degenerate triangle. This is aided by building up a device-accessible array of boolean flags:
```cpp
    // Compute a device bool array of validity flags
    axom::Array<bool> is_valid_d(axom::ArrayOptions::Uninitialized {},
                                 bbox_v.size(),
                                 kernel_allocator);
    auto is_valid_v = is_valid_d.view();

    axom::for_all<ExecSpace>(
      totalTriangles,
      AXOM_LAMBDA(axom::IndexType i) { is_valid_v[i] = !tris_v[i].degenerate(); });
```

Note that our input is a CSR-like array, which is not generally a great representation for GPU-like processors. As we iterate through this representation, we generate a more performant representation for our output in the form of a pair of corresponding arrays, in which each potential candidate pair is represented at position ``idx`` of the ``indices`` and ``validCandidates`` arrays. 

Our kernel therefore has the following implementation:
```cpp
  // Expand candidate list into corresponding arrays of indices
  // Only keep results where candidate has a greater index than triangle
  // and when both are non-degenerate
  IndexArray indices_d = ...;
  IndexArray validCandidates_d = ...;
  axom::IndexType numCandidates {};
  {
    IndexArray numValidCandidates_d(1, 1, kernel_allocator);
    numValidCandidates_d.fill(0);
    auto* numValidCandidates_p = numValidCandidates_d.data();

    // Keep pairs of valid triangles whose bounding boxes overlap
    axom::for_all<ExecSpace>(
      triMesh.numTriangles(),
      AXOM_LAMBDA(axom::IndexType i) {
        for(int j = 0; j < counts_v[i]; j++)
        {
          const axom::IndexType potential = candidates_v[offsets_v[i] + j];
          if(i < potential && is_valid_v[i] && is_valid_v[potential])
          {
            const auto idx = axom::atomicAdd<axom::auto_atomic>(numValidCandidates_p,
                                                                axom::IndexType {1});
            indices_v[idx] = i;
            validCandidates_v[idx] = potential;
          }
        }
      });

    axom::copy(&numCandidates,numValidCandidates_p, sizeof(axom::IndexType));
  }
```

### Third phase
At this point, we have a pair of corresponding arrays, allowing us to run a fairly straightforward kernel to check for actual intersections. As in the previous kernel, we use the ``axom::atomicAdd`` reduction to obtain a unique location for each intersection pair:
```cpp
  // Iterate through valid candidates to find actual intersections
  IndexArray intersect_d[2] = {IndexArray(axom::ArrayOptions::Uninitialized {},
                                          numCandidates,
                                          kernel_allocator),
                               IndexArray(axom::ArrayOptions::Uninitialized {},
                                          numCandidates,
                                          kernel_allocator)};
  axom::IndexType numIntersections {};
  {
    auto intersect1_v = intersect_d[0].view();
    auto intersect2_v = intersect_d[1].view();

    IndexArray numIntersections_d(1, 1, kernel_allocator);
    auto* numIntersections_p = numIntersections_d.data();

    auto indices_v = indices_d.view();
    auto validCandidates_v = validCandidates_d.view();

    // Perform triangle-triangle tests
    axom::for_all<ExecSpace>(
      numCandidates,
      AXOM_LAMBDA(axom::IndexType i) {
        constexpr bool includeBoundaries = false;
        const auto index = indices_v[i];
        const auto candidate = validCandidates_v[i];
        if(axom::primal::intersect(tris_v[index],
                                   tris_v[candidate],
                                   includeBoundaries,
                                   tol))
        {
          const auto idx = axom::atomicAdd<axom::auto_atomic>(numIntersections_p, axom::IndexType {1});
          intersect1_v[idx] = index;
          intersect2_v[idx] = candidate;
        }
      });

    axom::copy(&numIntersections, numIntersections_p, sizeof(axom::IndexType));
  }
  intersect_d[0].resize(numIntersections);
  intersect_d[1].resize(numIntersections);
```

The only remaining step is to copy the results from the ``intersect_d`` device arrays to the ``axom::Array`` of pairs of indices on the host, as we did in the previous example.

### Results
In this section we'll analyze our performance results. As we'll see, the spatial index considerably accelerates our query, and we are now able to efficiently perform self intersection tests on much larger triangle meshes.

The following uses ``car.stl``, a mesh with nearly 2 million triangles.

#### Comparing sequential and OpenMP results on a Sapphire Rapids node

```shell
>./bin/lesson_04_device_spatial_indexes -i car.stl --policy seq -v

[lesson_04: INFO] Loading the mesh took 0.154 seconds. 
[lesson_04: INFO] Vertex welding took 1.22 seconds. 
[lesson_04: INFO] After welding, mesh has 973,913 vertices and 1,947,910 triangles. 
[lesson_04: INFO] Running naive intersection algorithm in execution Space: [SEQ_EXEC] 
[lesson_04: INFO] 0: Initializing BVH took 0.518 seconds. 
[lesson_04: INFO] 1: Querying candidate bounding boxes took 3.06 seconds. 
[lesson_04: INFO] 2: Filtering invalid candidates took 0.195 seconds. 
[lesson_04: INFO] 3: Finding actual intersections took 0.82 seconds. 

[lesson_04: INFO] Stats for self-intersection query
    -- Number of mesh triangles 1,947,910
    -- Total possible candidates 947,648,911
    -- Candidates from BVH query 40,823,308
    -- Potential candidates after filtering 19,437,142
    -- Actual intersections 1,188
     
[lesson_04: INFO] Computing intersections with a BVH tree took 4.61 seconds. 
[lesson_04: INFO] Mesh had 1,188 intersection pairs 
```

```shell
>./bin/lesson_04_device_spatial_indexes -i car.stl --policy omp -v

[lesson_04: INFO] Loading the mesh took 0.992 seconds. 
[lesson_04: INFO] Vertex welding took 1.17 seconds. 
[lesson_04: INFO] After welding, mesh has 973,913 vertices and 1,947,910 triangles. 
[lesson_04: INFO] Running naive intersection algorithm in execution Space: [OMP_EXEC] 
[lesson_04: INFO] 0: Initializing BVH took 0.544 seconds. 
[lesson_04: INFO] 1: Querying candidate bounding boxes took 0.124 seconds. 
[lesson_04: INFO] 2: Filtering invalid candidates took 3.77 seconds. 
[lesson_04: INFO] 3: Finding actual intersections took 0.0225 seconds. 

[lesson_04: INFO] Stats for self-intersection query
    -- Number of mesh triangles 1,947,910
    -- Total possible candidates 947,648,911
    -- Candidates from BVH query 40,823,308
    -- Potential candidates after filtering 19,437,142
    -- Actual intersections 1,188
     
[lesson_04: INFO] Computing intersections with a BVH tree took 4.51 seconds. 
[lesson_04: INFO] Mesh had 1,188 intersection pairs 
```

#### Same run on an NVidia V100 GPU with a Power9 CPU

```shell
# sequential
>./bin/lesson_04_device_spatial_indexes -i car.stl --policy seq -v

[lesson_04: INFO] Loading the mesh took 0.434 seconds. 
[lesson_04: INFO] Vertex welding took  9.6 seconds. 
[lesson_04: INFO] After welding, mesh has 973,913 vertices and 1,947,910 triangles. 
[lesson_04: INFO] Running naive intersection algorithm in execution Space: [SEQ_EXEC] 
[lesson_04: INFO] 0: Initializing BVH took 0.779 seconds. 
[lesson_04: INFO] 1: Querying candidate bounding boxes took 10.5 seconds. 
[lesson_04: INFO] 2: Filtering invalid candidates took 0.39 seconds. 
[lesson_04: INFO] 3: Finding actual intersections took 2.18 seconds. 

[lesson_04: INFO] Stats for self-intersection query
    -- Number of mesh triangles 1,947,910
    -- Total possible candidates 947,648,911
    -- Candidates from BVH query 40,823,308
    -- Potential candidates after filtering 19,437,142
    -- Actual intersections 1,188
     
[lesson_04: INFO] Computing intersections with a BVH tree took 13.9 seconds. 
[lesson_04: INFO] Mesh had 1,188 intersection pairs 
```

```shell 
# cuda
>./bin/lesson_04_device_spatial_indexes -i car.stl --policy cuda -v
       
[lesson_04: INFO] Loading the mesh took 0.143 seconds. 
[lesson_04: INFO] Vertex welding took 9.61 seconds. 
[lesson_04: INFO] After welding, mesh has 973,913 vertices and 1,947,910 triangles. 
[lesson_04: INFO] Running naive intersection algorithm in execution Space: [CUDA_EXEC] 
[lesson_04: INFO] 0: Initializing BVH took 0.0265 seconds. 
[lesson_04: INFO] 1: Querying candidate bounding boxes took 0.0277 seconds. 
[lesson_04: INFO] 2: Filtering invalid candidates took 0.176 seconds. 
[lesson_04: INFO] 3: Finding actual intersections took 0.00603 seconds. 

[lesson_04: INFO] Stats for self-intersection query
    -- Number of mesh triangles 1,947,910
    -- Total possible candidates 947,648,911
    -- Candidates from BVH query 40,823,308
    -- Potential candidates after filtering 19,437,142
    -- Actual intersections 1,188
     
[lesson_04: INFO] Computing intersections with a BVH tree took 0.425 seconds. 
[lesson_04: INFO] Mesh had 1,188 intersection pairs 
```
We observe that kernel in phase 2, which operates on a ragged CSR-like array has the worst performance due to load imbalance among the threads.

We also observe that the "vertex welding" step is now a bottleneck for our self-intersection algorithm. This would be the next algorithm to consider as we continue to port and optimize this algorithm.

### Conclusion
In this tutorial, we showcased several of Axom's capabilities in the development of an efficient performance portable spatial query to find the self-intersections within a triangle mesh.

Among other things,
* We utilized Axom's ``slic`` logger for simplified logging and debugging throughout our application
* We explored Axom's device-compatible containers, like ``axom::Array`` and ``axom::ArrayView`` and how they work with Axom's abstractions over Umpire and RAJA.
* We utilized Axom's computational geometry primitives and operators for ``primal`` and its spatial indexes from ``spin``

As we saw, it is often not sufficient to simply port existing algorithms to the GPU using abstractions like RAJA and Umpire. We need to consider improved algorithms.

### Next steps
The code in this tutorial was adapted from Axom's [``mesh_tester`` application](https://github.com/LLNL/axom/blob/develop/src/tools/mesh_tester.cpp), which performs additional checks on triangle meshes, e.g. to determine if the mesh is watertight, and also provides implements the self-intersection algorithm against two other spatial indexes: a [``UniformGrid``](https://github.com/LLNL/axom/blob/develop/src/axom/spin/UniformGrid.hpp) and an [``ImplicitGrid``](https://github.com/LLNL/axom/blob/develop/src/axom/spin/ImplicitGrid.hpp).
