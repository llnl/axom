# Lesson 02: Naive self-intersection tests on our triangle mesh

In this lesson, we will implement a naive self-intersection algorithm on pairs of triangles within our triangle mesh.

Since the STL format is a "soup of triangles", we begin by adding a "welding" capability for the triangles in our mesh. Along the way, we utilize Axom's internal timer class.

## An aside about STL meshes

Due to its simplicity, the [STL format](https://en.wikipedia.org/wiki/STL_(file_format)) is a popular way of encoding and exchanging triangle meshes. It encodes a triangle mesh as triplets of coordinates in space in either an ASCII or a binary format. 

For example, a tetrahedron might be provided in the ASCII format as:
```
solid STL 
  facet normal -1 1 1
    outer loop
      vertex   1  1  1
      vertex  -1  1 -1
      vertex  -1 -1  1
    endloop
  endfacet
  facet normal  1 -1  1
    outer loop
      vertex   1  1  1
      vertex  -1 -1  1
      vertex   1 -1 -1
    endloop
  endfacet
  facet normal  1 1 -1
    outer loop
      vertex   1  1  1
      vertex   1 -1 -1
      vertex  -1  1 -1
    endloop
  endfacet
  facet normal -1 -1 -1
    outer loop
      vertex   1 -1 -1
      vertex  -1 -1  1
      vertex  -1  1 -1
    endloop
  endfacet
endsolid
```

Since the STL format has no means of ensuring that the coordinates of vertices in vertex-adjacent or edge-adjacent triangles actually match, we typically need to perform some preprocessing to "weld" or "unify" nearby vertices of an STL mesh. Otherwise, our intersection queries are likely to have false positives. However, one must also be careful when choosing the welding threshold, since in the extreme, we could end up welding all vertices together, yielding a mesh with degenerate triangles.

<details>
If your markdown renderer supports the ``stl`` tag, this might might look like:

```stl
solid STL 
  facet normal -1 1 1
    outer loop
      vertex   1  1  1
      vertex  -1  1 -1
      vertex  -1 -1  1
    endloop
  endfacet
  facet normal  1 -1  1
    outer loop
      vertex   1  1  1
      vertex  -1 -1  1
      vertex   1 -1 -1
    endloop
  endfacet
  facet normal  1 1 -1
    outer loop
      vertex   1  1  1
      vertex   1 -1 -1
      vertex  -1  1 -1
    endloop
  endfacet
  facet normal -1 -1 -1
    outer loop
      vertex   1 -1 -1
      vertex  -1 -1  1
      vertex  -1  1 -1
    endloop
  endfacet
endsolid
```
</details>

## Welding the triangles in our mesh and timing the operation

Axom provides an internal "welding" capability that treats triangles within a given threshold distance as identical. For this tutorial, we will expose the welding threshold as a command line argument ``--weld-threshold`` with a default value of `1e-6` and treat the following function from Axom's ``quest`` component as a black box:
```cpp
axom::quest::weldTriMeshVertices(&surface_mesh, weldThreshold);
```

We note that Axom provides a ``Timer`` class (based on ``std::chrono``) to assist in benchmarking sections of code. We can use this to see how long it takes to weld the triangle mesh, and use this to determine if it is sufficiently fast for our purposes: 
```cpp
axom::utilities::Timer timer;

timer.start();
axom::quest::weldTriMeshVertices(&surface_mesh, weldThreshold);
timer.stop();

SLIC_INFO(axom::fmt::format("Vertex welding took {:.3} seconds.",
                            timer.elapsedTimeInSec()));
```

## Primal: Computational geometry primitives and operations
Axom's ``primal`` component provides robust computational geometry primitives as well as operators between pairs of primitives.
In this lesson, we will be using the following primitives, which are defined as templated classes:
* ``primal::Point`` --  represents the coordinates of a point in space
* ``primal::Triangle`` -- represents a triangle as a triplet of ``primal::Point``s
* ``primal::BoundingBox`` -- represents an axis-aligned slab bounding a region of space

Since computational geometry primitives are typically used in performance-critical code sections, ``primal``'s primitives and operators are typically templated on the spatial dimension as well as the coordinate types. Thus, e.g., we have ``primal::Point<double, 3>`` for a point with 64-bit floating point coordinates in three-dimensional space. 

To simplify code readability, we often use type aliases when the types are known. For example, if our application always uses 3D triangles with floating point coordinate types, we might define:
```cpp
using BoundingBoxType = axom::primal::BoundingBox<double, 3>;
```

We will also be using the following operators:
* ``primal::intersect(BoundingBoxType, BoundingBoxType) -> bool`` -- checks whether a pair of bounding boxes intersect
* ``primal::intersect(TriangleType, TriangleType) -> bool`` -- checks whether a pair of triangles intersect
* ``primal::compute_bounding_box(TriangleType) -> BoundingBoxType`` -- computes and returns the axis-aligned bounding box for a given triangle

## Precomputing additional information about our triangle mesh
In this section, we extend our ``TriangleMesh`` class to include additional information about the triangles in our mesh. In particular, several of our algorithms require the bounding boxes of the triangles in the mesh, or of the entire mesh, so we'll precompute that as an ``axom::Array<BoundingBox>``:
```cpp
// in makeTriangleMesh(const std::string& stl_mesh_path, double weldThreshold) -> TriangleMesh 
  triMesh.m_triangleBoundingBoxes.reserve(numCells);
  for(const auto& tri : triMesh.triangles())
  {
    triMesh.m_triangleBoundingBoxes.emplace_back(
      axom::primal::compute_bounding_box(tri));
    triMesh.m_meshBoundingBox.addBox(triMesh.m_triangleBoundingBoxes.back());
  }

  SLIC_INFO(
    axom::fmt::format("Mesh bounding box is {}.", triMesh.meshBoundingBox()));
```

We expose the bounding boxes for the triangles and the overall mesh with the following helper functions:
```cpp
// New public functions in our local TriangleMesh class
  BoundingBox& meshBoundingBox();
  const BoundingBox& meshBoundingBox() const;

  axom::Array<BoundingBox>& triangleBoundingBoxes();
  const axom::Array<BoundingBox>& triangleBoundingBoxes() const;
```

Similarly, our algorithms will not apply to <i>degenerate</i> triangles, i.e. triangles with 0 area. As such, we precompute an array of booleans, using ``primal::Triangle::degenerate()`` to save later computation and expose these as:
```cpp
// More new public functions in our local TriangleMesh class
  axom::IndexType numDegenerateTriangles() const;
  bool isTriangleDegenerate(axom::IndexType idx) const;
```

> :clapper: Show updated ``TriangleMesh`` class and ``makeTriangleMesh()`` helper function

## A first self-intersection algorithm
We are now ready to present our first self-intersection algorithm!

We will use a nested ``for`` loop and check every pair of triangles for intersections. To prepare for future lessons, we will pass the intersection predicate in as a lambda, which uses ``primal::intersect(Triangle, Triangle) -> bool``. We will then return the results as an array of index pairs: ``axom::Array< std::pair<axom::IndexType, axom::IndexType>>``.

The basic structure of our algorithm is thus:
```cpp
using IndexPair = std::pair<axom::IndexType, axom::IndexType>;

template <typename IntersectionLambda>
axom::Array<IndexPair> naiveFindIntersections(const TriangleMesh& triMesh,
                                              IntersectionLambda&& trianglesIntersect,
                                              bool verboseOutput = false)
{
  axom::Array<IndexPair> intersectionPairs;

  const int numTriangles = triMesh.numTriangles();
  for(axom::IndexType idx1 = 0; idx1 < numTriangles - 1; ++idx1)
  {
    for(axom::IndexType idx2 = idx1 + 1; idx2 < numTriangles; ++idx2)
    {
      if(trianglesIntersect(idx1, idx2))
      {
        intersectionPairs.emplace_back(std::make_pair(idx1, idx2));
      }
    }
  }

  return intersectionPairs;
}
```
where the ``trianglesIntersect`` lambda is defined as:
```cpp
auto checkIntersect = [=,
                       tol = params.intersectionThreshold,
                       &mesh](axom::IndexType idx1, axom::IndexType idx2) {
  constexpr bool includeBoundaries = false;  // only use triangle interiors
  const auto& tris = mesh.triangles();

  return axom::primal::intersect(tris[idx1], tris[idx2], includeBoundaries, tol);
};
```
and the function is called as:
```cpp
auto intersectionPairs = naiveFindIntersections(mesh, checkIntersect, params.isVerbose());
```

> :information_source: Our ``intersect(Triangle, Triangle)`` operation is modeled after the state-of-the-art ["Faster Triangle-Triangle Intersection"](https://hal.inria.fr/inria-00072100/) algorithm of Olivier Devillers and Phillipe Guigue. It uses a decision-tree based on orientation primitives for its intersection predicates.

> :memo:  Since triangle intersection is a symmetric operation, the inner loop of our algorithm over ``idx2`` is initialized to begin at ``idx1 + 1`` (rather than ``0``) to avoid redundant tests.

> :information_source:  The above algorithm has complexity ``O(n^2)`` since every triangle must be tested against every other triangle.

### Wrinkle: ``primal::intersect`` only works for non-degenerate triangles
We need to modify our algorithm to skip degenerate triangles. Since we already precomputed this, we can modify our ``naiveFindIntersections`` algorithm as follows:
```cpp
template <typename IntersectionLambda>
axom::Array<IndexPair> naiveFindIntersections(const TriangleMesh& triMesh,
                                              IntersectionLambda&& trianglesIntersect,
                                              bool verboseOutput = false)
{
  axom::Array<IndexPair> intersectionPairs;

  const int numTriangles = triMesh.numTriangles();

  for(axom::IndexType idx1 = 0; idx1 < numTriangles - 1; ++idx1)
  {
    if(triMesh.isTriangleDegenerate(idx1))  // <-- HERE
    {
      continue;
    }

    SLIC_INFO_IF(
      verboseOutput && idx1 % 100 == 0,
      axom::fmt::format(axom::utilities::locale(), "Outer index {:L}", idx1));

    for(axom::IndexType idx2 = idx1 + 1; idx2 < numTriangles; ++idx2)
    {
      if(triMesh.isTriangleDegenerate(idx2)) // <-- HERE
      {
        continue;
      }

      if(trianglesIntersect(idx1, idx2))
      {
        intersectionPairs.emplace_back(std::make_pair(idx1, idx2));
      }
    }
  }

  return intersectionPairs;
}
```

We also took the opportunity to add some additional output when the user requests verbose output. The following snippet outputs a message in the outer loop every 100 iterations using the ``SLIC_INFO_IF`` macro:
```cpp
    SLIC_INFO_IF(
      verboseOutput && idx1 % 100 == 0,
      axom::fmt::format(axom::utilities::locale(), "Outer index {:L}", idx1));
```

## Optimization: Use the triangle bounding boxes to accelerate the query
We observe that the above "naive" algorithm can be quite inefficient since it must check every pair of triangles. We will address this in subsequent lessons using a spatial index to limit which triangles must be compared against each other. 

For now, we'll attempt to accelerate the local intersection tests by using a less expensive proxy for intersection -- checking whether the triangle bounding boxes intersect -- before running the more expensive triangle-based intersection tests. This will not improve the algorithmic complexity, but will (hopefully) reduce the runtime of the algorithm in practice.

We add a command line option ``--use-bounding-boxes`` and define a new lambda based on ``primal::intersect(BoundingBoxType, BoundingBoxType)`` to use it:
```cpp
  auto checkIntersectWithBoundingBoxes =
    [=, tol = params.intersectionThreshold, &mesh](axom::IndexType idx1,
                                                   axom::IndexType idx2) {
      constexpr bool includeBoundaries = false;  // only use triangle interiors
      const auto& tris = mesh.triangles();
      const auto& bboxes = mesh.triangleBoundingBoxes();

      return axom::primal::intersect(bboxes[idx1], bboxes[idx2]) &&
        axom::primal::intersect(tris[idx1], tris[idx2], includeBoundaries, tol);
    };
```
We can call our ``naiveFindIntersections()`` algorithm with this new lambda without any additional changes:
```cpp
// in main()
  auto intersectionPairs = params.useBoundingBoxes
    ? naiveFindIntersections(mesh,
                             checkIntersectWithBoundingBoxes,
                             params.isVerbose())
    : naiveFindIntersections(mesh, checkIntersect, params.isVerbose());
```

> :clapper: Look at code for completed example

## Running the code
So, how did our algorithms do?

Let's run the code (configured in "Release" mode) on a mesh with a few known self-intersections and compare the runtimes.

```shell
>./bin/lesson_02_naive_self_intersections -i ../stl_meshes/plane_simp_problems.stl 

[lesson_02: INFO] 
     Parsed parameters:
      * STL mesh: '../stl_meshes/plane_simp_problems.stl'
      * Threshold for welding: 1e-06
      * Skip welding: false
      * Threshold for intersections: 1e-08
      * verbose logging: false
      * use bounding boxes to accelerate query: false
       
[lesson_02: INFO] Reading file: '../stl_meshes/plane_simp_problems.stl'...
 
[lesson_02: INFO] Loading the mesh took 0.0704 seconds. 
[lesson_02: INFO] Vertex welding took 0.0131 seconds. 
[lesson_02: INFO] After welding, mesh has 15,006 vertices and 29,998 triangles. 
[lesson_02: INFO] Mesh has 0 degenerate triangles. 
[lesson_02: INFO] Mesh bounding box is { min:(-419.473,-519.348,-33.8787); max:(419.575,845.692,198.809); range:<839.047,1365.04,232.687> }. 
[lesson_02: INFO] Computing intersections without bounding boxes took 8.43 seconds. 
[lesson_02: INFO] Mesh had 5 intersection pairs 
```

```shell
>./bin/lesson_02_naive_self_intersections -i ../stl_meshes/plane_simp_problems.stl --use-bounding-boxes

[lesson_02: INFO] 
     Parsed parameters:
      * STL mesh: '../stl_meshes/plane_simp_problems.stl'
      * Threshold for welding: 1e-06
      * Skip welding: false
      * Threshold for intersections: 1e-08
      * verbose logging: false
      * use bounding boxes to accelerate query: true
       
[lesson_02: INFO] Reading file: '../stl_meshes/plane_simp_problems.stl'...
 
[lesson_02: INFO] Loading the mesh took 0.0699 seconds. 
[lesson_02: INFO] Vertex welding took 0.0131 seconds. 
[lesson_02: INFO] After welding, mesh has 15,006 vertices and 29,998 triangles. 
[lesson_02: INFO] Mesh has 0 degenerate triangles. 
[lesson_02: INFO] Mesh bounding box is { min:(-419.473,-519.348,-33.8787); max:(419.575,845.692,198.809); range:<839.047,1365.04,232.687> }. 
[lesson_02: INFO] Computing intersections with bounding boxes took 0.69 seconds. 
[lesson_02: INFO] Mesh had 5 intersection pairs 
```

So, for a mesh with around 30,000 triangles, the default intersection query took ``8.43`` seconds when not using bounding boxes. This was reduced to around ``.69`` seconds when using bounding boxes.

When we try this with a mesh with an order of magnitude more triangles (~250K triangles), our optimized "naive" algorithm begins to take unreasonable amounts of time (>45 seconds in this case):
```shell
>./bin/lesson_02_naive_self_intersections -i ../stl_meshes/golfBall.stl  --use-bounding-boxes
...
[lesson_02: INFO] Loading the mesh took 0.759 seconds. 
[lesson_02: INFO] Vertex welding took 0.121 seconds. 
[lesson_02: INFO] After welding, mesh has 122,882 vertices and 245,760 triangles. 
[lesson_02: INFO] Mesh has 0 degenerate triangles. 
[lesson_02: INFO] Mesh bounding box is { min:(-0.499487,-0.499529,-0.499867); max:(0.499495,0.499529,0.499867); range:<0.998982,0.999059,0.999735> }. 
[lesson_02: INFO] Computing intersections with bounding boxes took 46.5 seconds. 
[lesson_02: INFO] Mesh had 0 intersection pairs 
```

### Next time:
In the next lesson, we will port our naive algorithm to other execution and memory spaces with the help of the [RAJA](https://github.com/LLNL/raja) and [Umpire](https://github.com/LLNL/umpire) libraries.