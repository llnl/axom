# Lesson 03: Naive self-intersections... now on the device!

In this lesson, we will port our naive self-intersection algorithm from the previous lesson to different execution spaces, including serial (sequential), threaded (OpenMP) and/or GPU-based (CUDA) backends with the help of the RAJA and Umpire libraries. These libraries allow a single implementation to run on different execution spaces.

> :memo:  We are still using (an adaptation of) our naive algorithm in this lesson, so the performance will improve but will still scale quadratically in the number of triangles in the mesh.

> :warning:  This lesson requires Axom to be configured with RAJA and Umpire support. It will also strongly benefit from OpenMP and/or GPU support.

## A note about our thread-safe implementation

We begin with the naive self-intersection algorithm that we developed in the previous lesson.

To support a thread-safe implementation of our algorithm, our approach will be to split our algorithm into two passes: 
* In the first pass, we will count the number of intersections. 
* If we find any intersections, we allocate an array with sufficient storage to hold all intersection pairs. 
* Our second pass will find all intersections and insert their indices into a buffer.

We note that this approach essentially duplicates the work when there are intersections. This could be improved in the average case, e.g. by estimating the number of expected intersections and inserting up to that many in a first pass, and only running a second pass when necessary. We will, however, not pursue this approach. Rather, in the next lesson, we will develop a better algorithm based on spatial indexes.

## Adding execution policies

We begin by exposing several execution policies in the application.

We use Axom's built-in execution policies from the ``core`` component.

Our application will require RAJA, and can optionally be build against OpenMP and/or CUDA, so we'll guard against the availability of the latter two:

```cpp
#include "axom/config.hpp"
#include "axom/core.hpp"

using seq_exec = axom::SEQ_EXEC;

#if defined(AXOM_USE_OPENMP)
  using omp_exec = axom::OMP_EXEC;
#else
  using omp_exec = seq_exec;
#endif

#if defined(AXOM_USE_CUDA)
  constexpr int BLK_SZ = 256;
  using cuda_exec = axom::CUDA_EXEC<BLK_SZ>;
#else
  using cuda_exec = seq_exec;
#endif
```

We also add user support via the following ``RuntimePolicy`` enum:
```cpp
using RuntimePolicy = axom::runtime_policy::Policy;
```
which contains the following:

```cpp
namespace axom {
namespace runtime_policy {
  enum class Policy
  {
    seq  = 0,   // always available
    omp  = 1,   // only available when Axom is configured with OpenMP
    cuda = 2,   // only available when Axom is configured with cuda
    hip  = 3    // only available when Axom is configured with hip
  };
}
}
```
and add a ``RuntimePolicy`` member to  our ``Input`` class:
```cpp
  RuntimePolicy policy {RuntimePolicy::seq};
```
and use ``CLI11``'s ``CheckedTransformer`` to ensure that the user can only input valid choices for our configuration.

We expose this as a new  command-line argument ``--policy`` (with shortcut ``-p``).

> :clapper: Show this in the code if the audience is interested in the details

## Calling our find intersections algorithm

Since the execution space is templated, and we would to expose a runtime choice to our application users, we add a switch statement to ``main()`` to call the desired implementation:

```cpp
  // Check for self-intersections; results are returned as an array of index pairs
  axom::Array<IndexPair> intersectionPairs;

  switch(params.policy)
  {
  case RuntimePolicy::omp:
#ifdef AXOM_USE_OPENMP
    intersectionPairs =
      naiveFindIntersections<omp_exec>(mesh,
                                       params.intersectionThreshold,
                                       params.useBoundingBoxes,
                                       params.isVerbose());
#endif
    break;
  case RuntimePolicy::cuda:
#ifdef AXOM_USE_CUDA
    intersectionPairs =
      naiveFindIntersections<cuda_exec>(mesh,
                                        params.intersectionThreshold,
                                        params.useBoundingBoxes,
                                        params.isVerbose());
#endif
    break;
  default:  // RuntimePolicy::seq
    intersectionPairs =
      naiveFindIntersections<seq_exec>(mesh,
                                       params.intersectionThreshold,
                                       params.useBoundingBoxes,
                                       params.isVerbose());
    break;
  }
```
> :memo:  For simplicity, we moved our ``checkIntersection`` lambda inside the ``naiveFindIntersections()`` algorithm.

> :bulb:  This interplay between compile-time templated code used in our internal implementations and runtime-selected parameters that we expose to users is a common theme in Axom development and usage.

## Algorithm setup: memory spaces

Our algorithm requires memory to be allocated in the correct space for each kernel. As such, we distinguish between "host" and "kernel" memory allocators. The latter might be a ``host`` memory allocator, e.g. for sequential and OpenMP execution spaces. It might also be ``unified`` or ``device`` memory if run on GPUs.

We get the integer ID for the ``host_allocator`` and ``kernel_allocator`` as follows:
```cpp
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

  // Get ids of necessary allocators
  const int host_allocator =
    axom::getUmpireResourceAllocatorID(umpire::resource::Host);

  const int kernel_allocator = on_device
    ? axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    : axom::execution_space<ExecSpace>::allocatorID();

```
Note that the above calls utilize the function's ``ExecSpace`` template parameter as well as Axom's ``execution_space`` traits classes.

When allocating an ``axom::Array``, we pass this allocator ID to the appropriate constructor to allocate or copy data.

## Algorithm setup: creating a list of non-degenerate triangles.

We begin our algorithm with a slightly different approach to handling degenerate triangles than the previous example. Rather than checking for degeneracies as we iterate through the triangle pairs, we create a new array holding the indices of the valid (i.e. non-degenerate) triangles and copy the results to the device, if necessary: 
```cpp
  const auto totalTriangles = triMesh.numTriangles();
  IndexArray valid_h(0, totalTriangles, host_allocator);
  for(axom::IndexType idx = 0; idx < totalTriangles; ++idx)
  {
    if(!triMesh.isTriangleDegenerate(idx))
    {
      valid_h.push_back(idx);
    }
  }
  auto valid_d = on_device ? IndexArray(valid_h, kernel_allocator) : IndexArray();
  auto valid_v = on_device ? valid_d.view() : valid_h.view();
```

> :warning: We assume that this is a relatively inexpensive operation, so we perform this serially on the CPU for now and copy the data over to the device, as necessary. If/When this assumption changes, we can revisit that decision.


## Aside: ``axom::Array`` and ``axom::ArrayView``

We are now ready to discuss the ``axom::Array`` class, and its counterpart, ``axom::ArrayView``. 

``axom::Array`` is a drop-in replacement for ``std::vector`` which can be used with device code. In contrast to ``std::vector``, it:
* supports Umpire allocators for its memory allocations
* has functions that are host/device decorated for execution within GPU kernels
* does not require its members to be initialized on copy. This can be useful, e.g. when allocating memory that will be initialized immediately in a kernel.
* supports compile-time (non-ragged) multidimensional indexing <i>*[Not demonstrated in this tutorial]</i>
* has a `view()` function that returns an ``axom::ArrayView`` over the same memory

``axom::Array`` is a dynamic array class, and since (a) we do not (currently) support memory allocations within our GPU kernels, and (b) our host/device decorated kernels can only capture variables by value, we cannot use ``axom::Array`` within our RAJA kernels. 

Rather, we pass in a lightweight view class ``axom::ArrayView`` -- similar to a ``std::span`` -- that contains a pointer to the underlying data and the number of elements within an ``axom::Array``. ArrayView supports a subset of the operations of ``axom::Array``; it allows users to modify the data in the array, but does not allow modifications that might change the size of the array.

##### Naming conventions
Our implementation will use the following variable naming convention for arrays, views and pointers:
* We use the suffix ``_h`` for variables that live in the host space, e.g. ``foo_h``
* We use the suffix ``_d`` for variables that live in the device space, e.g. ``foo_d``
* We use the suffix ``_v`` for views of an array, e.g. ``foo_v``
* We use the suffix ``_p`` for raw pointers, e.g. ``foo_p``

## First pass: count the number of intersections

Now that we know the number of valid triangles and their indices, we are ready to apply the first pass of our algorithm, which counts the number of pairs of intersecting triangles.

We perform this test in parallel and use a reduction variable to atomically increment our intersection count:

```cpp
// Phase I: Find the number of intersecting pairs of triangles
  const auto validCount = valid_v.size();
  RAJA::RangeSegment row_range(0, validCount);
  RAJA::RangeSegment col_range(0, validCount);

  using KERNEL_POL = typename axom::internal::nested_for_exec<ExecSpace>::loop2d_policy;

  axom::ReduceSum<ExecSpace, int> numIntersect(0);

  RAJA::kernel<KERNEL_POL>(
    RAJA::make_tuple(col_range, row_range),
    AXOM_LAMBDA(int col, int row) {
      if(row < col && trianglesIntersect(valid_v[row], valid_v[col]))
      {
        numIntersect += 1;
      }
    });
```
> :memo: ``numIntersect`` is a reduction variable that sums the results over each kernel iteration. We access its value outside the kernel via the ``.get()`` member function.

The ``trianglesIntersect`` lambda is similar to the one in our previous lesson, although it now uses views into the triangle array and the bounding box array, and its caller is expected to only pass in valid (i.e. non-degenerate) triangle indices:
```cpp
  // lambda to check if two triangles w/ given indices intersect
  auto trianglesIntersect = AXOM_LAMBDA(axom::IndexType idx1, axom::IndexType idx2)
  {
    constexpr bool includeBoundaries = false;  // only use triangle interiors

    return (!useBoundingBoxes || axom::primal::intersect(bbox_v[idx1], bbox_v[idx2])) &&
      axom::primal::intersect(tris_v[idx1], tris_v[idx2], includeBoundaries, tol);
  };
```

## Second pass: storing the indices of the intersections

If the first pass did not find any intersections, we can stop execution here. Otherwise, we need to take another pass through the data to extract the indices of the intersecting pairs of triangles.

We first allocate an array (in the proper space), without initializing its members.
```cpp
    const auto numIndices = 2 * numIntersect.get();
    IndexArray intersections_d(axom::ArrayOptions::Uninitialized {},
                               numIndices,
                               kernel_allocator);
```

We use the counter, along with an ``axom::atomicAdd`` to keep track of the insertion location in the following algorithm:
```cpp
    auto intersections_v = intersections_d.view();
    RAJA::kernel<KERNEL_POL>(
      RAJA::make_tuple(col_range, row_range),
      AXOM_LAMBDA(int col, int row) {
        if(row < col && trianglesIntersect(valid_v[row], valid_v[col]))
        {
          const auto idx = axom::atomicAdd<axom::auto_atomic>(counter_p, axom::IndexType {2});
          intersections_v[idx + 0] = valid_v[row];
          intersections_v[idx + 1] = valid_v[col];
        }
      });
```
Our algorithm concludes by copying the intersection indices into an array of pairs of integers, while accounting for the former living in device memory space:
```cpp
    // copy intersections back to host, if necessary
    auto intersections_h =
      on_device ? IndexArray(intersections_d, host_allocator) : IndexArray();
    
    // in either event, intersections_h_v will be a valid host array
    auto interections_h_v =
      on_device ? intersections_h.view() : intersections_d.view();

    // copy the results into the return vector
    for(axom::IndexType idx = 0; idx < numIndices; idx += 2)
    {
      intersectionPairs.emplace_back(
        std::make_pair(interections_h_v[idx], interections_h_v[idx + 1]));
    }
```

## Comparing the performance
Let's see how our RAJA port performs.

#### Comparing the sequential RAJA-fied algorithm against our previous one
We first compare the performance of the sequential algorithm to our previous algorithm on a mesh with 30K triangles of which 5 pairs of triangles intersect:
```shell
# results from previous example
> ./bin/lesson_02_naive_self_intersections -i plane_simp_problems.stl --use-bounding-boxes

[lesson_02: INFO] Loading the mesh took 0.0706 seconds. 
[lesson_02: INFO] Vertex welding took 0.0126 seconds. 
[lesson_02: INFO] After welding, mesh has 15,006 vertices and 29,998 triangles. 
[lesson_02: INFO] Mesh has 0 degenerate triangles. 

[lesson_02: INFO] Computing intersections with bounding boxes took 0.712 seconds. 
[lesson_02: INFO] Mesh had 5 intersection pairs 
[lesson_02: INFO] Intersecting pairs:
	(488, 10975)
	(10973, 10975)
	(10975, 10998)
	(10975, 11000)
	(10975, 11001)
```

```shell
# results from current example
>./bin/lesson_03_device_self_intersections -i plane_simp_problems.stl --use-bounding-boxes -v
       
[lesson_03: INFO] Loading the mesh took 0.0694 seconds. 
[lesson_03: INFO] Vertex welding took 0.0129 seconds. 
[lesson_03: INFO] After welding, mesh has 15,006 vertices and 29,998 triangles. 
[lesson_03: INFO] Mesh has 0 degenerate triangles.

[lesson_03: INFO] Running naive intersection algorithm in execution Space: [SEQ_EXEC] 
[lesson_03: INFO] Computing non-degenerates took 6.48e-05 seconds. 
[lesson_03: INFO] Finding intersection count took 0.744 seconds. 
[lesson_03: INFO] Inserting intersection pairs took 0.701 seconds. 
[lesson_03: INFO] Copying back to array took 3.78e-06 seconds. 

[lesson_03: INFO] Computing intersections with bounding boxes took 1.44 seconds. 
[lesson_03: INFO] Mesh had 5 intersection pairs 
[lesson_03: INFO] Intersecting pairs:
	(488, 10975)
	(10973, 10975)
	(10975, 10998)
	(10975, 11000)
	(10975, 11001)
```
So, as expected, the sequential execution takes approximately twice as long to run on meshes that have self-intersections with this new formulation. However, we can now run the algorithm in parallel.

#### Comparing our algorithm in different execution spaces

##### Sequential vs. OpenMP on a Sapphire Rapids node
Here is an execution using 112 OpenMP threads on an Intel Sapphire Rapids node using a mesh with 370K triangles:
```shell
# sequential
> ./bin/lesson_03_device_self_intersections -i plane.stl --use-bounding-boxes -v --policy seq

[lesson_03: INFO] 
     Parsed parameters:
      * Runtime execution policy: 'seq'
...
[lesson_03: INFO] After welding, mesh has 186,708 vertices and 373,412 triangles. 
[lesson_03: INFO] Running naive intersection algorithm in execution Space: [SEQ_EXEC] 
[lesson_03: INFO] Computing intersections with bounding boxes took  145 seconds. 
[lesson_03: INFO] Mesh had 0 intersection pairs   

# threaded
> NUM_OMP_THREADS=112 ./bin/lesson_03_device_self_intersections -i plane.stl --use-bounding-boxes -v ---policy omp

[lesson_03: INFO] 
     Parsed parameters:
      * Runtime execution policy: 'omp'
       
[lesson_03: INFO] After welding, mesh has 186,708 vertices and 373,412 triangles. 
[lesson_03: INFO] Running naive intersection algorithm in execution Space: [OMP_EXEC] 
[lesson_03: INFO] Computing intersections with bounding boxes took 4.9 seconds. 
[lesson_03: INFO] Mesh had 0 intersection pairs 
```

##### Execution on a V100 GPU
Finally, here's a run on a single NVidia V100 GPU
```shell
>./bin/lesson_03_device_self_intersections -i plane.stl --use-bounding-boxes -v --policy cuda

[lesson_03: INFO] 
     Parsed parameters:
      * Runtime execution policy: 'cuda'
       
[lesson_03: INFO] After welding, mesh has 186,708 vertices and 373,412 triangles. 
[lesson_03: INFO] Running naive intersection algorithm in execution Space: [CUDA_EXEC] 
[lesson_03: INFO] Computing intersections with bounding boxes took 11.4 seconds. 
[lesson_03: INFO] Mesh had 0 intersection pairs 
```

> :clapper:  Let's run a few examples on different meshes

### Next time:
In this lesson, we ported our naive self-intersection algorithm to different execution space with the help of RAJA and Umpire. This accelerated our compute time, but did not help with the algorithmic complexity of the algorithm. In the next lesson, we will introduce a better algorithm based on a spatial index that, for reasonable meshes, should considerably reduce the expected work per triangle.
