---
name: cpp-unit-testing
description: Instructions for creating unit tests for Axom
---

This repo’s unit tests:
- Use GoogleTest (`TEST`, `TEST_F`, etc.).
- Live close to the code they cover in a `tests/` subdirectory.
- Are registered with CTest via the CMake macro `axom_add_test(...)`.

## Where tests go

Place unit tests next to the component they cover:
- `src/axom/<component>/tests/` (recommended for library code)
- `src/thirdparty/tests/` (third-party smoke tests)

Each `tests/` directory should contain:
- One or more `*.cpp` test files
- A `CMakeLists.txt` that builds test executables with `axom_add_executable` and registers them with `axom_add_test`

Also ensure the parent component adds the tests subdirectory when tests are enabled:

```cmake
if(AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()
```

## Writing GoogleTest unit tests

- Include GoogleTest: `#include <gtest/gtest.h>`
- Prefer fast, deterministic, CPU-only tests by default.
- Use fixtures (`TEST_F`) when setup/teardown is shared.
- Avoid depending on external files unless the test is explicitly a data-driven/integration test.

## Registering tests with `axom_add_test`

Axom tests typically create one test executable per source file with `axom_add_executable`, then register it with CTest using `axom_add_test`.
The executable name usually appends `_test`; the CTest name usually omits that suffix.

Minimal `tests/CMakeLists.txt` template:

```cmake
set(my_component_test_depends
  my_component
  gtest)

set(my_component_test_sources
  feature_a.cpp
  feature_b.cpp)

foreach(test ${my_component_test_sources})
  get_filename_component(test_name ${test} NAME_WE)
  axom_add_executable(NAME       ${test_name}_test
                      SOURCES    ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${my_component_test_depends}
                      FOLDER     axom/my_component/tests)
  axom_add_test(NAME    ${test_name}
                COMMAND ${test_name}_test)
endforeach()
```

Notes:
- Keep test source basenames unique across the project (CMake targets are global).
- Add only the libraries you need to `DEPENDS_ON` (the library under test + `gtest` + any required Axom/TPL targets).

### MPI / OpenMP tests

If a test requires MPI, add the MPI dependency and set `NUM_MPI_TASKS`:

```cmake
axom_add_executable(NAME       my_parallel_test
                    SOURCES    my_parallel_test.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON gtest blt::mpi my_component
                    FOLDER     axom/my_component/tests)
axom_add_test(NAME          my_parallel
              COMMAND       my_parallel_test
              NUM_MPI_TASKS 4)
```

If a test requires OpenMP, set `NUM_OMP_THREADS` accordingly.

### GPU-enabled unit tests

If you need a CUDA- or HIP-enabled unit test, gate it behind the matching Axom option and add the matching BLT GPU dependency:

- CUDA: `AXOM_ENABLE_CUDA` with `blt::cuda` for CUDA-compiled sources, or the CUDA runtime target when only runtime API linkage is needed.
- HIP: `AXOM_ENABLE_HIP` with `blt::hip`.

```cmake
if(AXOM_ENABLE_HIP)
  axom_add_executable(NAME       my_hip_test
                      SOURCES    my_hip_test.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON gtest blt::hip my_component
                      FOLDER     axom/my_component/tests)
  axom_add_test(NAME    my_hip
                COMMAND my_hip_test)
endif()
```

Use the same shape for CUDA tests:

```cmake
if(AXOM_ENABLE_CUDA)
  axom_add_executable(NAME       my_cuda_test
                      SOURCES    my_cuda_test.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON gtest blt::cuda my_component
                      FOLDER     axom/my_component/tests)
  axom_add_test(NAME    my_cuda
                COMMAND my_cuda_test)
endif()
```

For tests that are always built but need GPU linkage in GPU builds, append the dependency conditionally before creating the target:

```cmake
set(my_component_test_depends
  my_component
  gtest)
blt_list_append(TO my_component_test_depends ELEMENTS blt::cuda IF AXOM_ENABLE_CUDA)
blt_list_append(TO my_component_test_depends ELEMENTS blt::hip IF AXOM_ENABLE_HIP)
```

Use the repository's existing local convention when modifying an existing `CMakeLists.txt`: some test directories use `ENABLE_HIP` for build-configuration workarounds such as `axom_force_release_for_target(...)`, while dependency lists commonly use `AXOM_ENABLE_HIP`.

## Running tests

- Configure with tests enabled: `-DENABLE_TESTS=ON -DAXOM_ENABLE_TESTS=ON`
- Run all tests: `ctest`
- Run a single test by name: `ctest -R test_state_manager`

## Examples in this repo

- Macro definition: `src/cmake/AxomMacros.cmake`
- Typical component tests: `src/axom/inlet/tests/CMakeLists.txt`
- Third-party smoke tests: `src/thirdparty/tests/CMakeLists.txt`
