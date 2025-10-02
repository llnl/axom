# Lesson 00: Details about our Axom configuration

In this lesson, we will develop a simple application that uses an installed version of Axom. 
Specifically, we will print the version of Axom and some of its configuration properties.

> :information_source: Our examples use the ``C++17`` standard since that's Axom's current version.

We will also be using the [fmt](https://fmt.dev/latest/index.html) library for string formatting, which Axom uses internally and vendors for its users.
The vast majority of ``fmt`` has been included in recent C++ standards, starting with ``C++20``. See also the [fmt cheat sheet](https://hackingcpp.com/cpp/libs/fmt.html).

## Configuring our application

We will be using CMake to configure our examples against a pre-installed copy of Axom configured with several third-party libraries (TPLs) including the following RADIUSS libraries:
* [Conduit](https://github.com/LLNL/conduit/)
* [MFEM](https://mfem.org/)
* [RAJA](https://github.com/LLNL/raja)
* [Umpire](https://github.com/LLNL/umpire)

Axom's build system is based on [BLT](https://github.com/LLNL/blt) -- a streamlined CMake-based foundation for <b>B</b>uilding, <b>L</b>inking and <b>T</b>esting large-scale high performance computing (HPC) applications.

Since Axom users typically develop HPC applications that need to run on several different platforms ("hosts") with several different compilers, we will be using a CMake cache file containing several configuration variables, including the compiler and MPI information as well as compilation flags and paths to third-party libraries. We refer to this file as a "host-config".

Our application can be configured and built using the following commands:
```shell
> cd <path/to/tutorial/files>
> mkdir build
> cd build
> cmake -C ../host-config.cmake ..
> make -j16
```

> :memo:  BLT provides a few lightweight "smoke" tests to help diagnose configuration issues. If we configure our application with the ``ENABLE_TESTS`` CMake option, the above command will build these tests and we can run them using ``ctest``.

The following CMake snippet imports the provided Axom installation into our project's build system:
```cmake
include(CMakeFindDependencyMacro)
find_dependency(axom REQUIRED
                NO_DEFAULT_PATH 
                PATHS ${AXOM_DIR}/lib/cmake)
```
It assumes that the user has provided the path to our Axom installation via the CMake variable ``AXOM_DIR`` and provides the top-level ``axom`` CMake target as well as several other targets in the ``axom::`` namespace.

This snippet generates an executable for our application using blt's [blt_add_executable](https://llnl-blt.readthedocs.io/en/develop/api/target.html#blt-add-executable) macro:
```cmake
blt_add_executable(NAME       lesson_00_check_axom_configuration 
                   SOURCES    lesson_00/check_axom_configuration.cpp
                   DEPENDS_ON axom axom::fmt)
```

> :clapper: Look at ``CMakeLists.txt`` at the root of the tutorial directory.


## Including Axom components

To simplify usage of Axom components in applications, Axom generates per-component include files listing the available header files for a given component.
For example, to include Axom's ``core`` component, an application would include the following line
```cpp
#include "axom/core.hpp"
```

Axom also generates a config header file containing most of the compiler defines related to the configured Axom. This file is located in ``axom/config.hpp``. For example, to guard MPI-related code, we would use the ``AXOM_USE_MPI`` as follows:
```cpp
#include "axom/config.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

...
```

> :clapper: Look at ``<axom/config.hpp>`` and ``<axom/core.hpp>`` in the installation directory

## Checking the Axom version

We will use ``axom::getVersion()`` to print the version of Axom that we're using. 
This function is provided in Axom's core component.
```cpp
#include "axom/core.hpp"

std::cout << axom::fmt::format("Version: {}", axom::getVersion()) << "\n\n";
```

This includes the major, minor and patch version of Axom as well as the git SHA (when available).
The above command might produce something like:
```
Version: v0.9.0-f5b5b5d66
```

The version is also available as compiler defines in ``axom/config.hpp``. The following corresponds to the above version:
```cpp
#define AXOM_VERSION_MAJOR 0
#define AXOM_VERSION_MINOR 9
#define AXOM_VERSION_PATCH 0
#define AXOM_VERSION_FULL  "v0.9.0"
```

## Checking details about the Axom configuration

Another useful function for checking features of Axom's configuration is ``axom::about()``, which is available in Axom's ``core`` component.
```cpp
axom::about();
```


## Running the application

After building the application, it will be located in the ``bin`` directory of our build system.
We can run the application using:
```shell
> ./bin/lesson_00_check_axom_configuration
```

Its output should look something like:
```
Axom information:
  AXOM_VERSION_FULL: v0.9.0
  AXOM_VERSION_MAJOR: 0
  AXOM_VERSION_MINOR: 9
  AXOM_VERSION_PATCH: 0
  AXOM_GIT_SHA: f5b5b5d66
Compiler Settings: 
  C++ Standard: c++17
  Size of axom::IndexType: 4
Active programming models: { mpi;openmp }
Available components: { core;inlet;klee;lumberjack;mint;primal;quest;sidre;slam;slic;spin }
Active built-in dependencies: { CLI11;fmt;sol;sparsehash }
Active external dependencies: { adiak;caliper;conduit;hdf5;lua;mfem;raja;umpire }
```

> :clapper: Run the example for this lesson


### Next time:
In the next lesson, we will add a logger and command line parser and load in a triangle mesh. 
