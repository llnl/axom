# Lesson 01: Load an STL mesh

In this lesson, we will begin developing our application by loading a triangle mesh in the STL format and printing out some of its properties.

We use the [CLI11](https://github.com/CLIUtils/CLI11) library for parsing command line options. Axom vendors a copy of ``CLI11`` for internal use and provides it to downstream users. It is available via the `axom::cli11` CMake target.

This example uses the following Axom components:
* ``core`` -- has a dynamic array class ``axom::Array``
* ``slic`` -- provides logging functionality
* ``quest`` -- has a utility function for loading the STL mesh
* ``mint`` -- provides a class for triangle meshes

## Using a SLIC-based logger

We begin by setting up a logger for our application using Axom's `slic` component ([Simple Logging Interface Code](https://axom.readthedocs.io/en/develop/axom/slic/docs/sphinx/index.html)).

`slic` uses a C-style API over a singleton instance and also provides a few basic logging macros that applications can either use or extend/build upon. 

For our application, we will use an RAII-based wrapper struct:

```cpp
struct BasicLogger
{
  BasicLogger()
  {
    namespace slic = axom::slic;

    // Initialize the SLIC logger
    slic::initialize();
    slic::setLoggingMsgLevel(slic::message::Debug);

    // Customize logging levels and formatting
    const std::string slicFormatStr = "[lesson_01: <LEVEL>] <MESSAGE> \n";

    slic::addStreamToMsgLevel(new slic::GenericOutputStream(&std::cerr),
                              slic::message::Error);
    slic::addStreamToMsgLevel(
      new slic::GenericOutputStream(&std::cerr, slicFormatStr),
      slic::message::Warning);

    auto* compactStream =
      new slic::GenericOutputStream(&std::cout, slicFormatStr);
    slic::addStreamToMsgLevel(compactStream, slic::message::Info);
    slic::addStreamToMsgLevel(compactStream, slic::message::Debug);
  }

  ~BasicLogger() { axom::slic::finalize(); }
};
```

The constructor calls ``slic::initialize()``, sets the logging level to ``slic::message::Debug`` and then creates log streams for the different message levels. The destructor cleans up SLIC's memory with `slic::finalize()`.

> :memo: Axom provides a similar basic wrapper in the ``slic::SimpleLogger`` class.

> :bulb: If your application requires MPI-based logging, consider using the ``LumberjackStream`` log streams from Axom's [lumberjack component](https://axom.readthedocs.io/en/develop/axom/lumberjack/docs/sphinx/index.html) when setting up your logger.

> :memo:  `slic` has an internal formatter that understands the following keywords: `<LEVEL>`, `<MESSAGE>`, `<TAG>`, `<FILE>`, `<LINE>`, `<RANK>`, `<RANK_COUNT>` and `<TIMESTAMP>`. Our `BasicLogger` only uses `<LEVEL>` and `<MESSAGE>`.

We now add a `BasicLogger` instance to the application with
```cpp
  // Initialize logger; use RAII so it will finalize at the end of the application
  BasicLogger logger;
```

`slic`'s logging levels are:
* ``slic::message::Error``: Used for logging error messages. By default prints a stacktrace and exits the application.
* ``slic::message::Warning``: Used for logging warning messages. By default, this prints a warning message and continues executing the application.
* ``slic::message::Info``: Used for logging informational messages.
* ``slic::message::Debug``: Used for logging debug messages. By default, these are compiled out in ``Release`` configurations.
 
We are now free to use `slic`-based logging macros throughout our application, for example:
```cpp
  // Use slic macros to log some messages
  SLIC_INFO("A first message!");
  SLIC_WARNING("A warning message!");
  SLIC_DEBUG("A debug message!");
```
These produce the following output:
```shell
[lesson_01: INFO] A first message! 
[lesson_01: WARNING] A warning message! 
[lesson_01: DEBUG] A debug message! 
```

## Parsing command line input with CLI11

Axom provides the [CLI11](https://github.com/CLIUtils/CLI11) library to downstream users.

We can use this by including the `CLI11` header
```cpp
#include "axom/CLI11.hpp"
```

Our application defines a struct `Input` for parsing our input and for holding the results. 
For this lesson, our struct will hold a path to a triangle mesh file and a flag to indicate if we want to print verbose output:
```cpp
struct Input
{
  std::string mesh_file {""};
  bool verboseOutput {false};

  void parse(int argc, char** argv, axom::CLI::App& app);
  bool isVerbose() const { return verboseOutput; }
};
```

The bulk of this struct is in the implementation of the `parse()` function:
```cpp
void Input::parse(int argc, char** argv, axom::CLI::App& app)
{
  app.add_option("-i, --infile", mesh_file)
    ->description("The input STL mesh file")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_flag("-v,--verbose", verboseOutput)
    ->description("Increase logging verbosity")
    ->capture_default_str();

  app.get_formatter()->column_width(40);

  app.parse(argc, argv);  // Could throw an exception

  // Output parsed information
  SLIC_INFO(axom::fmt::format(R"(
     Parsed parameters:
      * STL mesh: '{}'
      * verbose logging: {}
      )",
                              mesh_file,
                              verboseOutput));
}
```
We will add parameters to the `Input` struct and logic for parsing to `Input::parse()` as necessary throughout that tutorial.

> :memo:  For the input file, `CLI11`'s `check(axom::CLI::ExistingFile)` ensures that the provided file exists, and its `required()` function for options ensures that the user provides a value.

We next add the `CLI::App` to our application's main function:
```cpp
  // Parse the command line arguments
  Input params;
  {
    axom::CLI::App app {"STL mesh loader"};
    try
    {
      params.parse(argc, argv, app);
    }
    catch(const axom::CLI::ParseError& e)
    {
      return app.exit(e);
    }
  }
```

## Setting up our triangle mesh

Our final step in this lesson will be to load the triangle mesh. Axom has built-in support for loading STL meshes via the ``quest::STLReader`` class, which loads the mesh into a ``mint::UnstructuredMesh`` instance. For the purposes of this tutorial, we will begin with a (simpler) array of triangles as our basic mesh data structure. We abstract this into our local ``TriangleMesh`` struct:
```cpp
struct TriangleMesh
{
  using Point = axom::primal::Point<double, 3>;
  using Triangle = axom::primal::Triangle<double, 3>;

  axom::IndexType numTriangles() const { return m_triangles.size(); }
  axom::Array<Triangle>& triangles() { return m_triangles; }
  const axom::Array<Triangle>& triangles() const { return m_triangles; }

  axom::Array<Triangle> m_triangles;
};
```
This class uses three internal types:

  ``axom::Array``
  : A drop-in replacement for ``std::vector`` which supports usage on host (i.e. CPU) and device (i.e. GPU) execution spaces

  ``axom::primal::Point``
  : A geometric primitive representing a point in space

  ``axom::primal::Triangle``
  : A geometric primitive representing a triangle defined by three ``axom::primal::Point`` instances

We will discuss these in more depth in subsequent lessons.

We use a helper function to load the triangle mesh:
```cpp
TriangleMesh makeTriangleMesh(const std::string& stl_mesh_path)
{
  TriangleMesh triMesh;

  // load STL mesh into a mint unstructured mesh
  auto surface_mesh =
    std::make_unique<axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>>(
      3,
      axom::mint::TRIANGLE);
  {
    auto reader = std::make_unique<axom::quest::STLReader>();
    reader->setFileName(stl_mesh_path);
    reader->read();
    reader->getMesh(surface_mesh.get());
  }

  // extract triangles into an axom::Array
  const int numCells = surface_mesh->getNumberOfCells();
  triMesh.m_triangles.reserve(numCells);
  for(int i = 0; i < numCells; ++i)
  {
    std::array<axom::IndexType, 3> triCell;
    TriangleMesh::Triangle tri;

    surface_mesh->getCellNodeIDs(i, triCell.data());
    surface_mesh->getNode(triCell[0], tri[0].data());
    surface_mesh->getNode(triCell[1], tri[1].data());
    surface_mesh->getNode(triCell[2], tri[2].data());

    triMesh.m_triangles.emplace_back(tri);
  }

  return triMesh;
}
```
The above function loads an STL mesh into a ``mint::UnstructuredMesh``, extracts its triangles into an ``axom::Array`` and returns the initialized ``TriangleMesh`` instance.

We call this from ``main()`` using:
```cpp
  TriangleMesh mesh = makeTriangleMesh(params.mesh_file);
```

Finally, we can output some properties of the parsed mesh using:
```cpp
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Parsed STL mesh has {:L} triangles.",
                              mesh.numTriangles()));
```
> :bulb: The ``locale`` helps with the formatting of the output for improved readability. Axom's default locale is `"en_US.UTF-8"`.

## Running the application
The [axom_data](https://github.com/LLNL/axom_data) repository provides several STL meshes that we will use in this tutorial. We will also use several larger meshes that are not part of the ``axom_data`` repo.

> :memo: Axom uses `axom_data` as a submodule in `<axom_root>/data`

#### Run application without passing an input mesh
```shell
# Missing path to STL mesh
>./bin/lesson_01_load_stl_mesh

[lesson_01: INFO] A first message! 
[lesson_01: WARNING] A warning message! 
[lesson_01: DEBUG] A debug message! 

--infile is required
Run with --help for more information.
```

#### Print usage with `-h`
```shell
# Print usage info
>./bin/lesson_01_load_stl_mesh -h

[lesson_01: INFO] A first message! 
[lesson_01: WARNING] A warning message! 
[lesson_01: DEBUG] A debug message! 
STL mesh loader
Usage: ./bin/lesson_01_load_stl_mesh [OPTIONS]

Options:
  -h,--help                             Print this help message and exit
  -i,--infile TEXT:FILE REQUIRED        The input STL mesh file
  -v,--verbose                          Increase logging verbosity
```

#### Run application with valid STL mesh
```shell
# Run with valid path to STL mesh
>./bin/lesson_01_load_stl_mesh -i ../stl_meshes/sphere.stl 

[lesson_01: INFO] A first message! 
[lesson_01: WARNING] A warning message! 
[lesson_01: DEBUG] A debug message! 

[lesson_01: INFO] 
     Parsed parameters:
      * STL mesh: '../stl_meshes/sphere.stl'
      * verbose logging: false
       
[lesson_01: INFO] Reading file: '../stl_meshes/sphere.stl'...
 
[lesson_01: INFO] Parsed STL mesh has 1,280 triangles. 
```

### Next time:
In the next lesson, we will check the mesh triangles for self-intersections.
