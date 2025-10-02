// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: load_stl_mesh.cpp
///
/// This example loads in an STL mesh and prints some properties.
/// It uses slic for logging and CLI11 to parse command line parameters.
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "../patch/hip_patch.hpp"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <memory>
#include <array>

//-----------------------------------------------------------------------------
/// Basic RAII utility class for initializing and finalizing slic logger
//-----------------------------------------------------------------------------
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

    slic::addStreamToMsgLevel(new slic::GenericOutputStream(&std::cerr), slic::message::Error);
    slic::addStreamToMsgLevel(new slic::GenericOutputStream(&std::cerr, slicFormatStr),
                              slic::message::Warning);

    auto* compactStream = new slic::GenericOutputStream(&std::cout, slicFormatStr);
    slic::addStreamToMsgLevel(compactStream, slic::message::Info);
    slic::addStreamToMsgLevel(compactStream, slic::message::Debug);
  }

  ~BasicLogger() { axom::slic::finalize(); }
};

//-----------------------------------------------------------------------------
/// Struct to help with parsing and storing command line args
//-----------------------------------------------------------------------------
struct Input
{
  std::string mesh_file {""};
  bool verboseOutput {false};

  void parse(int argc, char** argv, axom::CLI::App& app);
  bool isVerbose() const { return verboseOutput; }
};

void Input::parse(int argc, char** argv, axom::CLI::App& app)
{
  app.add_option("-i, --infile", mesh_file)
    ->description("The input STL mesh file")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_flag("-v,--verbose", verboseOutput)
    ->description("Increase logging verbosity?")
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

//-----------------------------------------------------------------------------
/// Basic triangle mesh to be used in our application
//-----------------------------------------------------------------------------
struct TriangleMesh
{
  using Point = axom::primal::Point<double, 3>;
  using Triangle = axom::primal::Triangle<double, 3>;

  axom::IndexType numTriangles() const { return m_triangles.size(); }
  axom::Array<Triangle>& triangles() { return m_triangles; }
  const axom::Array<Triangle>& triangles() const { return m_triangles; }

  axom::Array<Triangle> m_triangles;
};

TriangleMesh makeTriangleMesh(const std::string& stl_mesh_path)
{
  TriangleMesh triMesh;

  // load STL mesh into a mint unstructured mesh
  auto surface_mesh =
    std::make_unique<axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>>(3, axom::mint::TRIANGLE);
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

int main(int argc, char** argv)
{
  // Initialize logger; use RAII so it will finalize at the end of the application
  BasicLogger logger;

  // Use slic macros to log some messages
  SLIC_INFO("A first message!");
  SLIC_WARNING("A warning message!");
  SLIC_DEBUG("A debug message!");

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

  // Load STL mesh into local TriangleMesh struct
  SLIC_INFO(axom::fmt::format("Reading file: '{}'...\n", params.mesh_file));
  TriangleMesh mesh = makeTriangleMesh(params.mesh_file);

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Parsed STL mesh has {:L} triangles.",
                              mesh.numTriangles()));

  return 0;
}
